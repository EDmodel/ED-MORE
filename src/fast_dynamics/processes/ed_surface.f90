! SPDX-License-Identifier: Apache-2.0
!---~---
!   ed_surface -- the ground's surface: the bare-ground skin and the water layer (snow, ice
! and ponded water, one store), blended by the water layer's cover, each giving the soil-top
! boundary and the canopy-air fluxes (docs/science/ground_surface_layers.md):
!   * ground_surface_exchange -- the blend at one canopy-air and soil-top state;
!   * ground_surface_fluxes -- the bare ground's sensible and latent fluxes;
!   * water_layer_stage -- the water layer before a step: the water that reaches it, its
!     free water out to the water network, and a film's temperature; water_layer_cover;
!   * water_layer_commit -- the water layer after a step: its exchange, its free water back
!     and the condensate, and its ice matrix;
!   * snow_blend_ground_optics / ground_tile_radiation -- the ground's optics under the
!     layer, and the radiation each tile absorbs at its own;
!   * surface_stack_balance -- the stack's nodes, the water layer, the litter and the moss,
!     relaxed over the step together, with the heat the water the water network moves through
!     them carries; surface_node_balance -- one node's, its neighbours held.
!---~---
module ed_surface
   use ed_kinds,                only : wp, ik
   use ed_constants,            only : t_3ple, tiny_num, cp_air, cp_ice, cp_liq, cp_vap, stefan_boltzmann, &
                                       latent_heat_fusion, tsupercool_liq, rho_h2o
   use ed_therm_lib,            only : internal_energy_to_temp, sat_specific_humidity,                        &
                                       sat_specific_humidity_temp_deriv, enthalpy_vapor,           &
                                       internal_energy_liquid, precip_enthalpy
   use ed_patch_state_vars,     only : ground_store_t, water_layer_state_t, ground_store_present, store_move,  &
                                       store_take_exchange
   use ed_ground_surface_params, only : water_material_t, N_POROUS
   use ed_linear_algebra,       only : thomas_solve
   use ed_numerics,             only : root_bracket_t, regula_falsi_point, illinois_narrow
   use ed_soil_hydr_lib,        only : soil_layer_curve_t
   use ed_ground_surface_lib,   only : water_layer_cover, water_layer_albedo, water_layer_depth,              &
                                       water_layer_free_water, water_layer_resistance, water_layer_bulk_density, &
                                       matrix_depth_update, porous_energy_state, film_on_porous_layer,          &
                                       water_layer_optics_shares, snow_matrix_conductivity
   implicit none
   private

   public :: ground_surface_exchange, ground_surface_fluxes
   public :: water_layer_stage, water_layer_commit, water_layer_cover_of, water_layer_ice, water_layer_exposed
   public :: snow_blend_ground_optics
   public :: ground_tile_radiation
   public :: surface_env_t, porous_node_env_t, water_stage_out_t, surface_node_t, node_flux_t, store_exchange_t
   public :: STORE_WATER, N_SURFACE_STORE
   public :: surface_node_balance, surface_stack_balance, conductance_in_series
   public :: STACK_NODES_MAX

   !---~---
   !   The stack's solve (surface_stack_balance): the most nodes it holds, the water layer and the
   ! porous layers; Newton's most iterations before the Gauss-Seidel sweeps, and their most sweeps;
   ! and the temperature [K] within which a node's balance is met, the soil heat solve's.
   !---~---
   integer(ik), parameter :: STACK_NODES_MAX  = 3_ik
   integer(ik), parameter :: STACK_NEWTON_MAX = 12_ik
   integer(ik), parameter :: STACK_SWEEP_MAX  = 500_ik
   real(wp),    parameter :: STACK_TEMP_TOL   = 1.0e-8_wp
   !--- The most steps of a porous node's own root on its freezing curve (curve_node_root).
   integer(ik), parameter :: STACK_ROOT_MAX   = 60_ik

   !---~---
   !   The ground surface's stores that exchange with the canopy air apart from the soil: the
   ! water layer (STORE_WATER) and the porous layers, porous layer k's store STORE_WATER + k. Each
   ! one's exchange is one entry of a store_exchange_t, which every record of the march carries:
   ! the energy the store gains and the vapour it gives the canopy air, as rates in an evaluation
   ! and as amounts once the march has summed them. The update commits each store's amounts to it.
   !---~---
   integer(ik), parameter :: STORE_WATER = 1_ik, N_SURFACE_STORE = 1_ik + N_POROUS
   type :: store_exchange_t
      real(wp) :: energy_in(N_SURFACE_STORE)  = 0.0_wp   ! [W/m2] or [J/m2]       energy into the store
      real(wp) :: vapour_out(N_SURFACE_STORE) = 0.0_wp   ! [kg/m2/s] or [kg/m2]   vapour from it to the canopy air
      real(wp) :: buried_heat = 0.0_wp                   ! [W/m2] or [J/m2]       from the water layer into the buried tissue
   end type store_exchange_t

   !---~---
   !   One porous layer as the stack's balance reads it over a step: its dry heat capacity (none
   ! where there is no layer), its depth and its freezing curve (its retention curve at its
   ! depth), on which its energy holds its water's phase; the water it keeps over the step, the
   ! water network's end water, its own evaporation gone; its energy, temperature and liquid
   ! fraction at the step's start; its thermal resistance top to bottom; its evaporation, the
   ! water network's, held over the step per unit ground area; and the water network's net flow
   ! over the step out of its bottom, into the next layer below or the soil (< 0: up into it).
   !---~---
   type :: porous_node_env_t
      real(wp) :: dry_cap    = 0.0_wp                   ! [J/m2/K]
      real(wp) :: depth      = 0.0_wp                   ! [m]
      type(soil_layer_curve_t) :: curve
      real(wp) :: water      = 0.0_wp                   ! [kg/m2]
      real(wp) :: energy     = 0.0_wp                   ! [J/m2]
      real(wp) :: temp       = t_3ple                   ! [K]
      real(wp) :: fliq       = 1.0_wp                   ! [-]
      real(wp) :: resistance = 0.0_wp                   ! [m2 K/W]
      real(wp) :: vapour     = 0.0_wp                   ! [kg/m2/s]
      real(wp) :: down       = 0.0_wp                   ! [kg/m2/s]
   end type porous_node_env_t

   !---~---
   !   The ground surface's surroundings, read by the stack's balance and the bare ground alike:
   ! the radiation each tile absorbs and its conductance to the canopy air, the canopy air's
   ! and the soil top's state, the thermal resistance below the surface stores, the bare
   ! ground's evaporation, and the pack it exchanges with, the pack's cover and parameters, and
   ! the step its balance relaxes over (the ground's stage sets the water layer and its cover).
   !---~---
   type :: surface_env_t
      real(wp) :: rnet_soil = 0.0_wp                    ! [W/m2] the bare soil's net all-wave radiation per unit of its area
      real(wp) :: sw_snow   = 0.0_wp                    ! [W/m2] the shortwave the snow absorbs per unit of its area
      real(wp) :: lw_down   = 0.0_wp                    ! [W/m2] the longwave reaching the ground
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp ! [K],[kg/kg] the canopy air's temperature and specific humidity
      real(wp) :: ggnet    = 0.0_wp                     ! [m/s] ground to canopy-air conductance, heat and vapour alike
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp ! [kg/m3],[Pa] the canopy air's density and pressure
      real(wp) :: cp_cas   = cp_air                     ! [J/kg/K] its specific heat (cas_heat_capacity); dry until set
      real(wp) :: t_soil_top = 273.16_wp                ! [K] top soil node's temperature
      !---~---
      !   The thermal resistance from the ground's surface down to the top soil node, its depth
      ! over the top layer's conductivity, set with the step's record (build_frozen); and the
      ! shortwave the ground the water layer leaves absorbs per unit of its area, the porous
      ! layers' or the bare soil's, whose node forms its longwave at its own temperature.
      !---~---
      real(wp) :: r_soil_top = 0.0_wp                   ! [m2 K/W]
      real(wp) :: sw_ground  = 0.0_wp                   ! [W/m2]
      real(wp) :: soil_evap  = 0.0_wp                   ! [kg/m2/s] the bare ground's evaporation, already area-weighted
      !---~---
      !   The porous layers above the soil, top down, litter over moss, each a node of the stack,
      ! and their emissivity, the ground's.
      !---~---
      type(porous_node_env_t) :: porous(N_POROUS)
      real(wp) :: ground_emiss   = 0.95_wp              ! [-]
      !---~---
      !   The water layer's free water over the step: the water network's net flow from it into
      ! the ground below (< 0: up into it), its runoff, its temperature at the stage, and what it
      ! held there and holds at the step's end. Where the water layer is a layer it is part of the
      ! layer's node, which exchanges; where the store is a film, the free water does not exchange
      ! and the water that leaves it leaves at its stage temperature. The runoff, the free water
      ! the pond's capacity spills, leaves at that temperature either way.
      !---~---
      real(wp) :: free_down = 0.0_wp                    ! [kg/m2/s]
      real(wp) :: runoff    = 0.0_wp                    ! [kg/m2/s]
      real(wp) :: free_temp = t_3ple                    ! [K]
      !---~---
      !   The plant tissue a snow pack buries, which the water layer's node holds at its own
      ! temperature (docs/science/vegetation_energy_dynamics.md, "Buried tissue"): its heat
      ! capacity and its heat at the step's start, from zero kelvin.
      !---~---
      real(wp) :: buried_cap    = 0.0_wp                ! [J/m2/K]
      real(wp) :: buried_energy = 0.0_wp                ! [J/m2]
      type(ground_store_t) :: free_start                ! [kg/m2], [J/m2]
      real(wp) :: free_end  = 0.0_wp                    ! [kg/m2]
      type(water_layer_state_t) :: pack                 ! the water layer over the step, its free water apart
      type(water_layer_state_t) :: exposed              ! the water layer and its free water at the stage
      real(wp) :: snowfac    = 0.0_wp                   ! [-] its cover (water_layer_cover_of; 0 for a film)
      type(water_material_t) :: params                  ! the water layer's material
      real(wp) :: dt         = 0.0_wp                   ! [s] the step
   end type surface_env_t

   !---~---
   !   What the water layer's stage returns, per unit ground area: whether what it holds over the
   ! step is a layer (from min_depth deep), which exchanges, or a film or a skim, which does not;
   ! its free water and that water's energy, which the water network moves over the step, and
   ! its temperature; the enthalpy the throughfall and the shed water brought across the column's
   ! boundary; the heat a film took from the node it lies on; and the store's ice at the stage's
   ! start, the snow it took, and the ice it holds over the step, which the update's matrix law
   ! and the ice-melt export read.
   !---~---
   type :: water_stage_out_t
      logical  :: is_layer    = .false.
      type(ground_store_t) :: free                      ! [kg/m2], [J/m2]
      real(wp) :: free_temp   = t_3ple                  ! [K]
      real(wp) :: inflow_enth = 0.0_wp                  ! [J/m2]
      real(wp) :: film_heat   = 0.0_wp                  ! [J/m2] from the node it lies on, into the film
      real(wp) :: ice_start   = 0.0_wp                  ! [kg/m2]
      real(wp) :: snowfall    = 0.0_wp                  ! [kg/m2]
      real(wp) :: ice_held    = 0.0_wp                  ! [kg/m2]
   end type water_stage_out_t

   !---~---
   !   One node of the ground's surface, the water layer, the litter or the moss, as its balance
   ! over a step reads it: its water over the step (a porous layer's what it keeps, its held
   ! vapour apart: node_linearization) and its energy at the step's start, with its
   ! temperature and liquid fraction there, the point its exchange is linearized about, and the law
   ! its energy holds its water's phase by: pure water's for the water layer, the soil's freezing
   ! curve for a porous layer (docs/science/soil_biophysics.md eqs. (14)-(15)); the share
   ! of the ground it exchanges with the canopy air over, the net radiation it absorbs per unit of
   ! that share at its temperature, its vapour, either the canopy air's (live, from its own
   ! saturation) or the water network's rate held over the step, per unit ground area; its
   ! conductances to the node above and the node below, per unit ground area, with their
   ! temperatures; and the water the water network moves through it over the step, as steady
   ! rates: what arrives from the node above and from the node below, each at its giver's
   ! temperature, the heat the water from a store of held temperature (the free water above, the
   ! soil below) brings it, and what leaves it, at its own temperature.
   !---~---
   type :: surface_node_t
      real(wp) :: dry_cap     = 0.0_wp      ! [J/m2/K] heat capacity of its dry matter (none for water)
      logical  :: on_curve    = .false.     ! its phase on a freezing curve (a porous layer), else pure water's
      real(wp) :: depth       = 0.0_wp      ! [m]      a porous layer's, with its curve
      type(soil_layer_curve_t) :: curve
      real(wp) :: water       = 0.0_wp      ! [kg/m2]  liquid and ice
      real(wp) :: energy      = 0.0_wp      ! [J/m2]   internal energy
      real(wp) :: temp        = t_3ple      ! [K]
      real(wp) :: fliq        = 0.0_wp      ! [-]
      real(wp) :: area        = 0.0_wp      ! [-]      the share of the ground it exchanges over
      real(wp) :: rad         = 0.0_wp      ! [W/m2]   net all-wave radiation per unit of that share
      real(wp) :: emiss       = 1.0_wp      ! [-]      its emissivity
      logical  :: vapour_held = .false.     ! its vapour is the water network's
      real(wp) :: vapour_rate = 0.0_wp      ! [kg/m2/s] that rate (the network gives a porous layer no dew)
      real(wp) :: g_above     = 0.0_wp      ! [W/m2/K] to the node above
      real(wp) :: t_above     = t_3ple      ! [K]      the node above's temperature
      real(wp) :: g_below     = 0.0_wp      ! [W/m2/K] to the node below
      real(wp) :: t_below     = t_3ple      ! [K]      the node below's temperature
      real(wp) :: m_above     = 0.0_wp      ! [kg/m2/s] water from the node above
      real(wp) :: m_below     = 0.0_wp      ! [kg/m2/s] water from the node below
      real(wp) :: q_in        = 0.0_wp      ! [W/m2]   the heat water from a store of held temperature brings
      real(wp) :: m_out       = 0.0_wp      ! [kg/m2/s] water leaving it
   end type surface_node_t

   !---~---
   !   What a node's balance returns, per unit ground area, positive away from the node: its net
   ! radiation, its sensible and latent heat and vapour to the canopy air, its conduction to the
   ! node below, and its net gain; and the temperature it ends at.
   !---~---
   type :: node_flux_t
      real(wp) :: t_surf    = t_3ple                    ! [K]
      real(wp) :: rad       = 0.0_wp                    ! [W/m2]
      real(wp) :: h         = 0.0_wp                    ! [W/m2]
      real(wp) :: le        = 0.0_wp                    ! [W/m2]
      real(wp) :: w_flux    = 0.0_wp                    ! [kg/m2/s]
      real(wp) :: g_below   = 0.0_wp                    ! [W/m2]
      real(wp) :: energy_in = 0.0_wp                    ! [W/m2]
   end type node_flux_t

   !---~---
   !   A node's exchange with the canopy air and its emission, linearized about its temperature at
   ! the step's start: the sensible heat and its slope, its vapour live (and its slope) or held,
   ! the vapour's latent heat less the water's own energy, the energy a kilogram of its held vapour
   ! takes from its balance and that energy's slope, and its emission's slope.
   !---~---
   type :: node_lin_t
      real(wp) :: h_n = 0.0_wp, dh_dt = 0.0_wp           ! [W/m2], [W/m2/K]
      real(wp) :: w_n = 0.0_wp, dw_dt = 0.0_wp           ! [kg/m2/s], [kg/m2/s/K] per unit of the node's share
      real(wp) :: vh  = 0.0_wp                           ! [kg/m2/s] held, per unit ground area
      real(wp) :: l_vap = 0.0_wp                         ! [J/kg]
      real(wp) :: l_held = 0.0_wp, dl_held = 0.0_wp      ! [J/kg], [J/kg/K]
      real(wp) :: demit_dt = 0.0_wp                      ! [W/m2/K]
   end type node_lin_t

contains

   !--- Ground / surface skin

   !---~---
   !   ground_surface_exchange -- the ground surface at one canopy-air and soil-top state
   ! (docs/science/ground_surface_layers.md §3): its layers top down, the water layer over its
   ! cover, the litter and the moss, each a thermal node, relaxed over the step together
   ! (surface_stack_balance) over the soil top; and the bare soil over what the water layer leaves
   ! where no porous layer lies. The top porous
   ! layer exchanges over what the water layer leaves, the one beneath it with its neighbours
   ! only. A node's temperature is its top face's, the one it emits and exchanges at, so the heat
   ! it passes down crosses its whole depth to the next node's top face, and the bottom node's the
   ! top soil node's upper part (r_soil_top) too; the water layer passes it over its cover only.
   ! Each node forms its longwave at its own temperature. The water the water network moved
   ! through the stack over the step carries its heat at its giver's temperature, a node's at its
   ! own over the step: from the free water, which keeps its stage temperature, down through the
   ! porous layers into the soil, or up from the soil at its top's temperature (advect_through).
   ! The fluxes are zero for a layer that is absent, so with no layer this is exactly bare
   ! ground. soil_evap already carries its share of the ground (ground_evap_env_t%snow_free_frac),
   ! so le_soil is added whole; the soil top pays for it.
   !---~---
   pure subroutine ground_surface_exchange(env, h_ground, le_ground, g_top, ground_rad, stores, iterations, fallback)
      type(surface_env_t),    intent(in)  :: env          !        the ground surface's surroundings
      real(wp),               intent(out) :: h_ground     ! [W/m2] sensible heat to the canopy air
      real(wp),               intent(out) :: le_ground    ! [W/m2] latent heat to the canopy air
      real(wp),               intent(out) :: g_top        ! [W/m2] net energy into the soil top, its water's heat included
      real(wp),               intent(out) :: ground_rad   ! [W/m2] net radiation, the layers' and the bare soil's
      type(store_exchange_t), intent(out) :: stores       !        each surface store's exchange (0 without it)
      integer(ik),            intent(out) :: iterations   !        the stack's Newton iterations
      integer(ik),            intent(out) :: fallback     !        1 if the stack fell back to Gauss-Seidel sweeps
      type(surface_node_t) :: node(STACK_NODES_MAX)
      type(node_flux_t)    :: flux(STACK_NODES_MAX)
      real(wp)    :: resistance(STACK_NODES_MAX), contact(STACK_NODES_MAX), cover, h_bare, le_soil
      integer(ik) :: store(STACK_NODES_MAX), n, i, k, top_porous
      logical     :: water_node   ! the water layer is the stack's top node
      type(water_layer_state_t) :: exposed
      stores = store_exchange_t()
      cover  = max(0.0_wp, min(1.0_wp, env%snowfac))
      n = 0_ik
      water_node = .false.
      if (env%pack%water + env%free_start%water > 0.0_wp .and. cover > 0.0_wp) then
         !---~---
         !   The water layer's node: what it held over the step and its free water, at their start
         ! energy and temperature (water_layer_exposed), holding at the step's end what the water
         ! network leaves of the free water, and the plant tissue the pack buries at its
         ! temperature (vegetation_energy_dynamics.md); the water it gives the ground leaves at the
         ! node's temperature (advect_through).
         !---~---
         exposed = water_layer_exposed(env%pack, env%free_start)
         n = n + 1_ik
         water_node        = .true.
         store(n)          = STORE_WATER
         node(n)%water     = env%pack%water + env%free_end
         node(n)%energy    = exposed%energy + env%buried_energy
         node(n)%dry_cap   = env%buried_cap
         node(n)%temp      = exposed%temp
         node(n)%fliq      = exposed%fliq
         node(n)%area      = cover
         node(n)%emiss     = env%params%emissivity
         node(n)%rad       = env%sw_snow
         resistance(n)     = water_layer_resistance(layer_depth(exposed), exposed%matrix_depth,                  &
                                                    matrix_conductivity(exposed, env%params))
         contact(n)        = cover
      end if
      top_porous = n + 1_ik
      do k = 1_ik, N_POROUS
         if (.not. (env%porous(k)%dry_cap > 0.0_wp)) cycle
         n = n + 1_ik
         store(n)            = STORE_WATER + k
         node(n)%dry_cap     = env%porous(k)%dry_cap
         node(n)%on_curve    = .true.
         node(n)%depth       = env%porous(k)%depth
         node(n)%curve       = env%porous(k)%curve
         node(n)%water       = env%porous(k)%water
         node(n)%energy      = env%porous(k)%energy
         node(n)%temp        = env%porous(k)%temp
         node(n)%fliq        = env%porous(k)%fliq
         if (n == top_porous) node(n)%area = 1.0_wp - cover
         node(n)%emiss       = env%ground_emiss
         node(n)%rad         = env%sw_ground
         node(n)%vapour_held = .true.
         node(n)%vapour_rate = env%porous(k)%vapour
         resistance(n)       = env%porous(k)%resistance
         contact(n)          = 1.0_wp
      end do
      do i = 1_ik, n
         node(i)%rad = node(i)%rad + node(i)%emiss * (env%lw_down - stefan_boltzmann * node(i)%temp ** 4)
         if (i < n) then
            node(i)%g_below = contact(i) / max(resistance(i), tiny_num)
         else
            node(i)%g_below = contact(i) * conductance_in_series(resistance(i), env%r_soil_top)
            node(i)%t_below = env%t_soil_top
         end if
         if (i > 1_ik) node(i)%g_above = node(i - 1_ik)%g_below
      end do
      call advect_through(env, water_node, top_porous, n, node)
      call surface_stack_balance(node, n, env, flux, iterations, fallback)
      call ground_surface_fluxes(env%t_soil_top, env%can_temp, env%ggnet, env%rho_air, env%cp_cas, env%soil_evap,    &
                                 h_bare, le_soil)
      h_ground   = 0.0_wp
      le_ground  = le_soil
      ground_rad = 0.0_wp
      g_top      = -le_soil
      do i = 1_ik, n
         h_ground   = h_ground   + flux(i)%h
         le_ground  = le_ground  + flux(i)%le
         ground_rad = ground_rad + flux(i)%rad
         stores%energy_in(store(i))  = flux(i)%energy_in
         stores%vapour_out(store(i)) = flux(i)%w_flux
      end do
      if (n > 0_ik) g_top = g_top + flux(n)%g_below
      if (top_porous > n) then
         h_ground   = h_ground   + (1.0_wp - cover) * h_bare
         ground_rad = ground_rad + (1.0_wp - cover) * env%rnet_soil
         g_top      = g_top      + (1.0_wp - cover) * (env%rnet_soil - h_bare)
      end if
      call advected_heat(env, water_node, top_porous, n, flux, g_top, stores%energy_in(STORE_WATER))
      !--- The buried tissue's share of the water layer's node's heat, a matched pair with the layer.
      if (water_node .and. env%buried_cap > 0.0_wp) then
         stores%buried_heat = (env%buried_cap * flux(1)%t_surf - env%buried_energy) / env%dt
         stores%energy_in(STORE_WATER) = stores%energy_in(STORE_WATER) - stores%buried_heat
      end if
   end subroutine ground_surface_exchange

   !---~---
   !   The water the water network moved through the stack over the step, as the nodes read it
   ! (surface_node_t): down a chain from the free water through the porous layers (nodes
   ! top_porous .. n) into the soil, each face's net flow. Water a face carries from a node leaves
   ! it at its own temperature (m_out) and reaches the node across at that temperature (m_above,
   ! m_below); water from the soil, whose temperature is held, brings the node a fixed heat (q_in).
   ! Where the water layer is the stack's top node (water_node), the free water is part of it: what
   ! it gives the ground leaves it at its temperature, and what rises into it arrives at its giver's;
   ! its runoff leaves at the free water's stage temperature, a fixed heat. Where it is not, the
   ! free water keeps its stage temperature, and the water it gives the top porous layer brings
   ! that layer a fixed heat.
   !---~---
   pure subroutine advect_through(env, water_node, top_porous, n, node)
      type(surface_env_t),  intent(in)    :: env
      logical,              intent(in)    :: water_node      ! node 1 is the water layer's
      integer(ik),          intent(in)    :: top_porous, n   ! the porous layers are nodes top_porous .. n
      type(surface_node_t), intent(inout) :: node(:)
      real(wp)    :: q
      integer(ik) :: i, k
      q = env%free_down                       ! down the face above the top porous layer
      if (water_node) then
         node(1)%m_out = max(q, 0.0_wp)
         node(1)%q_in  = -env%runoff * internal_energy_liquid(env%free_temp)
         if (n < top_porous) node(1)%q_in = node(1)%q_in + max(-q, 0.0_wp) * internal_energy_liquid(env%t_soil_top)
      end if
      i = top_porous - 1_ik
      do k = 1_ik, N_POROUS
         if (.not. (env%porous(k)%dry_cap > 0.0_wp)) cycle
         i = i + 1_ik
         if (i == top_porous .and. .not. water_node) then
            node(i)%q_in  = node(i)%q_in + max(q, 0.0_wp) * internal_energy_liquid(env%free_temp)
            node(i)%m_out = node(i)%m_out + max(-q, 0.0_wp)
         else
            node(i)%m_above      = max(q, 0.0_wp)
            node(i - 1_ik)%m_below = max(-q, 0.0_wp)
            node(i)%m_out        = node(i)%m_out + max(-q, 0.0_wp)
            if (i > top_porous) node(i - 1_ik)%m_out = node(i - 1_ik)%m_out + max(q, 0.0_wp)
         end if
         q = env%porous(k)%down
      end do
      !--- The bottom porous layer's face into the soil: water up from the soil at its top's temperature.
      if (n >= top_porous) then
         node(n)%m_out = node(n)%m_out + max(q, 0.0_wp)
         node(n)%q_in  = node(n)%q_in + max(-q, 0.0_wp) * internal_energy_liquid(env%t_soil_top)
      end if
   end subroutine advect_through

   !---~---
   !   The heat the water the water network moved carried at the stack's solution (advect_through's
   ! chain): into the soil top, the bottom porous layer's down at its temperature, or the free
   ! water's, at the water layer's node's or its stage temperature, where no porous layer lies,
   ! less what rose from the soil at its top's; and, where the water layer is no node, the free
   ! water's own change, the water and its runoff leaving it at its stage temperature and the
   ! water rising into it at its giver's, which the water layer's store takes (a node's own
   ! balance holds them).
   !---~---
   pure subroutine advected_heat(env, water_node, top_porous, n, flux, g_top, free_heat)
      type(surface_env_t), intent(in)    :: env
      logical,             intent(in)    :: water_node
      integer(ik),         intent(in)    :: top_porous, n
      type(node_flux_t),   intent(in)    :: flux(:)
      real(wp),            intent(inout) :: g_top       ! [W/m2] into the soil top
      real(wp),            intent(inout) :: free_heat   ! [W/m2] into the water layer's store
      real(wp)    :: q_soil, t_giver, t_first
      integer(ik) :: k
      q_soil  = env%free_down
      t_giver = env%free_temp
      if (water_node) t_giver = flux(1)%t_surf
      do k = 1_ik, N_POROUS
         if (env%porous(k)%dry_cap > 0.0_wp) q_soil = env%porous(k)%down
      end do
      t_first = env%t_soil_top
      if (n >= top_porous) then
         t_giver = flux(n)%t_surf
         t_first = flux(top_porous)%t_surf
      end if
      g_top = g_top + max(q_soil, 0.0_wp) * internal_energy_liquid(t_giver)                                     &
                    - max(-q_soil, 0.0_wp) * internal_energy_liquid(env%t_soil_top)
      if (water_node) return
      free_heat = free_heat - (max(env%free_down, 0.0_wp) + env%runoff) * internal_energy_liquid(env%free_temp)     &
                            + max(-env%free_down, 0.0_wp) * internal_energy_liquid(t_first)
   end subroutine advected_heat

   !---~---
   !   Bare-ground fluxes to the canopy air, instantaneous: sensible through ggnet at the
   ! canopy air's density and specific heat, latent from the soil evaporation (already
   ! area-weighted) at the vapour enthalpy the canopy air receives. ground_surface_exchange
   ! assembles G_top and the snow blend.
   !---~---
   elemental pure subroutine ground_surface_fluxes(t_ground, t_cas, ggnet, rho, cp_cas, soil_evap, h_bare, le_soil)
      real(wp), intent(in)  :: t_ground, t_cas, ggnet       ! [K], [K], [m/s]
      real(wp), intent(in)  :: rho, cp_cas, soil_evap       ! [kg/m3], [J/kg/K], [kg/m2/s]
      real(wp), intent(out) :: h_bare, le_soil              ! [W/m2], [W/m2]
      h_bare  = sensible_heat_slope(ggnet, rho, cp_cas) * (t_ground - t_cas)
      le_soil = soil_evap * enthalpy_vapor(t_ground)
   end subroutine ground_surface_fluxes

   !---~---
   !   The sensible heat a ground surface store gives the canopy air per kelvin it is warmer
   ! [W/m2/K]: its conductance to the canopy air times the air's heat capacity per volume. The
   ! bare soil, the water layer and the porous layers all exchange through it.
   !---~---
   elemental real(wp) function sensible_heat_slope(ggnet, rho, cp_cas) result(slope)
      real(wp), intent(in) :: ggnet, rho, cp_cas   ! [m/s], [kg/m3], [J/kg/K]
      slope = ggnet * rho * cp_cas
   end function sensible_heat_slope

   !--- [W/m2/K] How a surface's thermal emission grows with its temperature: d(eps sigma T^4)/dT.
   elemental real(wp) function emission_slope(emiss, temp) result(slope)
      real(wp), intent(in) :: emiss, temp   ! [-], [K]
      slope = 4.0_wp * emiss * stefan_boltzmann * temp ** 3
   end function emission_slope

   !--- [W/m2/K] The conductance of two thermal resistances [m2 K/W] in series.
   elemental real(wp) function conductance_in_series(r_upper, r_lower) result(g)
      real(wp), intent(in) :: r_upper, r_lower
      g = 1.0_wp / max(r_upper + r_lower, tiny_num)
   end function conductance_in_series

   !---~---
   !   The ground's optics under a pack, for the two-stream: its albedo and emissivity blended
   ! toward the pack's by its cover, the one the pack's balance exchanges over, so a partial pack
   ! is continuous, and the temperature at which the blended emissivity emits what the soil and
   ! the pack emit at their own (docs/science/ground_surface_layers.md, "The ground's radiation").
   ! Unchanged without cover.
   !---~---
   pure subroutine snow_blend_ground_optics(snow, fc, params, albedo_vis, albedo_nir, emiss, temp)
      type(water_layer_state_t), intent(in) :: snow
      real(wp),            intent(in)    :: fc           ! [-] the water layer's cover (water_layer_cover_of)
      type(water_material_t), intent(in)    :: params
      real(wp),            intent(inout) :: albedo_vis   ! [-] VIS albedo: the bare soil's on entry
      real(wp),            intent(inout) :: albedo_nir   ! [-] NIR albedo: the bare soil's on entry
      real(wp),            intent(inout) :: emiss        ! [-] thermal emissivity: the bare soil's on entry
      real(wp),            intent(inout) :: temp         ! [K] emission temperature: the top soil node's on entry
      real(wp) :: snow_vis, snow_nir, emit_soil, emit_snow
      if (.not. (fc > 0.0_wp)) return
      call snow_albedos(snow, params, snow_vis, snow_nir)
      emit_soil  = emiss * temp ** 4                                ! [K4] emission over sigma
      emit_snow  = params%emissivity * snow%temp ** 4
      albedo_vis = (1.0_wp - fc) * albedo_vis + fc * snow_vis
      albedo_nir = (1.0_wp - fc) * albedo_nir + fc * snow_nir
      emiss      = (1.0_wp - fc) * emiss + fc * params%emissivity
      temp       = (((1.0_wp - fc) * emit_soil + fc * emit_snow) / emiss) ** 0.25_wp
   end subroutine snow_blend_ground_optics

   !---~---
   !   The water layer's own shortwave albedos (water_layer_albedo): the water standing above its
   ! ice matrix at standing water's, the matrix at snow's, fresh to aged by the layer's liquid
   ! fraction, as far as its density lies from solid ice's, else at ice's
   ! (water_layer_optics_shares). A layer this step's snowfall has yet to make is fresh snow.
   !---~---
   pure subroutine snow_albedos(snow, params, albedo_vis, albedo_nir)
      type(water_layer_state_t), intent(in) :: snow
      type(water_material_t), intent(in)  :: params
      real(wp),            intent(out) :: albedo_vis, albedo_nir   ! [-]
      real(wp) :: fl, standing, porous_share
      fl           = 0.0_wp
      standing     = 0.0_wp
      porous_share = 1.0_wp
      if (ground_store_present(snow)) then
         fl = snow%fliq
         call water_layer_optics_shares(snow%matrix_depth, water_layer_ice(snow), snow%fliq * snow%water,            &
                                        params%fresh_density, standing, porous_share)
      end if
      albedo_vis = water_layer_albedo(standing, porous_share, fl, params%albedo_vis_fresh, params%albedo_vis_aged,   &
                                      params%albedo_vis_ice, params%albedo_vis_liquid)
      albedo_nir = water_layer_albedo(standing, porous_share, fl, params%albedo_nir_fresh, params%albedo_nir_aged,   &
                                      params%albedo_nir_ice, params%albedo_nir_liquid)
   end subroutine snow_albedos

   !---~---
   !   [W/m/K] The water layer's ice matrix's conductivity at the matrix's density, its ice over its
   ! depth (snow_matrix_conductivity); the fresh snow's where it has no matrix, whose depth, nil,
   ! then carries no resistance.
   !---~---
   pure real(wp) function matrix_conductivity(layer, params) result(k)
      type(water_layer_state_t), intent(in) :: layer
      type(water_material_t),    intent(in) :: params
      k = params%conductivity_at_fresh
      if (layer%matrix_depth > 0.0_wp)                                                                          &
         k = snow_matrix_conductivity(water_layer_ice(layer) / layer%matrix_depth, params%conductivity_at_fresh,  &
                                      params%fresh_density)
   end function matrix_conductivity

   !---~---
   !   The radiation each tile of the ground absorbs per unit of its own area, from what
   ! reaches the ground in each band (docs/science/ground_surface_layers.md, "The ground's
   ! radiation"): the bare soil's net all-wave radiation at its albedo and emissivity and the
   ! top soil node's temperature, and the shortwave the snow absorbs at its own albedo. The
   ! snow's longwave is formed in its balance, at its own temperature. Weighted by the cover,
   ! the two give the ground absorption of the two-stream, whose optics
   ! snow_blend_ground_optics blends.
   !---~---
   pure subroutine ground_tile_radiation(dn_vis, dn_nir, dn_lw, soil_albedo_vis, soil_albedo_nir, soil_emiss,  &
                                         soil_temp, snow, params, rnet_soil, sw_snow, sw_ground)
      real(wp),            intent(in)  :: dn_vis, dn_nir, dn_lw   ! [W/m2] each band reaching the ground
      real(wp),            intent(in)  :: soil_albedo_vis, soil_albedo_nir, soil_emiss   ! [-] the bare soil's optics
      real(wp),            intent(in)  :: soil_temp               ! [K] the top soil node's temperature
      type(water_layer_state_t), intent(in) :: snow
      type(water_material_t), intent(in)  :: params
      real(wp),            intent(out) :: rnet_soil, sw_snow      ! [W/m2]
      real(wp),            intent(out) :: sw_ground               ! [W/m2] the bare ground's shortwave alone
      real(wp) :: snow_vis, snow_nir
      sw_ground = dn_vis * (1.0_wp - soil_albedo_vis) + dn_nir * (1.0_wp - soil_albedo_nir)
      rnet_soil = sw_ground + soil_emiss * (dn_lw - stefan_boltzmann * soil_temp ** 4)
      call snow_albedos(snow, params, snow_vis, snow_nir)
      sw_snow = dn_vis * (1.0_wp - snow_vis) + dn_nir * (1.0_wp - snow_nir)
   end subroutine ground_tile_radiation

   !---~---
   !   surface_node_balance -- one node's exchange over dt, per unit ground area, its neighbours held,
   ! from a linearized backward-Euler balance of its energy that is not committed
   ! (docs/science/ground_surface_layers.md §3): the balance's one root T, in closed form on the
   ! plateau at the triple point or else on the ice or liquid line, with its dry heat capacity on
   ! each, and the fluxes at T (node_fluxes). Over its share of the ground it takes its net radiation
   ! and gives sensible heat and vapour, a kilogram of which carries enthalpy_vapor less the water's
   ! own energy; it takes heat from the node above and gives it to the node below, implicit in its
   ! own temperature, so a node of little heat capacity relaxes without overshooting; and the water
   ! that passes through it brings its giver's heat and leaves with its own, so the water leaving
   ! over the step is part of what the balance warms or cools (node_storage). Its own vapour, once
   ! the linearization takes it past the water it holds, is held there and the balance solved again,
   ! so a node that loses all its water to the air ends with none of its energy. The stack's
   ! Gauss-Seidel sweeps solve each node by it.
   !---~---
   pure subroutine surface_node_balance(node, env, flux)
      type(surface_node_t), intent(in)  :: node
      type(surface_env_t),  intent(in)  :: env      ! the canopy air, the conductance to it and the step
      type(node_flux_t),    intent(out) :: flux
      type(node_lin_t) :: lin
      real(wp) :: a, dt, w, cd, e0, r0, drdt, u_freeze, u_melt, e_bal, t_new, fl_new, w_cap
      integer  :: pass
      a  = node%area
      dt = env%dt
      w  = node%water
      call node_storage(node, dt, cd, e0)
      call node_linearization(node, env, lin)
      w_cap = huge(1.0_wp)
      if (a > 0.0_wp .and. .not. node%vapour_held) w_cap = w / (dt * a)   ! the vapour that empties it over dt
      !---~---
      !   Solve e(T) = E^n + dt (r0 + drdt (T - T^n)): a porous node on its freezing curve by its own
      ! root (curve_node_root); the water layer on the plateau first, else on the ice or the liquid
      ! line, in closed form.
      !---~---
      do pass = 1, 2
         r0   = node_exchange(node, lin, node%temp) + link_heat(node, node%temp, node%t_above, node%t_below)
         drdt = -(a * (lin%demit_dt + lin%dh_dt + lin%l_vap * lin%dw_dt) + lin%dl_held * lin%vh             &
                  + node%g_above + node%g_below)                                                                ! < 0
         if (node%on_curve) then
            call curve_node_root(node, cd, e0, r0, drdt, dt, t_new, fl_new)
         else
            u_freeze = (cd + w * cp_ice) * t_3ple
            u_melt   = u_freeze + w * latent_heat_fusion
            t_new    = t_3ple
            e_bal    = e0 + dt * (r0 + drdt * (t_3ple - node%temp))
            fl_new   = min(1.0_wp, max(0.0_wp, (e_bal - u_freeze) / max(w * latent_heat_fusion, tiny_num)))
            if (e_bal < u_freeze) then
               t_new  = (e0 + dt * (r0 - drdt * node%temp)) / (cd + w * cp_ice - dt * drdt)
               fl_new = 0.0_wp
            else if (e_bal > u_melt) then
               t_new  = (e0 + w * cp_liq * tsupercool_liq + dt * (r0 - drdt * node%temp))                      &
                        / (cd + w * cp_liq - dt * drdt)
               fl_new = 1.0_wp
            end if
         end if
         if (lin%w_n + lin%dw_dt * (t_new - node%temp) <= w_cap) exit
         lin%w_n   = w_cap                                      ! held, and solved again
         lin%dw_dt = 0.0_wp
      end do
      call node_fluxes(node, lin, t_new, fl_new, flux)
   end subroutine surface_node_balance

   !---~---
   !   surface_stack_balance -- the stack's nodes, top down, relaxed over the step together: each
   ! node's linearized backward-Euler balance (surface_node_balance's) with the heat it takes
   ! from the node above and gives the node below, by conduction and with the water passing
   ! between them, implicit in both their temperatures, the bottom node's to the soil top at its
   ! temperature held. Newton on the nodes' energies, each node's temperature read off its energy
   ! on its own phase law (the ice line, the plateau at the triple point, the liquid line), so a node crossing a
   ! kink takes the next step from where it landed; the tridiagonal system's columns are
   ! diagonally dominant (the water a node gives is part of its own heat capacity), so it is never
   ! singular. A node whose live vapour would take more than it holds has that vapour held at
   ! what empties it, and the stack is solved again. Newton not converged within
   ! STACK_NEWTON_MAX, the nodes are swept top down, each its own balance with its neighbours held
   ! (surface_node_balance), until none moves by STACK_TEMP_TOL: nonlinear Gauss-Seidel, which
   ! converges here because each node's energy rises with its own temperature and falls with its
   ! neighbours'. The fluxes are formed once, at the solution's temperatures, so what one node
   ! gives the next takes. Returns the Newton iterations and whether it fell back (fallback = 1).
   !---~---
   pure subroutine surface_stack_balance(node, n, env, flux, iterations, fallback, newton_max)
      type(surface_node_t), intent(inout) :: node(:)       ! held vapour may be set
      integer(ik),          intent(in)    :: n
      type(surface_env_t),  intent(in)    :: env
      type(node_flux_t),    intent(out)   :: flux(:)
      integer(ik),          intent(out)   :: iterations
      integer(ik),          intent(out)   :: fallback
      integer(ik), optional, intent(in)   :: newton_max    ! Newton's most iterations (STACK_NEWTON_MAX)
      type(node_lin_t)  :: lin(n)
      type(node_flux_t) :: one
      real(wp)    :: e(n), e0(n), cd(n), t(n), fl(n), dtde(n), res(n), cap(n), slope(n)
      real(wp)    :: lower(n), diag(n), upper(n), de(n)
      real(wp)    :: moved
      integer(ik) :: i, it, sweep, it_max
      integer     :: pass
      logical     :: converged, capped
      iterations = 0_ik
      fallback   = 0_ik
      if (n == 0_ik) return
      it_max = STACK_NEWTON_MAX
      if (present(newton_max)) it_max = newton_max
      do i = 1_ik, n
         call node_storage(node(i), env%dt, cd(i), e0(i))
      end do
      do pass = 1, 2
         do i = 1_ik, n
            call node_linearization(node(i), env, lin(i))
            slope(i) = -(node(i)%area * (lin(i)%demit_dt + lin(i)%dh_dt + lin(i)%l_vap * lin(i)%dw_dt)         &
                         + lin(i)%dl_held * lin(i)%vh + node(i)%g_above + node(i)%g_below)
         end do
         !---~---
         !   Newton starts each node at its start temperature: a porous node's energy at the step's
         ! start holds the water its vapour takes, so it starts less that water's own energy.
         !---~---
         e = e0
         do i = 1_ik, n
            if (node(i)%on_curve .and. node(i)%vapour_held) e(i) = e0(i) - env%dt * node(i)%vapour_rate          &
                                                                   * water_specific_energy(node(i)%temp, node(i)%fliq)
         end do
         converged = .false.
         do it = 1_ik, it_max + 1_ik
            call stack_residual(node, lin, n, env%dt, e0, cd, e, t, fl, dtde, res, cap)
            if (all(abs(res) <= STACK_TEMP_TOL * (cap + env%dt * abs(slope)))) then
               converged = .true.
               exit
            end if
            if (it > it_max) exit
            iterations = iterations + 1_ik
            do i = 1_ik, n
               diag(i)  = 1.0_wp - env%dt * slope(i) * dtde(i)
               lower(i) = 0.0_wp
               upper(i) = 0.0_wp
               if (i > 1_ik) lower(i) = -env%dt * (node(i)%g_above + cp_liq * node(i)%m_above) * dtde(i - 1_ik)
               if (i < n)    upper(i) = -env%dt * (node(i)%g_below + cp_liq * node(i)%m_below) * dtde(i + 1_ik)
            end do
            call thomas_solve(lower, diag, upper, -res, de, n)
            e = e + de
         end do
         if (.not. converged) then
            fallback = 1_ik
            do sweep = 1_ik, STACK_SWEEP_MAX
               moved = 0.0_wp
               do i = 1_ik, n
                  if (i > 1_ik) node(i)%t_above = t(i - 1_ik)
                  if (i < n)    node(i)%t_below = t(i + 1_ik)
                  call surface_node_balance(node(i), env, one)
                  moved = max(moved, abs(one%t_surf - t(i)))
                  t(i)  = one%t_surf
                  e(i)  = e0(i) + env%dt * (node_exchange(node(i), lin(i), t(i))                                &
                                            + link_heat(node(i), t(i), neighbour_temp(i - 1_ik), neighbour_temp(i + 1_ik)))
               end do
               if (moved <= STACK_TEMP_TOL) exit
            end do
            call stack_residual(node, lin, n, env%dt, e0, cd, e, t, fl, dtde, res, cap)
         end if
         !--- A live vapour past what its node holds is held at what empties it, and the stack solved again.
         capped = .false.
         do i = 1_ik, n
            if (node(i)%vapour_held .or. .not. (node(i)%area > 0.0_wp)) cycle
            if (lin(i)%w_n + lin(i)%dw_dt * (t(i) - node(i)%temp) > node(i)%water / (env%dt * node(i)%area)) then
               node(i)%vapour_held = .true.
               node(i)%vapour_rate = node(i)%water / env%dt
               capped = .true.
            end if
         end do
         if (.not. capped) exit
      end do
      !--- The fluxes at the solution, each link's once.
      do i = 1_ik, n
         if (i > 1_ik) node(i)%t_above = t(i - 1_ik)
         if (i < n)    node(i)%t_below = t(i + 1_ik)
         call node_fluxes(node(i), lin(i), t(i), fl(i), flux(i))
      end do
   contains
      !--- [K] Node j's temperature in the sweep, its neighbour's held one beyond the stack's ends.
      pure real(wp) function neighbour_temp(j) result(tj)
         integer(ik), intent(in) :: j
         if (j < 1_ik) then
            tj = node(1)%t_above
         else if (j > n) then
            tj = node(n)%t_below
         else
            tj = t(j)
         end if
      end function neighbour_temp
   end subroutine surface_stack_balance

   !---~---
   !   The stack's residuals at the nodes' energies e: each node's temperature, liquid fraction and
   ! dT/dE on its phase law, its heat capacity there (node_phase, with the water it gives as part
   ! of its heat capacity, node_storage), and R = e - E^n - dt Phi, Phi its exchange
   ! (node_exchange) plus its links at those temperatures (link_heat), the top node's above and
   ! the bottom node's below at the held ones.
   !---~---
   pure subroutine stack_residual(node, lin, n, dt, e0, cd, e, t, fl, dtde, res, cap)
      type(surface_node_t), intent(in)  :: node(:)
      type(node_lin_t),     intent(in)  :: lin(:)
      integer(ik),          intent(in)  :: n
      real(wp),             intent(in)  :: dt, e0(:), cd(:), e(:)
      real(wp),             intent(out) :: t(:), fl(:), dtde(:), res(:), cap(:)
      real(wp)    :: t_up, t_down
      integer(ik) :: i
      do i = 1_ik, n
         call node_phase(node(i), e(i), cd(i), t(i), fl(i), dtde(i), cap(i))
      end do
      do i = 1_ik, n
         t_up   = node(i)%t_above
         t_down = node(i)%t_below
         if (i > 1_ik) t_up   = t(i - 1_ik)
         if (i < n)    t_down = t(i + 1_ik)
         res(i) = e(i) - e0(i) - dt * (node_exchange(node(i), lin(i), t(i)) + link_heat(node(i), t(i), t_up, t_down))
      end do
   end subroutine stack_residual

   !---~---
   !   A node's storage over the step [J/m2/K], [J/m2]: the water leaving it over the step goes at
   ! its own temperature, as liquid, so its heat over the step is part of the node's, a heat
   ! capacity m_out dt cp_liq beside its dry matter's (cd), and its energy at the step's start
   ! (e0) counts that water's energy from the liquid line's datum. The phase law then gives the
   ! node's temperature from e0 and what it takes over the step.
   !---~---
   pure subroutine node_storage(node, dt, cd, e0)
      type(surface_node_t), intent(in)  :: node
      real(wp),             intent(in)  :: dt
      real(wp),             intent(out) :: cd, e0
      real(wp) :: c_out
      c_out = node%m_out * dt * cp_liq
      cd    = node%dry_cap + c_out
      e0    = node%energy + c_out * tsupercool_liq
   end subroutine node_storage

   !---~---
   !   [W/m2] The heat a node takes over its links at temperature t, its neighbours at t_up and
   ! t_down: conduction from the node above and to the node below, and the water arriving from
   ! each at its giver's temperature and from the stores of held temperature (q_in).
   !---~---
   pure real(wp) function link_heat(node, t, t_up, t_down) result(q)
      type(surface_node_t), intent(in) :: node
      real(wp),             intent(in) :: t, t_up, t_down
      q = node%g_above * (t_up - t) - node%g_below * (t - t_down)                                               &
          + node%m_above * internal_energy_liquid(t_up) + node%m_below * internal_energy_liquid(t_down) + node%q_in
   end function link_heat

   !---~---
   !   [W/m2] A node's exchange with the canopy air and its radiation at temperature t, per unit
   ! ground area, on its linearization: over its share of the ground its net radiation less its
   ! sensible heat and its vapour's latent heat, and a held vapour's latent heat.
   !---~---
   pure real(wp) function node_exchange(node, lin, t) result(phi)
      type(surface_node_t), intent(in) :: node
      type(node_lin_t),     intent(in) :: lin
      real(wp),             intent(in) :: t
      real(wp) :: dtemp
      dtemp = t - node%temp
      phi = node%area * (node%rad - lin%demit_dt * dtemp - lin%h_n - lin%dh_dt * dtemp                          &
                         - lin%l_vap * (lin%w_n + lin%dw_dt * dtemp))                                           &
            - (lin%l_held + lin%dl_held * dtemp) * lin%vh
   end function node_exchange

   !---~---
   !   A node's linearization about its temperature at the step's start (node_lin_t): its
   ! sensible heat through the canopy air's conductance, its vapour live from its own saturation
   ! or held, the vapour's latent heat less the water's own energy (from the datum, so a trace of
   ! water stays exact), and its emission's slope. The water layer holds the water of its vapour
   ! over the step: the balance takes the vapour's latent heat and the water leaves with its own
   ! energy at the node's end state (water_specific_energy), held vapour too (a live vapour past
   ! what the layer holds is held at what empties it). A porous layer's held vapour is the water
   ! network's, its water apart from what the layer keeps over the step (node%water): the vapour
   ! takes its whole enthalpy at the layer's temperature from the layer's energy, and the layer's
   ! phase law reads what is left. So the vapour reaches the canopy air with its own enthalpy, and
   ! the layer's balance does not depend on which phase gave it, which on the freezing curve
   ! would jump by the latent heat of fusion at the curve's onset, where a stage's canopy-air
   ! balance, its ground's sources jumping, could have no root.
   !---~---
   pure subroutine node_linearization(node, env, lin)
      type(surface_node_t), intent(in)  :: node
      type(surface_env_t),  intent(in)  :: env
      type(node_lin_t),     intent(out) :: lin
      lin%dh_dt = sensible_heat_slope(env%ggnet, env%rho_air, env%cp_cas)
      lin%h_n   = lin%dh_dt * (node%temp - env%can_temp)
      if (node%vapour_held) then
         lin%vh = node%vapour_rate
      else
         lin%w_n   = env%ggnet * env%rho_air * (sat_specific_humidity(node%temp, env%press, node%fliq) - env%can_shv)
         lin%dw_dt = env%ggnet * env%rho_air * sat_specific_humidity_temp_deriv(node%temp, env%press, node%fliq)
      end if
      lin%demit_dt = emission_slope(node%emiss, node%temp)
      lin%l_vap    = enthalpy_vapor(node%temp) - water_specific_energy(node%temp, node%fliq)
      lin%l_held   = lin%l_vap
      lin%dl_held  = 0.0_wp
      if (node%on_curve) then
         lin%l_held  = enthalpy_vapor(node%temp)
         lin%dl_held = cp_vap
      end if
   end subroutine node_linearization

   !---~---
   !   A node's fluxes at its temperature t and liquid fraction fl, the node above at t_above and
   ! the one below at t_below, per unit ground area, positive away from the node: its net
   ! radiation, sensible heat, vapour and its latent heat with the leaving water's own energy, its
   ! conduction to the node below, and its net gain, the water passing through it included.
   !---~---
   pure subroutine node_fluxes(node, lin, t, fl, flux)
      type(surface_node_t), intent(in)  :: node
      type(node_lin_t),     intent(in)  :: lin
      real(wp),             intent(in)  :: t, fl
      type(node_flux_t),    intent(out) :: flux
      real(wp) :: dtemp
      dtemp = t - node%temp
      flux%t_surf    = t
      flux%rad       = node%area * (node%rad - lin%demit_dt * dtemp)
      flux%h         = node%area * (lin%h_n + lin%dh_dt * dtemp)
      flux%w_flux    = node%area * (lin%w_n + lin%dw_dt * dtemp) + lin%vh
      flux%le        = flux%w_flux * (lin%l_vap + water_specific_energy(t, fl))
      if (node%on_curve) flux%le = (flux%w_flux - lin%vh) * (lin%l_vap + water_specific_energy(t, fl))           &
                                   + lin%vh * (lin%l_held + lin%dl_held * dtemp)
      flux%g_below   = node%g_below * (t - node%t_below)
      flux%energy_in = flux%rad - flux%h - flux%le + link_heat(node, t, node%t_above, node%t_below)               &
                       - node%m_out * internal_energy_liquid(t)
   end subroutine node_fluxes

   !---~---
   !   A node's temperature, liquid fraction, dT/dE and heat capacity [J/m2/K] at an energy e
   ! [J/m2] with dry heat capacity cd [J/m2/K]: the water layer's on pure water's law
   ! (water_phase), a porous layer's on its freezing curve (porous_energy_state, per volume of its
   ! depth), its temperature fixed at t_3ple while water it stores beyond saturation melts.
   !---~---
   pure subroutine node_phase(node, e, cd, t, fl, dtde, cap)
      type(surface_node_t), intent(in)  :: node
      real(wp),             intent(in)  :: e, cd
      real(wp),             intent(out) :: t, fl, dtde, cap
      real(wp) :: de
      if (.not. node%on_curve) then
         call water_phase(e, node%water, cd, t, fl, dtde, cap)
         return
      end if
      call porous_energy_state(e / node%depth, node%water / (rho_h2o * node%depth), cd / node%depth, node%curve,  &
                               t, fl, de)
      cap  = de * node%depth
      dtde = 1.0_wp / max(cap, tiny_num)
      if (t == t_3ple .and. fl < 1.0_wp) dtde = 0.0_wp
   end subroutine node_phase

   !---~---
   !   A porous node's balance over the step, its neighbours held, on its freezing curve: the root e
   ! of g(e) = e - e0 - dt (r0 + drdt (T(e) - T^n)), and its temperature and liquid fraction there.
   ! g rises with e at slope 1 or more (drdt < 0, and T rises with e), so a bracket lies one
   ! residual from e0; Newton's step is taken where it falls inside the bracket, the Illinois
   ! method's point where it does not.
   !---~---
   pure subroutine curve_node_root(node, cd, e0, r0, drdt, dt, t, fl)
      type(surface_node_t), intent(in)  :: node
      real(wp),             intent(in)  :: cd, e0, r0, drdt, dt   ! [J/m2/K], [J/m2], [W/m2], [W/m2/K], [s]
      real(wp),             intent(out) :: t, fl                  ! [K], [-]
      type(root_bracket_t) :: br
      real(wp)    :: e, g, dtde, cap, e_new, t_end, fl_end, d_end, c_end
      integer(ik) :: it
      e = e0
      call node_balance_on_curve(node, cd, e0, r0, drdt, dt, e, g, t, fl, dtde, cap)
      if (abs(g) <= STACK_TEMP_TOL * (cap + dt * abs(drdt))) return
      br%a = min(e0, e0 - g)
      br%b = max(e0, e0 - g)
      call node_balance_on_curve(node, cd, e0, r0, drdt, dt, br%a, br%fa, t_end, fl_end, d_end, c_end)
      call node_balance_on_curve(node, cd, e0, r0, drdt, dt, br%b, br%fb, t_end, fl_end, d_end, c_end)
      do it = 1_ik, STACK_ROOT_MAX
         e_new = e - g / (1.0_wp - dt * drdt * dtde)
         if (.not. (e_new > br%a .and. e_new < br%b)) e_new = regula_falsi_point(br)
         e = e_new
         call node_balance_on_curve(node, cd, e0, r0, drdt, dt, e, g, t, fl, dtde, cap)
         if (abs(g) <= STACK_TEMP_TOL * (cap + dt * abs(drdt))) exit
         call illinois_narrow(br, e, g)
      end do
   end subroutine curve_node_root

   !--- curve_node_root's residual g at energy e, with the node's temperature, liquid fraction, dT/dE and heat capacity there.
   pure subroutine node_balance_on_curve(node, cd, e0, r0, drdt, dt, e, g, t, fl, dtde, cap)
      type(surface_node_t), intent(in)  :: node
      real(wp),             intent(in)  :: cd, e0, r0, drdt, dt, e
      real(wp),             intent(out) :: g, t, fl, dtde, cap
      call node_phase(node, e, cd, t, fl, dtde, cap)
      g = e - e0 - dt * (r0 + drdt * (t - node%temp))
   end subroutine node_balance_on_curve

   !---~---
   !   Water with dry heat (internal_energy_to_temp) at an energy: its temperature, liquid
   ! fraction, dT/dE (0 on the plateau at the triple point) and the heat capacity of its line, the ice line's on
   ! the plateau [J/m2/K].
   !---~---
   elemental pure subroutine water_phase(e, w, cd, t, fl, dtde, cap)
      real(wp), intent(in)  :: e, w, cd   ! [J/m2], [kg/m2], [J/m2/K]
      real(wp), intent(out) :: t, fl, dtde, cap
      call internal_energy_to_temp(e, w, cd, t, fl)
      cap  = cd + w * cp_ice
      dtde = 0.0_wp
      if (e < (cd + w * cp_ice) * t_3ple) then
         dtde = 1.0_wp / max(cap, tiny_num)
      else if (fl >= 1.0_wp) then
         cap  = cd + w * cp_liq
         dtde = 1.0_wp / max(cap, tiny_num)
      end if
   end subroutine water_phase

   !--- [J/kg] Water's specific internal energy at a temperature and liquid fraction (the inverter's datum).
   elemental function water_specific_energy(temp, fliq) result(u)
      real(wp), intent(in) :: temp, fliq   ! [K], [-]
      real(wp) :: u
      u = fliq * internal_energy_liquid(temp) + (1.0_wp - fliq) * cp_ice * temp
   end function water_specific_energy

   !--- The water layer: its stage before a step, its cover, and its commit after.

   !---~---
   !   water_layer_stage -- the water layer before a step (docs/science/ground_surface_layers.md §2):
   ! the throughfall's rain and snow and the water the tissues shed join it with their enthalpy,
   ! the snow adding its depth at the fresh density, and the matrix follows the ice the store then
   ! holds (matrix_depth_update: rain that freezes in cold snow fills its
   ! pores). From min_depth deep the store is a layer, its free water the liquid beyond what its
   ! ice holds; thinner it is a film, whose liquid is all free and which comes to the temperature
   ! of the node it lies on by a matched heat move (film_heat, from that node): the top porous
   ! layer (below, with its dry heat capacity, depth and freezing curve) and the film reach one
   ! temperature together (film_on_porous_layer),
   ! since a thin litter holds too little heat to melt a dusting alone; the top soil layer, whose
   ! heat capacity a film is far below, gives it its temperature and phase (t_below, fliq_below).
   ! The free water is set apart with its energy at the layer's temperature (free, free_temp): it
   ! is what the water network moves over the step. The rest is held. A store all liquid sets all
   ! it holds apart, as does one the water network left a round-off short of empty. A layer
   ! (is_layer), the two together, a pond's water as a snowpack's, exchange with the canopy air
   ! over the layer's cover (ground_surface_exchange); a film does not.
   !---~---
   pure subroutine water_layer_stage(rainfall, snowfall, t_air, shed, shed_enth, dt, params, min_depth, t_below,   &
                                     fliq_below, layer, st, below, below_dry_cap, below_depth, below_curve)
      real(wp),                  intent(in)    :: rainfall, snowfall   ! [kg/m2/s] the throughfall's rain and snow
      real(wp),                  intent(in)    :: t_air                ! [K] the air temperature precipitation is valued at
      real(wp),                  intent(in)    :: shed, shed_enth      ! [kg/m2/s], [W/m2] water the tissues shed
      real(wp),                  intent(in)    :: dt                   ! [s]
      type(water_material_t),    intent(in)    :: params
      real(wp),                  intent(in)    :: min_depth            ! [m] the depth from which it is a layer
      real(wp),                  intent(in)    :: t_below, fliq_below  ! [K], [-] the node it lies on
      type(water_layer_state_t), intent(inout) :: layer                ! the store in; what it holds over the step out
      type(water_stage_out_t),   intent(out)   :: st
      type(ground_store_t), optional, intent(in) :: below              ! a porous layer it lies on: its water and energy
      real(wp),             optional, intent(in) :: below_dry_cap      ! [J/m2/K] that layer's dry heat capacity
      real(wp),             optional, intent(in) :: below_depth        ! [m] its depth
      type(soil_layer_curve_t), optional, intent(in) :: below_curve    ! its freezing curve
      real(wp) :: ice, fl, t_film, e_film, holding, free
      st = water_stage_out_t()
      call diagnose_water_layer(layer)
      st%ice_start   = water_layer_ice(layer)
      st%snowfall    = max(0.0_wp, snowfall) * dt
      st%inflow_enth = (precip_enthalpy(max(0.0_wp, rainfall), max(0.0_wp, snowfall), t_air) + shed_enth) * dt
      layer%water  = layer%water + (max(0.0_wp, rainfall) + max(0.0_wp, snowfall) + shed) * dt
      layer%energy = layer%energy + st%inflow_enth
      call diagnose_water_layer(layer)
      ice = water_layer_ice(layer)
      layer%matrix_depth = matrix_depth_update(layer%matrix_depth, st%ice_start, ice, st%snowfall, params%fresh_density)
      st%is_layer = layer_depth(layer) >= min_depth
      holding = params%liquid_holding_frac
      if (.not. st%is_layer) then
         if (present(below)) then
            !--- A film on a porous layer: the two at one temperature, the layer on its freezing curve.
            call film_on_porous_layer(below%energy + layer%energy, layer%water, below%water, below_dry_cap,         &
                                      below_depth, below_curve, t_film, fl, e_film)
         else
            !--- A film on the soil: its temperature, and its phase there; at the triple point, its share of liquid.
            t_film = t_below
            fl = merge(1.0_wp, 0.0_wp, t_below > t_3ple)
            if (t_below == t_3ple) fl = fliq_below
         end if
         if (.not. present(below)) e_film = layer%water * water_specific_energy(t_film, fl)
         st%film_heat  = e_film - layer%energy
         layer%energy  = e_film
         call diagnose_water_layer(layer)
         layer%matrix_depth = matrix_depth_update(layer%matrix_depth, ice, water_layer_ice(layer), 0.0_wp,          &
                                                  params%fresh_density)
         holding = 0.0_wp
      end if
      st%free_temp = layer%temp
      if (layer%fliq >= 1.0_wp .or. layer%water <= 0.0_wp) then
         call store_move(layer, st%free, layer%water, layer%energy)
      else
         free = water_layer_free_water(layer%water, layer%fliq, holding)
         if (free > 0.0_wp) call store_move(layer, st%free, free, free * internal_energy_liquid(layer%temp))
      end if
      st%ice_held = water_layer_ice(layer)
      call diagnose_water_layer(layer)
   end subroutine water_layer_stage

   !---~---
   !   The water layer as it lies exposed at the stage: what it holds over the step and its free
   ! water together, at their energy, its temperature and liquid fraction read off it, over its
   ! ice matrix: the store whose cover, optics and conductance the ground's radiation and its
   ! stack read.
   !---~---
   pure function water_layer_exposed(held, free) result(layer)
      type(water_layer_state_t), intent(in) :: held
      type(ground_store_t),      intent(in) :: free
      type(water_layer_state_t) :: layer
      layer        = held
      layer%water  = held%water + free%water
      layer%energy = held%energy + free%energy
      call diagnose_water_layer(layer)
   end function water_layer_exposed

   !---~---
   !   [-] The water layer's cover over a step: Niu & Yang (2007)'s, over the depth and the bulk
   ! density of what it holds over the step (water_layer_cover); none for a film, which does not
   ! exchange.
   !---~---
   pure real(wp) function water_layer_cover_of(layer, is_layer, params, ground_rough) result(cover)
      type(water_layer_state_t), intent(in) :: layer
      logical,                   intent(in) :: is_layer
      type(water_material_t),    intent(in) :: params
      real(wp),                  intent(in) :: ground_rough   ! [m] the ground's roughness ([aerodynamics].soil_rough)
      real(wp) :: depth
      cover = 0.0_wp
      depth = layer_depth(layer)
      if (.not. (is_layer .and. depth > 0.0_wp)) return
      cover = water_layer_cover(depth, water_layer_bulk_density(layer%water, depth), params%cover_scale_factor,     &
                                params%cover_density_exponent, ground_rough)
   end function water_layer_cover_of

   !---~---
   !   water_layer_commit -- the water layer after a step: what it held over the step and its free
   ! water back take its exchange with the canopy air (the energy in, the vapour out, and the heat
   ! the free water's flows carried, all in energy_in) and lose the water the network's flows took
   ! from the free water, with the heat energy_in already holds; then the condensate joins with
   ! the enthalpy the canopy air gave up for it, its matrix follows its ice (matrix_depth_update
   ! from ice_held, the ice it held over the step: ice that melted or sublimated takes its share of
   ! the matrix, ice that formed fills its pores), and its temperature and liquid fraction are read
   ! off its energy.
   !---~---
   pure subroutine water_layer_commit(layer, energy_in, vapour_out, free, free_end, condensate, ice_held, params)
      type(water_layer_state_t), intent(inout) :: layer
      real(wp),                  intent(in)    :: energy_in, vapour_out   ! [J/m2] in, [kg/m2] out
      type(ground_store_t),      intent(in)    :: free                    ! the free water set apart at the stage
      real(wp),                  intent(in)    :: free_end                ! [kg/m2] the free water at the step's end
      type(ground_store_t),      intent(in)    :: condensate              ! dew and fog, and their enthalpy
      real(wp),                  intent(in)    :: ice_held                ! [kg/m2]
      type(water_material_t),    intent(in)    :: params
      type(ground_store_t) :: giver
      giver = free
      call store_move(giver, layer, giver%water, giver%energy)
      !---~---
      !   The vapour and the water the network took from the free water leave the layer together. The
      ! network can leave an empty pond a round-off past empty, which is no water.
      !---~---
      call store_take_exchange(layer, energy_in, vapour_out + (free%water - max(free_end, 0.0_wp)))   ! clamp-ok: round-off
      giver = condensate
      call store_move(giver, layer, giver%water, giver%energy)
      call diagnose_water_layer(layer)
      layer%matrix_depth = matrix_depth_update(layer%matrix_depth, ice_held, water_layer_ice(layer), 0.0_wp,        &
                                               params%fresh_density)
   end subroutine water_layer_commit

   !--- [kg/m2] The water layer's ice, from its water and its liquid fraction.
   elemental real(wp) function water_layer_ice(layer) result(ice)
      type(water_layer_state_t), intent(in) :: layer
      ice = (1.0_wp - layer%fliq) * layer%water
   end function water_layer_ice

   !--- [m] The water layer's depth (water_layer_depth): its matrix's, or its ice and liquid where more.
   elemental real(wp) function layer_depth(layer) result(depth)
      type(water_layer_state_t), intent(in) :: layer
      depth = water_layer_depth(layer%matrix_depth, water_layer_ice(layer), layer%fliq * layer%water)
   end function layer_depth

   !--- The water layer's temperature and liquid fraction from its energy; the triple point, frozen, when it holds no water.
   elemental subroutine diagnose_water_layer(layer)
      type(water_layer_state_t), intent(inout) :: layer
      if (layer%water > 0.0_wp) then
         call internal_energy_to_temp(layer%energy, layer%water, 0.0_wp, layer%temp, layer%fliq)
      else
         layer%temp = t_3ple
         layer%fliq = 0.0_wp
      end if
   end subroutine diagnose_water_layer

end module ed_surface
