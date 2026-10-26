! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_params_config -- the config blocks reach the records the model runs on.             !
!   1. THE BLOCK LOADERS (ed_config_io): [soil], [fast], [snow] and [aerodynamics] from a    !
!      small TOML file. String selectors map to their codes, scalar keys override, keys absent !
!      from a block keep their type default, a file with no [soil] block keeps the geothermal  !
!      bottom, the drainage boundary's two named ends, and the parameter record lists every key !
!      the loaders read, set or defaulted, and none they did not.                             !
!   2. THE [soil] COLUMN (ed_soil_derive::build_soil_profile): the default column is          !
!      SOIL_LAYER_TEMPLATE (SoilGrids' depths, then deeper) cut at the 2 m default depth; depth !
!      cuts the same list; a short remainder joins the layer above; layer_bottoms replaces the  !
!      list; and the texture, retention and thermal keys reach every layer of the profile       !
!      (docs/science/soil_biophysics.md, The soil layers). A mis-threaded key leaves the model  !
!      running on a default with every budget closed, so these assert that the keys turn.      !
!==========================================================================================!
program test_params_config
   use ed_test_assert,         only : check, check_true, test_report
   use ed_kinds,               only : wp, ik
   use ed_constants,           only : n_soil_layer_max
   use ed_biophysics_params,   only : snow_params_t, aero_params_t
   use ed_soil_params,         only : soil_params_t, soil_profile_t, SOIL_BC_AQUIFER, SOIL_BC_DRAIN,      &
                                      ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   use ed_soil_derive,         only : soil_layers, build_soil_profile
   use ed_soil_hydr_lib,       only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use ed_config,              only : ed_config_t
   use ed_toml,                only : toml_table_t, toml_parse_file, toml_write_record
   use ed_config_io,           only : load_fast_params, load_soil, load_snow_params, load_aero_params
   implicit none

   call test_block_loaders()
   call test_soil_column()

   call test_report('test_params_config')

contains

   !----- 1. The block loaders. -------------------------------------------------------------!
   subroutine test_block_loaders()
      character(len=*), parameter :: TOMLFILE = 'test_params_config_tmp.toml'
      type(toml_table_t)   :: tm
      type(soil_params_t)  :: s
      type(ed_config_t)    :: c
      type(snow_params_t)  :: sn
      type(aero_params_t)  :: a
      integer(ik) :: u
      logical     :: ok

      !----- Write a fixture that overrides a SUBSET of each block (leaving others to default). ----!
      open(newunit=u, file=TOMLFILE, status='replace', action='write')
      write(u,'(a)') '[soil]'
      write(u,'(a)') 'water_bottom_bc   = "aquifer"'
      write(u,'(a)') 'w_pond_max        = 9.0'
      write(u,'(a)') 'ground_albedo_nir = 0.25'
      write(u,'(a)') 'thermal_bottom_bc = "dirichlet"'
      write(u,'(a)') 'deep_temp         = 283.24'
      write(u,'(a)') ''
      write(u,'(a)') '[fast]'
      write(u,'(a)') 'water_atol              = 2.5e-5'
      write(u,'(a)') 'soil_energy_budget_atol = 0.05'
      write(u,'(a)') ''
      write(u,'(a)') '[snow]'
      write(u,'(a)') 'rho_snow   = 300.0'
      write(u,'(a)') 'snow_emiss = 0.95'
      write(u,'(a)') ''
      write(u,'(a)') '[aerodynamics]'
      write(u,'(a)') 'z0m_ratio = 0.10'
      write(u,'(a)') 'n_iter_mo = 6'
      write(u,'(a)') 'z0m_ration = 0.20'             ! a misspelling: the parameter record must not list it
      close(u)

      call toml_parse_file(TOMLFILE, tm, ok)
      call check_true('TOML fixture parsed', ok, 0.0_wp)

      call load_fast_params(tm, c)
      call load_soil(tm, s)
      call load_snow_params(tm, sn)
      call load_aero_params(tm, a)

      print '(a)', 'test_block_loaders:'

      !----- [soil] and [fast]: string selectors -> enums, scalar overrides, and absent-keeps-default. -!
      call check_true('soil.water_bottom_bc "aquifer" -> SOIL_BC_AQUIFER', s%water_bottom_bc == SOIL_BC_AQUIFER,   &
                      real(s%water_bottom_bc, wp))
      call check('fast.water_atol overridden', c%water_atol,     2.5e-5_wp, 1.0e-12_wp)
      call check('soil.w_pond_max overridden', s%w_pond_max, 9.0_wp, 0.0_wp)
      call check('soil.ground_albedo_nir overridden', s%ground_albedo_nir, 0.25_wp, 1.0e-12_wp)
      call check('soil.ground_albedo_vis absent -> default 0.15', s%ground_albedo_vis, 0.15_wp, 0.0_wp)
      call check('soil.ground_emissivity absent -> default 0.95', s%ground_emissivity, 0.95_wp, 0.0_wp)
      call check('soil.dsl_dmax absent -> default 0.015', s%dsl_dmax, 0.015_wp, 0.0_wp)

      !----- The soil heat's budget-closure threshold, for the debug halt. The soil-thermal solve is   !
      !      one backward-Euler step, so there are no substep keys. ---------------------------------------!
      call check('fast.soil_energy_budget_atol overridden', c%soil_energy_budget_atol, 0.05_wp, 1.0e-12_wp)
      !----- The [soil] thermal bottom BC: selector -> enum, anchor temperature read, and the         !
      !      anchor DEPTH left at its derived default when the key is absent. -------------------------!
      call check_true('soil.thermal_bottom_bc "dirichlet" -> ENERGY_BC_DIRICHLET',                    &
                      s%thermal_bottom_bc == ENERGY_BC_DIRICHLET, real(s%thermal_bottom_bc, wp))
      call check('soil.deep_temp overridden', s%deep_temp, 283.24_wp, 1.0e-9_wp)
      call check('soil.deep_depth absent -> derived default 3.145', s%deep_depth, 3.145_wp, 1.0e-12_wp)

      !----- [snow]: overrides + default. ---------------------------------------------------------!
      call check('snow.rho_snow overridden',   sn%rho_snow,   300.0_wp, 1.0e-9_wp)
      call check('snow.snow_emiss overridden', sn%snow_emiss, 0.95_wp,  1.0e-12_wp)
      call check('snow.k_snow absent -> default', sn%k_snow,  0.15_wp,  1.0e-12_wp)

      !----- [aerodynamics]: override (real + int) + default. -------------------------------------!
      call check('aero.z0m_ratio overridden',   a%z0m_ratio, 0.10_wp, 1.0e-12_wp)
      call check('aero.n_iter_mo overridden',   real(a%n_iter_mo, wp), 6.0_wp, 0.0_wp)
      call check('aero.d_ratio absent -> default', a%d_ratio, 0.63_wp, 1.0e-12_wp)

      !----- A table with NO [soil] block keeps the compiled default bottom boundary, geothermal. -----!
      block
         type(toml_table_t)  :: tm2
         type(soil_params_t) :: s2
         logical             :: ok2
         open(newunit=u, file='test_biophysics_params_empty.toml', status='replace', action='write')
         write(u,'(a)') '[run]'
         write(u,'(a)') 'dt_slow = "1d"'
         close(u)
         call toml_parse_file('test_biophysics_params_empty.toml', tm2, ok2)
         call check_true('empty-block TOML parsed', ok2, 0.0_wp)
         call load_soil(tm2, s2)
         call check_true('no [soil] block -> default ENERGY_BC_GEOTHERMAL',                           &
                         s2%thermal_bottom_bc == ENERGY_BC_GEOTHERMAL, real(s2%thermal_bottom_bc, wp))
         open(newunit=u, file='test_biophysics_params_empty.toml', status='old', action='write')
         close(u, status='delete')
      end block

      !----- The drainage boundary's two named ends (docs/science/soil_biophysics.md §Boundary         !
      !      conditions): free_drain is drainage at f = 1, the default, and bedrock at f = 0. ------------!
      block
         type(toml_table_t)  :: tb
         type(soil_params_t) :: s_free, s_rock, s_default
         logical             :: okb
         character(len=*), parameter :: NAMES(2) = [character(len=10) :: 'free_drain', 'bedrock']
         integer(ik) :: i
         s_default = soil_params_t()                   ! the type's default
         do i = 1_ik, 2_ik
            open(newunit=u, file='test_biophysics_params_bc.toml', status='replace', action='write')
            write(u,'(a)') '[soil]'
            write(u,'(3a)') 'water_bottom_bc = "', trim(NAMES(i)), '"'
            close(u)
            call toml_parse_file('test_biophysics_params_bc.toml', tb, okb)
            if (i == 1_ik) call load_soil(tb, s_free)
            if (i == 2_ik) call load_soil(tb, s_rock)
         end do
         open(newunit=u, file='test_biophysics_params_bc.toml', status='old', action='write')
         close(u, status='delete')
         call check_true('soil.water_bottom_bc "free_drain" -> drainage at f = 1, the default',          &
                         s_free%water_bottom_bc == SOIL_BC_DRAIN .and. s_free%drain_factor == 1.0_wp .and.     &
                         s_default%water_bottom_bc == SOIL_BC_DRAIN .and. s_default%drain_factor == 1.0_wp)
         call check_true('soil.water_bottom_bc "bedrock" -> drainage at f = 0',                          &
                         s_rock%water_bottom_bc == SOIL_BC_DRAIN .and. s_rock%drain_factor == 0.0_wp)
      end block

      !----- The parameter record: every key the loaders read, set or defaulted, and no key they did   !
      !      not read. A calibration trial checks the keys it set against it. ------------------------!
      block
         character(len=*), parameter :: RECFILE = 'test_parameter_record_tmp.csv'
         character(len=1024) :: line, z0m_line
         integer :: ios
         logical :: set_seen, default_seen, typo_seen, header_ok
         call toml_write_record(RECFILE)
         set_seen = .false.
         default_seen = .false.
         typo_seen = .false.
         z0m_line = '(none)'
         open(newunit=u, file=RECFILE, status='old', action='read')
         read(u, '(a)') line
         header_ok = trim(line) == 'source,key,index,present,value'
         do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, ',aerodynamics.z0m_ratio,') > 0) z0m_line = line
            if (index(line, ',aerodynamics.z0m_ratio,0,true,') > 0 .and. index(line, '1.0000000000000001E-001') > 0) &
               set_seen = .true.
            if (index(line, ',aerodynamics.d_ratio,0,false,') > 0) default_seen = .true.
            if (index(line, 'z0m_ration') > 0) typo_seen = .true.
         end do
         close(u, status='delete')
         call check_true('record: header', header_ok, 0.0_wp)
         if (.not. set_seen) write(*,'(2a)') '   z0m_ratio row: ', trim(z0m_line)
         call check_true('record: a key set in the file, with the value used to the last digit', set_seen, 0.0_wp)
         call check_true('record: a key left to its default is listed as a default', default_seen, 0.0_wp)
         call check_true('record: a misspelt key that nothing reads is absent', .not. typo_seen, 0.0_wp)
      end block

      !----- Clean up the fixture. ----------------------------------------------------------------!
      open(newunit=u, file=TOMLFILE, status='old', action='write')
      close(u, status='delete')

   end subroutine test_block_loaders

   !----- 2. The [soil] column, through the production build_soil_profile. ------------------!
   subroutine test_soil_column()
      type(soil_params_t)  :: sc
      type(soil_profile_t) :: p, p2
      real(wp)             :: bottoms(n_soil_layer_max)
      integer(ik)          :: n_layer, k
      real(wp), parameter  :: SOILGRIDS_BOTTOMS(5) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp]
      print '(a)', 'test_soil_column:'

      !----- The defaults: a 2 m van Genuchten column on the template layers. ---------------!
      call check_true('the default layers are the template', sc%n_layer_bottoms == 0_ik)
      call check('the default depth is 2 m', sc%depth, 2.0_wp, 0.0_wp)
      call check('the default theta_sat', sc%theta_sat, 0.43_wp, 0.0_wp)
      call check('the default ksat', sc%ksat, 2.89e-6_wp, 0.0_wp)
      call check_true('the default retention is van Genuchten', sc%retention == SOIL_RETENTION_VG)

      !----- The default column is the template cut at 2 m: seven layers, a 5 cm top layer,  !
      !      SoilGrids' bottoms to 1 m, and the base at -depth exactly. ---------------------!
      call build_soil_profile(sc, p)
      call check_true('the default 2 m column is 7 layers', p%n_active == 7_ik, real(p%n_active, wp))
      call check('the top layer is 5 cm', p%dz(1), 0.05_wp, 1.0e-12_wp)
      do k = 1_ik, 5_ik
         call check('a layer bottom is on SoilGrids'' depths', -p%soil_layer_z(k + 1_ik), SOILGRIDS_BOTTOMS(k), 1.0e-12_wp)
      end do
      call check('the column base is at -depth', p%soil_layer_z(p%n_active + 1_ik), -sc%depth, 1.0e-12_wp)

      !----- depth cuts the same list: 0.3 m is 5, 10 and 15 cm; 4 m is nine layers whose top is  !
      !      still 5 cm, the last running from 3 to 4 m. -------------------------------------!
      sc%depth = 0.3_wp
      call build_soil_profile(sc, p2)
      call check_true('a 0.3 m column is 5, 10 and 15 cm', p2%n_active == 3_ik .and. abs(p2%dz(3) - 0.15_wp) < 1.0e-12_wp)
      sc%depth = 4.0_wp
      call build_soil_profile(sc, p2)
      call check_true('a 4 m column is 9 layers', p2%n_active == 9_ik, real(p2%n_active, wp))
      call check('a deeper column keeps the 5 cm top layer', p2%dz(1), 0.05_wp, 1.0e-12_wp)
      call check('the 4 m column ends in a 3-4 m layer', p2%dz(9), 1.0_wp, 1.0e-12_wp)

      !----- A remainder under half of the list layer it falls in joins the layer above: 0.65 m  !
      !      ends in a 30-65 cm layer, while 0.85 m keeps a 60-85 cm layer of its own. ---------!
      sc%depth = 0.65_wp
      call soil_layers(sc, bottoms, n_layer)
      call check_true('a 5 cm remainder joins the layer above', n_layer == 4_ik .and.                      &
                      abs(bottoms(4) - 0.65_wp) < 1.0e-12_wp .and. abs(bottoms(3) - 0.30_wp) < 1.0e-12_wp)
      sc%depth = 0.85_wp
      call soil_layers(sc, bottoms, n_layer)
      call check_true('a 25 cm remainder, over half its layer, stays a layer',                              &
                      n_layer == 5_ik .and. abs(bottoms(5) - 0.85_wp) < 1.0e-12_wp)

      !----- layer_bottoms replaces the list, and depth still cuts it. -----------------------!
      sc%n_layer_bottoms = 3_ik
      sc%layer_bottoms(1:3) = [0.2_wp, 0.5_wp, 1.0_wp]
      sc%depth = 1.0_wp
      call build_soil_profile(sc, p2)
      call check_true('layer_bottoms reaches the column', p2%n_active == 3_ik .and. abs(p2%dz(1) - 0.2_wp) < 1.0e-12_wp)
      sc%depth = 0.8_wp
      call build_soil_profile(sc, p2)
      call check_true('depth cuts the layer_bottoms list',                                                &
                      p2%n_active == 3_ik .and. abs(p2%soil_layer_z(4) + 0.8_wp) < 1.0e-12_wp)
      sc%n_layer_bottoms = 0_ik
      sc%layer_bottoms = 0.0_wp
      sc%depth = 2.0_wp

      !----- Texture reaches the layers, and the derived field capacity follows it. -----------!
      sc%theta_sat = 0.55_wp
      sc%theta_res = 0.02_wp
      call build_soil_profile(sc, p2)
      call check('theta_sat reaches the layers', p2%theta_sat(1), 0.55_wp, 1.0e-14_wp)
      call check('theta_res reaches the layers', p2%theta_res(1), 0.02_wp, 1.0e-14_wp)
      call check_true('field capacity lies between residual content and saturation',                      &
                      p2%theta_fc(1) > p2%theta_res(1) .and. p2%theta_fc(1) < p2%theta_sat(1), p2%theta_fc(1))
      sc%theta_sat = 0.43_wp
      sc%theta_res = 0.078_wp

      !----- The retention family selects what curve_par_a/n mean: for Campbell, the air-entry    !
      !      suction psi_sat [m] (< 0) and the exponent b, as validate_config requires. ---------!
      sc%retention = SOIL_RETENTION_CAMPBELL
      sc%curve_par_a = -0.26_wp
      sc%curve_par_n = 5.65_wp
      call build_soil_profile(sc, p2)
      call check_true('the retention family reaches the column', p2%retention == SOIL_RETENTION_CAMPBELL)
      call check_true('Campbell field capacity lies between residual content and saturation',             &
                      p2%theta_fc(1) > p2%theta_res(1) .and. p2%theta_fc(1) < p2%theta_sat(1), p2%theta_fc(1))
      sc%retention = SOIL_RETENTION_VG
      sc%curve_par_a = 3.6_wp                                   ! the van Genuchten defaults
      sc%curve_par_n = 1.56_wp

      !----- The thermal keys reach every layer of the profile. -------------------------------!
      sc%solid_conductivity = 2.5_wp
      sc%dry_conductivity = 0.20_wp
      sc%dry_heat_capacity = 1.5e6_wp
      call build_soil_profile(sc, p2)
      call check('solid_conductivity reaches the top layer', p2%solid_conductivity(1), 2.5_wp, 1.0e-14_wp)
      call check('dry_conductivity reaches the top layer', p2%dry_conductivity(1), 0.20_wp, 1.0e-14_wp)
      call check('dry_heat_capacity reaches the bottom layer', p2%dry_heat_capacity(p2%n_active), 1.5e6_wp, 1.0e-8_wp)
   end subroutine test_soil_column

end program test_params_config
