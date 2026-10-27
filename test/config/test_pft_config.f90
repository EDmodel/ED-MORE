! SPDX-License-Identifier: Apache-2.0
!---~---
!   test_pft_config -- the PFT file's loader (load_ed_config's PFT half) and the per-PFT
! tables derive_parameters builds from it.
!   1. THE SHIPPED FILE: ed_config_pft.toml, loaded through ed_config_main.toml, gives
!      as many PFTs as wood_density lists; each [pft] optics key and each [hydraulics] key
!      reaches every PFT's field (one value broadcast, or one per PFT), and the two-stream
!      table holds each shipped PFT's optics.
!   2-5. On the test fixture: the optics, the clumping and the leaf angles, and the
!      hydraulics arrive per PFT and DIFFERENTIATE the PFTs; the ground optics come from
!      [soil]; a cohort is born with its water on its PFT's curve, and the slow tier caps
!      it there.
!   6. The PFT-parameter CSV has as many values as column headers.
!
!   The fourteen optics traits are what lets two plant functional types differ in how they
! intercept or scatter light -- the competition the two-stream solver exists to resolve.
! They are only worth anything if the values arrive per PFT, and that wiring fails
! silently: a mis-threaded index leaves every PFT on the first one's optics, the radiation
! still solves, the energy budget still closes, and the canopy is quietly one optical
! type. So this asserts the values arrive AND that PFTs differ, not merely that the code
! compiles. The same holds for the per-PFT hydraulics, down to a cohort's water at birth
! and the slow tier's shed of water above capacity.
!
!   It runs from the source tree, where ed_config_main.toml names ed_config_pft.toml,
! and writes only into the directory its first argument names.
!---~---
program test_pft_config
   use ed_kinds,            only : wp, ik
   use ed_config,           only : ed_config_t, derive_parameters
   use ed_pft_derive,       only : pft_leaf_psi_tlp
   use ed_config_io,        only : load_ed_config, write_pft_params_csv
   use ed_toml,             only : toml_table_t, toml_parse_file, toml_has, toml_is_array, toml_real, toml_real_array
   use ed_pft_params,       only : hydraulics_table_t, RAD_VIS, RAD_NIR, RAD_LW
   use ed_constants,        only : stefan_boltzmann
   use ed_met_types,        only : met_forcing_t
   use ed_fast_integ_types,     only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t,       &
                                   ensure_static_capacity, ensure_forcing_capacity
   use ed_fast_forcing,     only : patch_forcing
   use ed_plant_hydr_lib,   only : conductance_retained_fraction, tissue_psi_tlp,  &
                                   tissue_psi_to_water
   use ed_polygon_state_vars, only : polygon_t
   use ed_patch_state_vars, only : patch_physical_state
   use ed_soil_params,    only : soil_profile_t
   use ed_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use ed_slow_state, only : shed_excess_tissue_water
   use ed_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(ed_config_t)  :: cfg
   real(wp) :: om1, om2
   character(len=1024) :: work_dir   !--- where the CSV of section 6 is written

   call banner('the PFT file: its loader, and its per-PFT tables')
   call get_command_argument(1, work_dir)
   call check(len_trim(work_dir) > 0, 'the work directory was passed on the command line')

   !--- 1. The shipped PFT file through its loader.
   call test_shipped_pft_file()

   cfg = build_test_config()

   !---~---
   !   2. Uniform traits => a uniform table. The control for test 3. The tables are the
   ! configuration's, derived from its traits (derive_parameters), so each change below
   ! derives them again.
   !---~---
   call check(cfg%pft%optics%n_pft >= 2_ik, 'fixture must have at least two PFTs to differentiate')
   call check_close(cfg%pft%optics%omega_leaf(RAD_VIS,1), cfg%pft%optics%omega_leaf(RAD_VIS,2), 1.0e-14_wp,       &
                    'identical traits must give identical optics')

   !---~---
   !   3. Differentiate ONE PFT and the table must follow it -- only that one.
   ! Single-scatter albedo is reflectance + transmittance, so raising the reflectance
   ! of PFT 2 must raise ITS omega and leave PFT 1 alone.
   !---~---
   om1 = cfg%pft%optics%omega_leaf(RAD_VIS,1)
   cfg%pft%leaf_reflect_vis(2) = 0.30_wp
   call derive_parameters(cfg)
   om2 = cfg%pft%optics%omega_leaf(RAD_VIS,2)
   call check_close(cfg%pft%optics%omega_leaf(RAD_VIS,1), om1, 1.0e-14_wp,                                    &
                    'changing PFT 2 must not move PFT 1 (per-PFT indexing)')
   call check_close(om2, 0.30_wp + cfg%pft%leaf_transmit_vis(2), 1.0e-12_wp,                               &
                    'leaf VIS omega must be reflectance + transmittance of THAT PFT')
   call check(om2 > om1, 'a more reflective PFT must have the larger single-scatter albedo')

   !--- 4. Longwave is configured as EMISSIVITY: omega = 1 - emissivity, transmittance 0.
   cfg%pft%leaf_emissivity(2) = 0.80_wp
   call derive_parameters(cfg)
   call check_close(cfg%pft%optics%omega_leaf(RAD_LW,2), 1.0_wp - 0.80_wp, 1.0e-12_wp,                        &
                    'leaf LW omega must be 1 - emissivity (a leaf is opaque in the thermal band)')

   !---~---
   !   4b. The ground optics come from [soil]: a bare patch reflects each shortwave band by its
   ! albedo and returns the sky's longwave by 1 - emissivity, plus its own emission.
   !---~---
   cfg%soil%ground_albedo_vis = 0.11_wp
   cfg%soil%ground_albedo_nir = 0.22_wp
   cfg%soil%ground_emissivity = 0.93_wp
   block
      type(fast_static_state_t)  :: bare
      type(fast_dynamic_state_t) :: dyn
      type(fast_forcing_t)       :: forcing
      type(met_forcing_t)        :: met
      real(wp)                   :: lw_up
      call ensure_static_capacity(bare, 0_ik)
      call ensure_forcing_capacity(forcing, 0_ik)
      dyn%soil_energy%soil_temp(1) = 290.0_wp
      call patch_forcing(cfg, bare, dyn, met, .false., forcing)
      call check_close(forcing%sw_up_vis, 0.11_wp * forcing%sw_in_vis, 1.0e-12_wp,                            &
                       'soil.ground_albedo_vis did not reach the radiation solver')
      call check_close(forcing%sw_up_nir, 0.22_wp * forcing%sw_in_nir, 1.0e-12_wp,                            &
                       'soil.ground_albedo_nir did not reach the radiation solver')
      lw_up = (1.0_wp - 0.93_wp) * met%lw_down + 0.93_wp * stefan_boltzmann * 290.0_wp**4
      call check_close(forcing%lw_up, lw_up, 1.0e-9_wp, 'soil.ground_emissivity did not reach the radiation solver')
   end block
   cfg%soil%ground_albedo_vis = 0.15_wp
   cfg%soil%ground_albedo_nir = 0.30_wp
   cfg%soil%ground_emissivity = 0.95_wp

   !--- 5. Clumping and the leaf-angle distribution are per PFT too.
   cfg%pft%leaf_clumping(2)   = 0.40_wp
   cfg%pft%leaf_angle_mean(2) = 20.0_wp        ! erectophile against PFT 1's 45 deg (spherical is 57.3)
   call derive_parameters(cfg)
   call check_close(cfg%pft%optics%clumping_leaf(2), 0.40_wp, 1.0e-14_wp, 'leaf clumping did not reach the table')
   call check_close(cfg%pft%optics%clumping_leaf(1), 0.80_wp, 1.0e-14_wp, 'PFT 1 clumping must be untouched')
   call check(abs(cfg%pft%optics%bf(2) - cfg%pft%optics%bf(1)) > 1.0e-6_wp,                            &
              'a different leaf-angle mean must give a different angle distribution')

   !---~---
   !   5b. PER-PFT HYDRAULIC TRAITS. Same shape as the optics checks above: each PFT's table
   ! entry carries its own [hydraulics] record, and a different record CHANGES the answer,
   ! not merely that it compiles.
   !---~---
   block
      type(hydraulics_table_t)   :: tab
      real(wp) :: plc1, plc2
      cfg%pft%hydro(2)%wood_psi50     = -0.8_wp   ! vs PFT 1's -2.0
      cfg%pft%hydro(2)%leaf_curve%pi0 = -2.5_wp   ! vs PFT 1's -1.5
      cfg%pft%hydro(2)%wood_kmax      =  2.0_wp
      call derive_parameters(cfg)
      tab = cfg%pft%hydraulics
      call check_close(tab%pft(1)%wood_psi50, cfg%pft%hydro(1)%wood_psi50, 0.0_wp,                &
                       'PFT 1 keeps its own psi50')
      call check_close(tab%pft(2)%wood_psi50, -0.8_wp,   0.0_wp, 'PFT 2 takes its own psi50')
      call check_close(tab%pft(2)%leaf_curve%pi0, -2.5_wp,   0.0_wp, 'PFT 2 takes its own leaf pi0')
      call check_close(tab%pft(2)%wood_kmax, 2.0_wp, 0.0_wp, 'PFT 2 takes its own conductivity')
      call check_close(tab%pft(1)%wood_kmax, cfg%pft%hydro(1)%wood_kmax, 0.0_wp,                  &
                       'PFT 1 keeps its own conductivity')
      call check_close(tab%pft(2)%vessel_curl, cfg%pft%hydro(2)%vessel_curl, 0.0_wp,              &
                       'every trait of the table is the PFT''s own')
      !--- And it MATTERS: at the same potential the more vulnerable xylem has lost more.
      plc1 = 1.0_wp - conductance_retained_fraction(-1.5_wp, tab%pft(1)%wood_psi50, tab%pft(1)%wood_kexp)
      plc2 = 1.0_wp - conductance_retained_fraction(-1.5_wp, tab%pft(2)%wood_psi50, tab%pft(2)%wood_kexp)
      call check(plc2 > plc1, 'a more vulnerable PFT loses more conductance at the same psi')
      !---~---
      !   The fast water-stress limb (the turgor-loss point to twice it) follows the PFT's OWN
      ! pressure-volume curve, the one the plant-water solver stores its leaf water on.
      !---~---
      call check_close(cfg%pft%leaf_photo%pft(2)%psi_open,                                   &
                       tissue_psi_tlp(-2.5_wp, cfg%pft%hydro(2)%leaf_curve%elastic_mod), 1.0e-12_wp,   &
                       'PFT 2 loses capacity from its own turgor-loss point')
      call check_close(cfg%pft%leaf_photo%pft(2)%psi_close,                                  &
                       2.0_wp * tissue_psi_tlp(-2.5_wp, cfg%pft%hydro(2)%leaf_curve%elastic_mod),        &
                       1.0e-12_wp, 'PFT 2 has lost it at twice its turgor-loss point')
      call check_close(cfg%pft%leaf_photo%pft(1)%psi_open,                                   &
                       tissue_psi_tlp(cfg%pft%hydro(1)%leaf_curve%pi0, cfg%pft%hydro(1)%leaf_curve%elastic_mod), &
                       1.0e-12_wp, 'PFT 1 keeps its own turgor-loss point')
      call check_close(pft_leaf_psi_tlp(cfg%pft, 2_ik), cfg%pft%leaf_photo%pft(2)%psi_open, 0.0_wp, &
                       'the phenology reads the same per-PFT turgor-loss point')
   end block

   !---~---
   !   A cohort is born with its tissue water on ITS PFT's curve, the one the fast loop reads
   ! the water back on, and the slow tier caps a store on that curve too.
   !---~---
   block
      type(hydraulics_table_t)   :: tab
      type(polygon_t) :: poly
      real(wp) :: w_expect, w_shed, t_birth, psi_birth
      integer(ik) :: i2
      type(soil_profile_t) :: soil
      cfg%pft%hydro(2)%leaf_curve%water_sat = 3.0_wp   ! vs PFT 1's 2.0
      call derive_parameters(cfg)
      tab = cfg%pft%hydraulics
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 20.0_wp)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 20.0_wp)
      call finalize_init(poly)
      i2 = merge(1_ik, 2_ik, poly%cohort%pft(1) == 2_ik)        ! the PFT-2 cohort, after the sort
      call patch_physical_state(poly%patch%cas(1), poly%patch%soil_water(1), poly%site%surface(1),          &
                                t_birth, psi_birth)
      w_expect = tissue_psi_to_water(psi_birth, tab%pft(2)%leaf_curve, poly%cohort%leaf_carbon(i2))
      call check_close(poly%cohort%leaf_water_mass(i2), w_expect, 1.0e-12_wp,                      &
                       'PFT 2 leaf water at birth is on its own curve')
      call check(abs(w_expect - tissue_psi_to_water(psi_birth, tab%pft(1)%leaf_curve, poly%cohort%leaf_carbon(i2))) &
                 > 1.0e-6_wp * w_expect, 'the two PFTs'' curves give different birth water (the fixture bites)')
      !--- Capped on its own ceiling too: a supersaturated store falls to 3.0 x leaf carbon.
      poly%cohort%leaf_water_mass(i2) = 10.0_wp * poly%cohort%leaf_carbon(i2)
      call shed_excess_tissue_water(poly, cfg, w_shed)
      call check_close(poly%cohort%leaf_water_mass(i2), 3.0_wp * poly%cohort%leaf_carbon(i2), 1.0e-12_wp, &
                       'PFT 2 leaf water is capped at its own saturation')
   end block

   !---~---
   !   6. The PFT-parameter CSV dump has as many values as it has column headers.
   ! A Fortran format SHORTER than its output list does not fail: it REVERTS to the last
   ! repeat group and keeps going, so an integer edit descriptor silently receives a real
   ! and prints its bit pattern, and the trailing columns disappear. Adding PFT traits
   ! without widening the format does exactly that, and nothing else in the suite notices,
   ! because the run completes and every budget still closes. Counting the two lines is
   ! the whole test.
   !---~---
   block
      character(len=:), allocatable :: csvfile
      character(len=8192) :: head_line, val_line
      integer(ik)         :: u, ios, nhead, nval
      csvfile = trim(work_dir)//'/test_pft_params_tmp.csv'
      call derive_parameters(cfg)                ! keep cfg consistent with the edits above
      call write_pft_params_csv(cfg, csvfile)
      open(newunit=u, file=csvfile, status='old', action='read', iostat=ios)
      call check(ios == 0, 'PFT parameter CSV must be written')
      read(u,'(a)') head_line
      read(u,'(a)') val_line
      close(u, status='delete')
      nhead = count_fields(head_line)
      nval  = count_fields(val_line)
      call check(nhead == nval, 'PFT CSV: header column count must equal the value count')
      if (nhead /= nval) print '(a,i0,a,i0)', '   header fields = ', nhead, ', value fields = ', nval
      !---~---
      !   And no field may be a bit pattern: every value is either an integer or es15.8, so a
      ! 19-digit run of digits is the format-reversion signature, not a real number.
      !---~---
      call check(index(val_line, '4605380978949069210') == 0,                                     &
                 'PFT CSV: a real printed through an integer descriptor (format reversion)')
   end block

   print '(a)', 'test_pft_config: all checks passed'

contains

   !---~---
   !   1. ed_config_main.toml and the PFT file it names, loaded as a run loads them, against
   ! the PFT file read on its own: the PFT count, every optics key of [pft] and every key of
   ! [hydraulics] for every PFT, and the two-stream table's optics for every PFT.
   !---~---
   subroutine test_shipped_pft_file()
      type(ed_config_t)   :: shipped
      type(toml_table_t)  :: t
      real(wp)    :: buf(64)
      integer(ik) :: n, ip
      logical     :: ok
      call load_ed_config('ed_config_main.toml', shipped)
      call toml_parse_file(trim(shipped%init%pft_config), t, ok)
      call check(ok, 'the PFT file ed_config_main.toml names parses on its own')
      call toml_real_array(t, 'pft.wood_density', buf, n)
      call check(n > 1_ik .and. shipped%pft%n == n, 'the loader gives as many PFTs as pft.wood_density lists')
      associate (p => shipped%pft)
         call check_key(t, 'pft.leaf_reflect_vis',  p%leaf_reflect_vis,  p%n)
         call check_key(t, 'pft.leaf_transmit_vis', p%leaf_transmit_vis, p%n)
         call check_key(t, 'pft.leaf_reflect_nir',  p%leaf_reflect_nir,  p%n)
         call check_key(t, 'pft.leaf_transmit_nir', p%leaf_transmit_nir, p%n)
         call check_key(t, 'pft.leaf_emissivity',   p%leaf_emissivity,   p%n)
         call check_key(t, 'pft.wood_reflect_vis',  p%wood_reflect_vis,  p%n)
         call check_key(t, 'pft.wood_transmit_vis', p%wood_transmit_vis, p%n)
         call check_key(t, 'pft.wood_reflect_nir',  p%wood_reflect_nir,  p%n)
         call check_key(t, 'pft.wood_transmit_nir', p%wood_transmit_nir, p%n)
         call check_key(t, 'pft.wood_emissivity',   p%wood_emissivity,   p%n)
         call check_key(t, 'pft.leaf_clumping',     p%leaf_clumping,     p%n)
         call check_key(t, 'pft.wood_clumping',     p%wood_clumping,     p%n)
         call check_key(t, 'pft.leaf_angle_mean',   p%leaf_angle_mean,   p%n)
         call check_key(t, 'pft.leaf_angle_std',    p%leaf_angle_std,    p%n)
         call check_key(t, 'hydraulics.leaf_pi0',             p%hydro(1:p%n)%leaf_curve%pi0,           p%n)
         call check_key(t, 'hydraulics.leaf_elastic_mod',     p%hydro(1:p%n)%leaf_curve%elastic_mod,   p%n)
         call check_key(t, 'hydraulics.leaf_apoplast_frac',   p%hydro(1:p%n)%leaf_curve%apoplast_frac, p%n)
         call check_key(t, 'hydraulics.leaf_water_sat',       p%hydro(1:p%n)%leaf_curve%water_sat,     p%n)
         call check_key(t, 'hydraulics.wood_pi0',             p%hydro(1:p%n)%wood_curve%pi0,           p%n)
         call check_key(t, 'hydraulics.wood_elastic_mod',     p%hydro(1:p%n)%wood_curve%elastic_mod,   p%n)
         call check_key(t, 'hydraulics.wood_apoplast_frac',   p%hydro(1:p%n)%wood_curve%apoplast_frac, p%n)
         call check_key(t, 'hydraulics.wood_water_sat',       p%hydro(1:p%n)%wood_curve%water_sat,     p%n)
         call check_key(t, 'hydraulics.wood_psi50',           p%hydro(1:p%n)%wood_psi50,               p%n)
         call check_key(t, 'hydraulics.wood_kexp',            p%hydro(1:p%n)%wood_kexp,                p%n)
         call check_key(t, 'hydraulics.wood_kmax',            p%hydro(1:p%n)%wood_kmax,                p%n)
         call check_key(t, 'hydraulics.vessel_curl',          p%hydro(1:p%n)%vessel_curl,              p%n)
         call check_key(t, 'hydraulics.root_beta',            p%hydro(1:p%n)%root_beta,                p%n)
         call check_key(t, 'hydraulics.root_depth_b1',        p%hydro(1:p%n)%root_depth_b1,            p%n)
         call check_key(t, 'hydraulics.root_depth_b2',        p%hydro(1:p%n)%root_depth_b2,            p%n)
         call check_key(t, 'hydraulics.specific_root_length', p%hydro(1:p%n)%specific_root_length,     p%n)
         call check_key(t, 'hydraulics.fine_root_radius',     p%hydro(1:p%n)%fine_root_radius,         p%n)
         !---~---
         !   The two-stream table: a leaf's single-scatter albedo is its reflectance plus its
         ! transmittance in each shortwave band and 1 - emissivity in the thermal one, and its
         ! clumping is its own, for every shipped PFT.
         !---~---
         do ip = 1_ik, p%n
            call check_close(p%optics%omega_leaf(RAD_VIS, ip), p%leaf_reflect_vis(ip) + p%leaf_transmit_vis(ip),   &
                             1.0e-12_wp, 'a shipped PFT''s VIS optics did not reach the two-stream table')
            call check_close(p%optics%omega_leaf(RAD_NIR, ip), p%leaf_reflect_nir(ip) + p%leaf_transmit_nir(ip),   &
                             1.0e-12_wp, 'a shipped PFT''s NIR optics did not reach the two-stream table')
            call check_close(p%optics%omega_leaf(RAD_LW, ip), 1.0_wp - p%leaf_emissivity(ip), 1.0e-12_wp,          &
                             'a shipped PFT''s emissivity did not reach the two-stream table')
            call check(p%optics%clumping_leaf(ip) == p%leaf_clumping(ip) .and.                                   &
                       p%optics%clumping_wood(ip) == p%wood_clumping(ip),                                       &
                       'a shipped PFT''s clumping did not reach the two-stream table')
         end do
      end associate
   end subroutine test_shipped_pft_file

   !---~---
   !   A per-PFT key of the PFT file against what the loader put in the configuration: a list
   ! gives each PFT its own value, one value gives every PFT the same.
   !---~---
   subroutine check_key(t, key, loaded, npft)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      real(wp),           intent(in) :: loaded(:)
      integer(ik),        intent(in) :: npft
      real(wp)    :: buf(64)
      integer(ik) :: n
      logical     :: same
      if (toml_is_array(t, key)) then
         call toml_real_array(t, key, buf, n)
         same = n == npft
         if (same) same = all(loaded(1:npft) == buf(1:npft))
      else
         same = toml_has(t, key)
         if (same) same = all(loaded(1:npft) == toml_real(t, key, 0.0_wp))
      end if
      call check(same, 'ed_config_pft.toml '//key//' did not reach every PFT')
   end subroutine check_key

   pure integer(ik) function count_fields(line) result(n)
      character(len=*), intent(in) :: line
      integer(ik) :: i
      n = 1_ik
      do i = 1_ik, len_trim(line)
         if (line(i:i) == ',') n = n + 1_ik
      end do
   end function count_fields

end program test_pft_config
