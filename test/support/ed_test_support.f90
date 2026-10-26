! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_test_support -- the in-code test CONFIGURATION fixture, plus a re-export of the        !
! assertion helpers so a test needs one `use` rather than two. The assertions themselves live !
! in ed_test_assert, which depends on nothing but a kind; see that module for why, and for    !
! the two assertion families (fatal and accumulating).                                     !
!==========================================================================================!
module ed_test_support
   use ed_kinds,        only : wp, ik
   use ed_constants,    only : day_sec
   use ed_time,         only : ed_time_t
   use ed_pft_params,   only : pft_table_t, PATH_C3, PATH_C4
   use ed_config,       only : ed_config_t, derive_parameters, INIT_BARE, alloc_pft_table
   use ed_fast_params,  only : INTEG_ESDIRK2
   use ed_leaf_params,     only : SM_MEDLYN, COLIM_QUADRATIC
   use ed_therm_lib,       only : TEMP_RESPONSE_PEAKED
   use ed_test_assert,   only : check, check_close, check_true, check_int, banner,             &
                                test_report, test_reset, test_failures
   implicit none
   private

   !----- Re-exported, so a test that links the full support library needs one `use`. -----!
   public :: check, check_close, check_true, check_int, banner
   public :: test_report, test_reset, test_failures
   public :: build_test_config, test_pft_table

contains

   !---------------------------------------------------------------------------------------!
   ! The tests' configuration FIXTURE: a complete configuration built in code, so that     !
   ! no test reads a TOML file. Every value it sets is the shipped config's                !
   ! (ed_config_main.toml, ed_config_pft.toml), dt_fast (900 s) and                        !
   ! growth_memory_days (365 days) among them, except on purpose:                          !
   !   * no files: the start is bare ground from no census, and nothing is written;        !
   !   * allometry and hydraulics differ by PFT. Where the PFT file gives one value        !
   !     for every PFT, PFT 1 keeps it and PFTs 2 and 3 step away from it (shorter         !
   !     and stouter, with less sapwood, tougher leaves, safer xylem and deeper            !
   !     roots), so a check that reads another PFT's record sees another number;           !
   !     where the file lists a value per PFT (wood_kmax), the fixture takes it.           !
   !     wood_kexp is 2, 1 and 3, so the Kirchhoff integral's two closed forms and         !
   !     its lookup table each have a PFT.                                                 !
   ! The blocks it does not set ([soil], [snow], [aerodynamics], [soil_carbon], the        !
   ! solver settings of [fast]) keep their types' initializers, which are the shipped      !
   ! values. A test that needs another value (a shorter dt_fast for an integrator          !
   ! scenario) sets it itself and says why. Pass dt_slow [s] to vary the slow step.        !
   !---------------------------------------------------------------------------------------!
   function build_test_config(dt_slow) result(cfg)
      real(wp), intent(in), optional :: dt_slow
      type(ed_config_t)              :: cfg

      cfg%dt_slow = day_sec
      if (present(dt_slow)) cfg%dt_slow = dt_slow
      cfg%fast_biophysics_on = .false.
      cfg%dt_fast            = 900.0_wp
      cfg%time_integrator    = INTEG_ESDIRK2
      cfg%start_time = ed_time_t(2000_ik,   1_ik, 1_ik)
      cfg%end_time   = ed_time_t(2250_ik,   1_ik, 1_ik)
      cfg%demography%demography_on = .true.
      cfg%demography%recruitment_on = .true.
      cfg%demography%cohort_fusion_on = .true.
      cfg%demography%cohort_fission_on = .true.
      cfg%demography%patch_disturbance_on = .true.
      cfg%demography%patch_fusion_on = .true.

      cfg%demography%max_cohort = 60_ik
      cfg%demography%n_cohort_fusion_iter = 6_ik
      cfg%demography%cohort_size_tol_min = 0.02_wp
      cfg%demography%cohort_size_tol_max = 0.10_wp
      cfg%demography%cohort_lai_cap = 1.0_wp
      cfg%demography%min_cohort_agb = 1.0e-6_wp
      cfg%demography%negligible_nplant = 1.0e-8_wp
      cfg%demography%split_eps = 1.0e-4_wp

      cfg%demography%n_height_layers = 16_ik
      cfg%demography%max_patch = 12_ik
      cfg%demography%n_patch_fusion_iter = 6_ik
      cfg%demography%patch_light_tol = 0.10_wp
      cfg%demography%patch_light_tol_max = 0.15_wp
      cfg%demography%patch_light_maxdev_factor = 1.5_wp
      cfg%demography%min_patch_area = 1.0e-4_wp

      cfg%demography%patch_disturbance_rate = 0.014_wp
      cfg%demography%disturbance_survive_height = 10.0_wp
      cfg%demography%growth_memory_days = 365.0_wp
      cfg%demography%min_recruit_size = 1.0e-2_wp
      cfg%demography%conservation_tol = 0.001_wp
      cfg%demography%light_ext_coeff_fusefiss = 0.5_wp

      cfg%init%init_mode    = INIT_BARE
      cfg%init%restart_file = ''
      cfg%init%census_file  = ''
      cfg%init%pft_config   = ''
      cfg%output%state_file_dir    = '.'
      cfg%output%state_file_prefix = 'test'
      cfg%output%state_write       = .false.
      cfg%gpp_ref = 0.3_wp

      !----- Leaf physiology: the model selectors and the biochemistry the PFTs share. ----!
      associate (lp => cfg%leaf_physiology)
         lp%stomatal_model     = SM_MEDLYN
         lp%temp_response_form = TEMP_RESPONSE_PEAKED
         lp%colimitation       = COLIM_QUADRATIC
         lp%use_boundary_layer = .false.
         lp%kc25    = 40.49_wp
         lp%ko25    = 27840.0_wp
         lp%gstar25 = 4.275_wp
         lp%ea_kc    = 79430.0_wp
         lp%ea_ko    = 36380.0_wp
         lp%ea_gstar = 37830.0_wp
         lp%ea_vcmax = 65330.0_wp
         lp%ea_jmax  = 43540.0_wp
         lp%ea_rd    = 46390.0_wp
         lp%hd_vcmax = 200000.0_wp
         lp%hd_jmax  = 200000.0_wp
         lp%hd_rd    = 200000.0_wp
         lp%ds_vcmax = 650.0_wp
         lp%ds_jmax  = 640.0_wp
         lp%ds_rd    = 490.0_wp
         lp%o2_mol_frac      = 0.209_wp
         lp%leaf_absorptance = 0.85_wp
         lp%phi_psii         = 0.74_wp
      end associate

      call alloc_pft_table(cfg%pft, 3_ik)
      associate (p => cfg%pft)
         p%wood_density = [ 0.40_wp, 0.60_wp, 0.85_wp ]
         p%reproduction_investment_fraction = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%recruit_carbon_efficiency = [ 1.0e-3_wp, 1.0e-3_wp, 1.0e-3_wp ]
         p%seed_rain_recruits = [ 0.01_wp, 0.01_wp, 0.01_wp ]
         p%include_pft = [ 1_ik, 1_ik, 1_ik ]
         p%min_cohort_height = 2.0_wp
         p%min_repro_height  = [ 13.912064431085485_wp, 13.912064431085485_wp, 13.912064431085485_wp ]
         !----- The mortality hazard and the plasticity slopes. ---------------------------!
         p%mort_gamma = [ 0.019815037439806053_wp, 0.0094_wp, 0.004953547651645506_wp ]
         p%mort_alpha = [ 0.07968043353213256_wp, 0.05_wp, 0.0335056513529311_wp ]
         p%mort_beta  = [ 16.71633880800162_wp, 18.72_wp, 20.631926375047502_wp ]
         p%kplastic_vcmax = [ -0.1459341188305866_wp, -0.11760191025318557_wp, -0.10943768264434321_wp ]
         p%kplastic_rd    = p%kplastic_vcmax
         p%kplastic_sla   = [ 0.15654750929099287_wp, 0.12615477645341722_wp, 0.11739678683665906_wp ]
         p%kplastic_leaf_turnover = [ -0.058535787713143994_wp, -0.015560662518427415_wp, 0.0_wp ]
         !----- Leaf-photosynthesis per-PFT traits (PFT 3 is C4). -------------------------!
         p%photosynthetic_pathway = [ PATH_C3, PATH_C3, PATH_C4 ]
         p%vcmax25          = [ 60.0_wp, 45.0_wp, 40.0_wp ]
         p%jmax_vcmax_ratio = [ 1.8_wp, 1.7_wp, 4.0_wp ]
         p%tpu_vcmax_ratio  = [ 0.167_wp, 0.167_wp, 0.167_wp ]
         p%rd_vcmax_ratio   = [ 0.015_wp, 0.015_wp, 0.025_wp ]
         p%kp25             = [ 1.0e-9_wp, 1.0e-9_wp, 0.7_wp ]
         p%stomatal_g0      = [ 0.01_wp, 0.01_wp, 0.04_wp ]
         p%stomatal_g1      = [ 4.0_wp, 3.0_wp, 1.6_wp ]
         p%stomatal_d0      = [ 1500.0_wp, 1500.0_wp, 1500.0_wp ]
         p%quantum_yield_c4 = [ 0.0_wp, 0.0_wp, 0.04_wp ]
         p%theta_j          = [ 0.70_wp, 0.70_wp, 0.85_wp ]
         p%theta_cj_c3      = [ 0.98_wp, 0.98_wp, 0.98_wp ]
         p%theta_ip_c3      = [ 0.95_wp, 0.95_wp, 0.95_wp ]
         p%theta_cj_c4      = [ 0.80_wp, 0.80_wp, 0.80_wp ]
         p%theta_ic_c4      = [ 0.95_wp, 0.95_wp, 0.95_wp ]
         p%katul_lambda25   = [ 600.0_wp, 800.0_wp, 350.0_wp ]
         p%wstress_lambda_exp = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%wstress_sref_stomata = [ 2.0_wp, 2.0_wp, 2.0_wp ]
         !----- Carbon-dynamics traits. ---------------------------------------------------!
         p%sla                    = [ 16.0_wp, 13.0_wp, 10.0_wp ]
         !----- Canopy-element geometry. --------------------------------------------------!
         p%leaf_width             = [ 0.04_wp, 0.04_wp, 0.04_wp ]
         p%branch_diameter        = [ 0.02_wp, 0.02_wp, 0.02_wp ]
         p%root_to_leaf_ratio     = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         !----- Non-leaf maintenance respiration. -----------------------------------------!
         p%is_woody               = [ .true., .true., .true. ]
         p%stem_resp_factor25     = [ 0.06_wp, 0.06_wp, 0.06_wp ]
         p%root_resp_factor25     = [ 0.30_wp, 0.30_wp, 0.30_wp ]
         !----- Canopy optics. ------------------------------------------------------------!
         p%leaf_reflect_vis  = [ 0.10_wp, 0.10_wp, 0.10_wp ]
         p%leaf_transmit_vis = [ 0.05_wp, 0.05_wp, 0.05_wp ]
         p%leaf_reflect_nir  = [ 0.45_wp, 0.45_wp, 0.45_wp ]
         p%leaf_transmit_nir = [ 0.25_wp, 0.25_wp, 0.25_wp ]
         p%leaf_emissivity   = [ 0.97_wp, 0.97_wp, 0.97_wp ]
         p%wood_reflect_vis  = [ 0.11_wp, 0.11_wp, 0.11_wp ]
         p%wood_transmit_vis = [ 0.001_wp, 0.001_wp, 0.001_wp ]
         p%wood_reflect_nir  = [ 0.25_wp, 0.25_wp, 0.25_wp ]
         p%wood_transmit_nir = [ 0.001_wp, 0.001_wp, 0.001_wp ]
         p%wood_emissivity   = [ 0.90_wp, 0.90_wp, 0.90_wp ]
         p%leaf_clumping     = [ 0.80_wp, 0.80_wp, 0.80_wp ]
         p%wood_clumping     = [ 0.50_wp, 0.50_wp, 0.50_wp ]
         p%leaf_angle_mean   = [ 45.0_wp, 45.0_wp, 45.0_wp ]
         p%leaf_angle_std    = [ 20.0_wp, 20.0_wp, 20.0_wp ]
         p%storage_cushion        = [ 3.0_wp, 3.0_wp, 3.0_wp ]
         p%growth_resp_factor     = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%leaf_turnover_rate_toc = [ 1.0_wp, 0.5_wp, 0.33_wp ]          ! [1/yr]
         p%fineroot_turnover_rate = [ 1.0_wp, 0.8_wp, 0.6_wp ]
         p%f_labile_leaf          = [ 0.80_wp, 0.70_wp, 0.60_wp ]
         p%f_labile_stem          = [ 0.05_wp, 0.05_wp, 0.05_wp ]
         p%struct_lignin_frac     = [ 0.20_wp, 0.25_wp, 0.30_wp ]
      end associate

      !----- Each PFT's allometry, all on the power-law height curve: PFT 1 on the         !
      !      shipped values, PFTs 2 and 3 shorter, stouter and with less sapwood. ---------!
      associate (a => cfg%pft%allom)
         a%hgt_max          = [ 46.0_wp, 42.0_wp, 38.0_wp ]
         a%b1Ht             = [ 1.139963_wp, 1.09_wp, 1.04_wp ]
         a%b2Ht             = [ 0.564899_wp, 0.570_wp, 0.575_wp ]
         a%agb_c1           = [ 0.03365_wp, 0.0320_wp, 0.0305_wp ]
         a%agb_c2           = [ 0.976_wp, 0.970_wp, 0.964_wp ]
         a%lai_b1           = [ 0.23384770_wp, 0.220_wp, 0.205_wp ]
         a%lai_b2           = [ 0.6410495_wp, 0.645_wp, 0.650_wp ]
         a%wai_b1           = [ 0.0096_wp, 0.0090_wp, 0.0084_wp ]
         a%wai_b2           = [ 2.0947_wp, 2.110_wp, 2.125_wp ]
         a%sapwood_area_b1  = [ 1.582_wp, 1.450_wp, 1.320_wp ]
         a%sapwood_area_b2  = [ 1.764_wp, 1.780_wp, 1.800_wp ]
         a%aboveground_frac = [ 0.70_wp, 0.72_wp, 0.74_wp ]
      end associate
      !----- ... and its hydraulics: PFT 1 on the shipped values, PFTs 2 and 3 with        !
      !      tougher leaves (a lower turgor-loss point), safer xylem and deeper roots;     !
      !      wood_kmax is the PFT file's own per-PFT fit to wood density. -----------------!
      associate (h => cfg%pft%hydro)
         h%leaf_curve%pi0           = [ -1.5_wp, -1.7_wp, -1.9_wp ]
         h%leaf_curve%elastic_mod   = [ 12.0_wp, 14.0_wp, 16.0_wp ]
         h%leaf_curve%apoplast_frac = [ 0.30_wp, 0.32_wp, 0.34_wp ]
         h%leaf_curve%water_sat     = [ 2.0_wp, 1.85_wp, 1.7_wp ]
         h%wood_curve%pi0           = [ -1.0_wp, -1.15_wp, -1.3_wp ]
         h%wood_curve%elastic_mod   = [ 8.0_wp, 10.0_wp, 12.0_wp ]
         h%wood_curve%apoplast_frac = [ 0.20_wp, 0.22_wp, 0.24_wp ]
         h%wood_curve%water_sat     = [ 1.0_wp, 0.85_wp, 0.70_wp ]
         h%wood_psi50           = [ -2.0_wp, -2.5_wp, -3.0_wp ]
         h%wood_kexp            = [ 2.0_wp, 1.0_wp, 3.0_wp ]
         h%wood_kmax            = [ 5.340397027887704_wp, 3.2683981418222547_wp, 1.7692398658114543_wp ]
         h%vessel_curl          = [ 1.5_wp, 1.6_wp, 1.7_wp ]
         h%root_beta            = [ 0.10_wp, 0.15_wp, 0.20_wp ]
         h%root_depth_b1        = [ 1.1140580_wp, 1.20_wp, 1.30_wp ]
         h%root_depth_b2        = [ 0.4223014_wp, 0.43_wp, 0.44_wp ]
         h%specific_root_length = [ 2.0e4_wp, 1.7e4_wp, 1.4e4_wp ]
         h%fine_root_radius     = [ 2.5e-4_wp, 2.75e-4_wp, 3.0e-4_wp ]
      end associate
      call derive_parameters(cfg)
   end function build_test_config

   !----- The test configuration's PFT table, for a fixture site that has no config of its own. --!
   function test_pft_table() result(pft)
      type(pft_table_t)   :: pft
      type(ed_config_t)   :: cfg
      cfg = build_test_config()
      pft = cfg%pft
   end function test_pft_table

end module ed_test_support
