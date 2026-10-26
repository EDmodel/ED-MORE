! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_maps_params -- the [maps] block, as a config leaf: the static maps a run reads its  !
! soil and its height allometry's site index from (docs/configuration.md, "Input maps").   !
! Each map is on the forcing grid; meds_maps reads it at the run's cells into a            !
! map_inputs_t. A region needs the three soil maps; a site may use them or set [soil] by   !
! hand, and a key a map provides is refused while the map is set.                          !
!==========================================================================================!
module meds_maps_params
   use meds_kinds,         only : wp, ik
   use meds_soil_hydr_lib, only : soil_texture_t
   implicit none
   private

   public :: maps_params_t, soil_maps_on, MAPS_PATH_LEN
   public :: map_inputs_t, MAP_BC_FREE_DRAIN, MAP_BC_LEAKY, MAP_BC_SEALED, MAP_BC_AQUIFER

   integer, parameter :: MAPS_PATH_LEN = 512

   type :: maps_params_t
      !----- The three soil maps (docs/configuration.md, "Input maps"). -------------------!
      character(len=MAPS_PATH_LEN) :: soil_depth_map     = ''
      character(len=MAPS_PATH_LEN) :: soil_texture_map   = ''
      character(len=MAPS_PATH_LEN) :: soil_chemistry_map = ''
      !----- Chave et al. (2014)'s stress index E, for height_allometry = "chave2014". ----!
      character(len=MAPS_PATH_LEN) :: chave_e_map        = ''
   end type maps_params_t

   !----- The bottom water boundary of the soil depth map's bc_class. ---------------------------!
   integer(ik), parameter :: MAP_BC_FREE_DRAIN = 1_ik   !< q = K(theta)
   integer(ik), parameter :: MAP_BC_LEAKY      = 2_ik   !< q = f_leak K(theta)
   integer(ik), parameter :: MAP_BC_SEALED     = 3_ik   !< q = 0
   integer(ik), parameter :: MAP_BC_AQUIFER    = 4_ik   !< a water table at the column base

   !----- What the [maps] maps hold at each cell of a run, read once (read_map_inputs). --!
   type :: map_inputs_t
      logical :: soil  = .false.                          !< the soil maps were read
      logical :: chave = .false.                          !< the chave_e map was read
      real(wp), allocatable :: layer_top(:), layer_bottom(:)   !< [m] the soil maps' layers, positive down
      real(wp), allocatable :: soil_depth(:)              !< [m] (cell) the soil column's depth
      real(wp), allocatable :: bc_class(:)                !< [-] (cell) its bottom water boundary, MAP_BC_*
      real(wp), allocatable :: f_leak(:)                  !< [-] (cell) a leaky base's drainage factor
      real(wp), allocatable :: deep_temp(:)               !< [K] (cell) the deep soil temperature
      type(soil_texture_t), allocatable :: texture(:,:)   !< (map layer, cell) each layer's soil
      real(wp), allocatable :: chave_e(:)                 !< [-] (cell) Chave et al. (2014)'s E
   end type map_inputs_t

contains

   !----- Whether the soil comes from the maps (the loader sets all three or none). -------!
   pure logical function soil_maps_on(s)
      type(maps_params_t), intent(in) :: s
      soil_maps_on = len_trim(s%soil_depth_map) > 0
   end function soil_maps_on

end module meds_maps_params
