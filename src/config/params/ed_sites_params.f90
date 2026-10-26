! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! ed_sites_params -- the [[sites]] tables, as a config leaf: one record per site of the    !
! polygon (docs/configuration.md, "Several sites"). A site's surface is built from it at   !
! setup (make_site_surface, through site_config); the record itself is never state.        !
!==========================================================================================!
module ed_sites_params
   use ed_kinds,                only : wp, ik
   use ed_soil_params, only :   soil_params_t
   implicit none
   private

   public :: site_spec_t

   !----- One site: its share of the polygon, terrain and census rows. A table that sets   !
   !      any of [soil]'s column keys gives the site its own ground, with [soil]'s value   !
   !      for each key it leaves out. -----------------------------------------------------!
   type :: site_spec_t
      real(wp)    :: area_fraction  = 0.0_wp      !< [-] its share of the polygon
      logical     :: has_elevation  = .false.     !< .false.: the polygon's elevation
      real(wp)    :: elevation      = 0.0_wp      !< [m]
      real(wp)    :: slope          = 0.0_wp      !< [degrees] from horizontal
      real(wp)    :: aspect         = 0.0_wp      !< [degrees] clockwise from north
      integer(ik) :: census_site_id = 0_ik        !< the census rows it reads; 0: not given
      logical     :: own_ground     = .false.     !< .false.: the polygon's ground
      type(soil_params_t) :: soil                 !< its ground: [soil] with the table's keys
   end type site_spec_t

end module ed_sites_params
