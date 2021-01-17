module pygsw_
   use iso_c_binding, only: c_double, c_double_complex, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr
   use gsw_mod_toolbox, only : gsw_SA_from_SP, gsw_pt0_from_t

   implicit none

contains

   subroutine calculate_pt(n, long, lat, z, t, sp, pt) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: calculate_pt
      integer(c_int), intent(in),    value  :: n
      real(c_double), intent(in),    target :: long(*), lat(*), z(*), t(*), sp(*)
      real(c_double), intent(inout), target :: pt(*)

      real(c_double), pointer :: long_(:), lat_(:), z_(:), t_(:), sp_(:), pt_(:)
      real(c_double), allocatable :: sa(:), p(:)

      allocate(sa(n), p(n))
      call c_f_pointer(c_loc(long), long_, (/n/))
      call c_f_pointer(c_loc(lat), lat_, (/n/))
      call c_f_pointer(c_loc(z), z_, (/n/))
      call c_f_pointer(c_loc(t), t_, (/n/))
      call c_f_pointer(c_loc(sp), sp_, (/n/))
      call c_f_pointer(c_loc(pt), pt_, (/n/))

      p(:) = -z_(:) !gsw_p_from_z?
      sa(:) = gsw_SA_from_SP(sp_(:), p(:), long_(:), lat_(:))
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]

      pt_(:) = gsw_pt0_from_t(sa(:), t_(:), p(:))
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
   end subroutine

end module pygsw_