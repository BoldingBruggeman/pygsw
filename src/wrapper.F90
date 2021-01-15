module pygsw_
   use iso_c_binding, only: c_double, c_double_complex, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr
   use gsw_mod_toolbox, only : gsw_SA_from_SP, gsw_pt0_from_t

   implicit none

contains

   subroutine calculate_pt(long, lat, n, z, t, sp, pt) bind(c)
!DIR$ ATTRIBUTES DLLEXPORT :: calculate_pt
      real(c_double), value, intent(in) :: long, lat
      integer(c_int), value, intent(in) :: n
      real(c_double), intent(in), target :: z, t, sp(*)
      real(c_double), intent(inout), target :: pt(*)

      real(c_double), pointer :: z_(:), t_(:), sp_(:), pt_(:)
      real(c_double) :: sa(n), p(n)

      call c_f_pointer(c_loc(z), z_, (/n/))
      call c_f_pointer(c_loc(t), t_, (/n/))
      call c_f_pointer(c_loc(sp), sp_, (/n/))
      call c_f_pointer(c_loc(pt), pt_, (/n/))

      p(:) = -z_(:) !gsw_p_from_z?
      sa(:) = gsw_SA_from_SP(sp_(:), p(:), long, lat)
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