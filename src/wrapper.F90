module pygsw_
   use iso_c_binding, only: c_double, c_double_complex, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr
   use gsw_mod_toolbox, only : gsw_SA_from_SP, gsw_pt0_from_t, gsw_ct_from_pt, gsw_pt_from_ct, gsw_grav !, gsw_rho, gsw_rho_alpha_beta
   use gsw_mod_teos10_constants, only: db2pa

   implicit none

   private

contains

   subroutine cgsw_SA_from_SP(n, sp, p, long, lat, sa) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sp(n), p(n), long(n), lat(n)
      real(c_double), intent(inout)      :: sa(n)
      sa(:) = gsw_sa_from_sp(sp(:), p(:), long(:), lat(:))
   end subroutine

   subroutine cgsw_pt0_from_t(n, sa, t, p, pt) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), t(n), p(n)
      real(c_double), intent(inout)      :: pt(n)
      pt(:) = gsw_pt0_from_t(sa(:), t(:), p(:))
   end subroutine

   subroutine cgsw_ct_from_pt(n, sa, pt, ct) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), pt(n)
      real(c_double), intent(inout)      :: ct(n)
      ct(:) = gsw_ct_from_pt(sa, pt)
   end subroutine

   subroutine cgsw_pt_from_ct(n, sa, ct, pt) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), ct(n)
      real(c_double), intent(inout)      :: pt(n)
      pt(:) = gsw_pt_from_ct(sa, ct)
   end subroutine

   subroutine cgsw_rho(n, sa, ct, p, rho) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), ct(n), p(n)
      real(c_double), intent(inout)      :: rho(n)
      !rho(:) = gsw_rho(sa, ct, p)
      rho(:) = 1.0_c_double / gsw_specvol(sa, ct, p)
   end subroutine

   subroutine cgsw_nsquared_3d(nx, ny, nz, mask, h, sa, ct, p, lat, n2) bind(c)
      ! This routine is essentially a 3D version of 1D (k-only) gsw_nsquared
      integer(c_int), intent(in), value :: nx, ny, nz
      integer(c_int), intent(in) :: mask(nx, ny)
      real(c_double), intent(in) :: h(nx, ny, nz), sa(nx, ny, nz), ct(nx, ny, nz), p(nx, ny, nz), lat(nx, ny)
      real(c_double), intent(inout) :: n2(nx, ny, nz-1)

      real(c_double), allocatable, dimension(:,:) :: dsa, sa_mid, dct, ct_mid, dp, p_mid, hfrac
      real(c_double), allocatable, dimension(:,:) :: rho_mid, alpha_mid, beta_mid !, grav_mid
      real(c_double), parameter :: grav_mid = 9.81_c_double
      integer :: k

      allocate (dsa(nx, ny), sa_mid(nx, ny), dct(nx, ny), ct_mid(nx, ny), dp(nx, ny), p_mid(nx, ny), step(nx, ny))
      allocate (rho_mid(nx, ny), alpha_mid(nx, ny), beta_mid(nx, ny)) !, grav_mid(nx, ny))

      do k=1,nz-1
         step = h(:,:,k) / (h(:,:,k) + h(:,:,k+1))
         !grav_mid(:,:) = (h(:,:,k+1) * gsw_grav(lat(:,:), p(:,:,k)) + h(:,:,k) * gsw_grav(lat(:,:), p(:,:,k+1))) / (h(:,:,k) + h(:,:,k+1))
         dsa(:,:) = (sa(:,:,k+1) - sa(:,:,k))
         sa_mid(:,:) = sa(:,:,k) + step * dsa
         dct(:,:) = (ct(:,:,k+1) - ct(:,:,k))
         ct_mid(:,:) = ct(:,:,k) + step * dct
         dp(:,:) = (p(:,:,k+1) - p(:,:,k))
         p_mid(:,:) = p(:,:,k) + step * dp
         call gsw_rho_alpha_beta(sa_mid, ct_mid, p_mid, rho_mid, alpha_mid, beta_mid)
         where (mask /= 0) n2(:,:,k) = (grav_mid**2) * (rho_mid / (db2pa * dp)) * &
                         (beta_mid * dsa - alpha_mid * dct)
      end do
   end subroutine

#include "../extern/GSW-Fortran/toolbox/gsw_rho_alpha_beta.f90"
#include "../extern/GSW-Fortran/toolbox/gsw_specvol.f90"

end module pygsw_