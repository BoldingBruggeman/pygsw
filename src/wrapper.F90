module pygsw_
   use iso_c_binding, only: c_double, c_double_complex, c_int, c_char, C_NULL_CHAR, c_f_pointer, c_loc, c_ptr
   use gsw_mod_toolbox, only : gsw_SA_from_SP, gsw_pt0_from_t, gsw_ct_from_pt, gsw_rho_alpha_beta, gsw_grav
   use gsw_mod_teos10_constants, only: db2pa

   implicit none

contains

   subroutine cgws_SA_from_SP(n, sp, p, long, lat, sa) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sp(n), p(n), long(n), lat(n)
      real(c_double), intent(inout)      :: sa(n)
      sa(:) = gsw_sa_from_sp(sp(:), p(:), long(:), lat(:))
   end subroutine

   subroutine cgws_pt0_from_t(n, sa, t, p, pt) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), t(n), p(n)
      real(c_double), intent(inout)      :: pt(n)
      pt(:) = gsw_pt0_from_t(sa(:), t(:), p(:))
   end subroutine

   subroutine cgws_ct_from_pt(n, sa, pt, ct) bind(c)
      integer(c_int), intent(in), value  :: n
      real(c_double), intent(in)         :: sa(n), pt(n)
      real(c_double), intent(inout)      :: ct(n)
      ct(:) = gsw_ct_from_pt(sa, pt)
   end subroutine

   subroutine cgws_nsquared_3d(nx, ny, nz, h, sa, ct, p, lat, n2) bind(c)
      ! This routine is essentially a 3D version of 1D (k-only) gsw_nsquared
      integer(c_int), intent(in), value :: nx, ny, nz
      real(c_double), intent(in) :: h(nx, ny, nz), sa(nx, ny, nz), ct(nx, ny, nz), p(nx, ny, nz), lat(nx, ny)
      real(c_double), intent(inout) :: n2(nx, ny, nz-1)

      real(c_double), allocatable, dimension(:,:) :: dsa, sa_mid, dct, ct_mid, dp, p_mid, h_sum_inv
      real(c_double), allocatable, dimension(:,:) :: rho_mid, alpha_mid, beta_mid, grav_mid
      integer :: k

      allocate (dsa(nx, ny), sa_mid(nx, ny), dct(nx, ny), ct_mid(nx, ny), dp(nx, ny), p_mid(nx, ny), h_sum_inv(nx, ny))
      allocate (rho_mid(nx, ny), alpha_mid(nx, ny), beta_mid(nx, ny), grav_mid(nx, ny))

      do k=1,nz-1
         h_sum_inv = 1._c_double / (h(:,:,k) + h(:,:,k+1))
         grav_mid(:,:) = h_sum_inv * (h(:,:,k+1) * gsw_grav(lat(:,:), p(:,:,k)) + h(:,:,k) * gsw_grav(lat(:,:), p(:,:,k+1)))
         dsa(:,:) = (sa(:,:,k+1) - sa(:,:,k))
         sa_mid(:,:) = h_sum_inv * (h(:,:,k+1) * sa(:,:,k) + h(:,:,k) * sa(:,:,k+1))
         dct(:,:) = (ct(:,:,k+1) - ct(:,:,k))
         ct_mid(:,:) = h_sum_inv * (h(:,:,k+1) * ct(:,:,k) + h(:,:,k) * ct(:,:,k+1))
         dp(:,:) = (p(:,:,k+1) - p(:,:,k))
         p_mid(:,:) = h_sum_inv * (h(:,:,k+1) * p(:,:,k) + h(:,:,k) * p(:,:,k+1))
         call gsw_rho_alpha_beta(sa_mid,ct_mid,p_mid,rho_mid,alpha_mid,beta_mid)
         n2(:,:,k) = (grav_mid**2) * (rho_mid/(db2pa*dp)) * &
                      (beta_mid*dsa - alpha_mid*dct)
      end do
   end subroutine

end module pygsw_