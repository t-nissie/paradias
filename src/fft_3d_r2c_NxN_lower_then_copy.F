! fft_3d_r2c_NxN_lower_then_copy.F -*-f90-*-
! Time-stamp: <2013-02-14 11:11:13 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"

subroutine fft_3d_r2c_NxN_lower_then_copy(Lx,Ly,Lz,N,r,c)
  implicit none
  integer    :: Lx,Ly,Lz,N
  real*8     :: r(0:Lx-1, 0:Ly-1, 0:Lz-1, N, N)
  complex*16 :: c(0:Lx/2, 0:Ly-1, 0:Lz-1, N, N)
  integer    :: alpha, beta
#if defined(HAVE_LIBFFTW3)
  integer ip
  integer*8 plan(N*(N+1)/2)
# include "fftw3.f"
  ip = 0
  do beta = 1, N
     do alpha = 1, beta
        ip = ip + 1
        call dfftw_plan_dft_r2c_3d(plan(ip), Lx, Ly, Lz, &
             & r(0,0,0,alpha,beta), c(0,0,0,alpha,beta), FFTW_ESTIMATE)
     end do
  end do
  do ip = 1, N*(N+1)/2
     call dfftw_execute(plan(ip))
  end do
  do ip = 1, N*(N+1)/2
     call dfftw_destroy_plan(plan(ip))
  end do
#elif defined(HAVE_LIBACML_MP)
  integer info
  real*8 comm(4*Lx+6*Ly+6*Lz+500)
  do beta = 1, 3
     do alpha = 1, beta
        call DZFFT3D(2, Lx, Ly, Lz, &
             & r(0,0,0,alpha,beta), c(0,0,0,alpha,beta), comm, info)
     end do
  end do
#else
  do beta = 1, N
     do alpha = 1, beta
        call fft_3d_r2c_or_c2r_out_of_place(.true., Lx, Ly, Lz, &
             & r(0,0,0,alpha,beta), c(0,0,0,alpha,beta))
     end do
  end do
#endif

  ! Copy triangle
  do beta = 1, N
     do alpha = 1, beta-1
        c(:,:,:,beta,alpha) = conjg(c(:,:,:,alpha,beta))
     end do
  end do
end subroutine fft_3d_r2c_NxN_lower_then_copy
