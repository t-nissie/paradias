! ewald_functions_module.F -*-f90-*-
! Time-stamp: <2013-03-16 12:39:02 t-nissie>
! Author: Takeshi NISHIMATSU
!!
#include "define.h"
module Dipole_Dipole_Functions_module
  implicit none
contains
  function ewald_kappa(epsilon_inf_det, volume, L)
    implicit none
    real*8                 ewald_kappa
    real*8,  intent(in) :: epsilon_inf_det
    real*8,  intent(in) :: volume
    integer, intent(in) :: L(3)
    integer, parameter  :: dim=3
    ! kappa for Ewald sum should be a function of epsilon_inf or epsilon_inf_det
    ewald_kappa = sqrt(M_PI) *  epsilon_inf_det ** (0.5d0/dim) &
         &                   /  volume          ** (1.0d0/dim) &
         &                   / (L(1)*L(2)*L(3)) ** (0.5d0/dim)
  end function ewald_kappa

  function dyad(a,b)   ! returns dyadic tensor
    implicit none
    real*8                dyad(3,3)
    real*8, intent(in) :: a(3), b(3)
    integer i,j
    do j = 1, 3
       do i = 1, 3
          dyad(i,j) = a(i)*b(j)
       end do
    end do
  end function dyad

  function determinant3x3(A)
    implicit none
    real*8             :: determinant3x3
    real*8, intent(in) :: A(3,3)
    determinant3x3 = &
         &  A(1,1)*A(2,2)*A(3,3) + &
         &  A(1,2)*A(2,3)*A(3,1) + &
         &  A(1,3)*A(2,1)*A(3,2) - &
         &  A(3,1)*A(2,2)*A(1,3) - &
         &  A(3,2)*A(2,3)*A(1,1) - &
         &  A(3,3)*A(2,1)*A(1,2)
  end function determinant3x3

  function inv3x3(A)
    implicit none
    real*8             :: inv3x3(3,3)
    real*8, intent(in) ::      A(3,3)
    integer info
    !real*8             ::  check(3,3)
    !integer j

    inv3x3(:,:) = A(:,:)   ! copy

    call dpotrf('U', 3, inv3x3(:,:), 3, info)   ! LAPACK
    call dpotri('U', 3, inv3x3(:,:), 3, info)   ! LAPACK

    inv3x3(2,1) = inv3x3(1,2)
    inv3x3(3,1) = inv3x3(1,3)   ! copy upper to lower
    inv3x3(3,2) = inv3x3(2,3)

    !write(UNIT_LOG,'(a)') 'A_inv = '
    !write(UNIT_LOG,'(3f22.16)') (inv3x3(j,:),j=1,3)
    !check = matmul(inv3x3(:,:), A(:,:))
    !write(UNIT_LOG,'(a)') 'check_inv3x3 = '
    !write(UNIT_LOG,'(3f22.16)') (check(j,:),j=1,3)
  end function inv3x3

  function decay1(y,yy)
    implicit none
    real*8                decay1
    real*8, intent(in) :: y,yy
    real*8 derfc
    decay1  =        derfc(y)/(y*yy) + M_2_SQRTPI * exp(-y**2) / yy
  end function decay1

  function decay2(y,yy)
    implicit none
    real*8                decay2
    real*8, intent(in) :: y,yy
    real*8 derfc
    decay2  =  3.0d0*derfc(y)/(y*yy) + M_2_SQRTPI * exp(-y**2) * (3.0d0/yy+2.0d0)
  end function decay2

  function B_decay(r,rr,kappa)
    implicit none
    real*8                B_decay
    real*8, intent(in) :: r,rr,kappa
    real*8 derfc
    B_decay  =  derfc(kappa*r)/(r*rr) + M_2_SQRTPI * kappa * exp(-kappa**2*rr) / rr
  end function B_decay

  function Crr(r,rr,kappa)
    implicit none
    real*8                Crr
    real*8, intent(in) :: r,rr,kappa
    real*8 derfc
    Crr=3*derfc(kappa*r)/(r*rr) + M_2_SQRTPI * kappa * exp(-kappa**2*rr) * (2*kappa**2+3/rr)
  end function Crr

  function C_decay(r,rr,kappa)
    implicit none
    real*8                C_decay
    real*8, intent(in) :: r,rr,kappa
    C_decay=Crr(r,rr,kappa)/rr
  end function C_decay
end module Dipole_Dipole_Functions_module
