! k2eigenvalues.F -*-f90-*-
! Time-stamp: <2013-03-17 13:02:06 t-nissie>
! Author: Takeshi NISHIMATSU
! TRICKY: because L(1)/2 in C_tild(0:L(1)/2,
!                                  0:L(2)-1,
!                                  0:L(3)-1, :,:)
! Parents: dispersion.F, lu_matrix.F
!!
#include "define.h"
subroutine k2eigenvalues(n, k, L, mat, eigenvalues, eigenvectors)
  implicit none
  !--- arguments -------------------------------------
  integer,    intent(in)  :: n
  real*8,     intent(in)  :: k(3)
  integer,    intent(in)  :: L(3)
  complex*16, intent(in)  :: mat(0:L(1)/2, &
       &                         0:L(2)-1, &
       &                         0:L(3)-1, n, n)
  real*8,     intent(out) :: eigenvalues(n)
  complex*16, intent(out) :: eigenvectors(n,n)
  !--- local -----------------------------------------
  real*8 k_reduced(3)
  real*8 k_plot(3)
  complex*16  work(LWORK)
  real*8     rwork(LWORK)
  integer i,info
  !character (len=50) fmt
  !--- END of variable definitions -------------------

  !write(UNIT_LOG,'(3f10.5)') k(:)

  k_reduced(:) = k(:)
  do i = 1, 3
     do while (k_reduced(i).le.-0.5d0 .or. k_reduced(i).gt.0.5d0)
        if (k_reduced(i).le.-0.5d0) k_reduced(i) = k_reduced(i) + 1.0d0
        if (k_reduced(i).gt. 0.5d0) k_reduced(i) = k_reduced(i) - 1.0d0
     end do
  end do

  if (k_reduced(1).lt.0.0d0) then
     k_plot(:) = -k_reduced(:)
  else
     k_plot(:) =  k_reduced(:)
  end if

  do i = 2, 3   !!!only 2 and 3!!!
     if (k_plot(i).lt.0.0d0)         k_plot(i) = k_plot(i) + 1.0d0
     if (k_plot(i).gt.0.999999999d0) k_plot(i) = 0.0d0
  end do

!   write(8,'(a,3i3,3f10.4)') '#',&
!        &                  NINT(k_plot(1)*L(1)),&
!        &                  NINT(k_plot(2)*L(2)),&
!        &                  NINT(k_plot(3)*L(3)), k_reduced

  eigenvectors(:,:) = mat(NINT(k_plot(1)*L(1)),&
       &                  NINT(k_plot(2)*L(2)),&
       &                  NINT(k_plot(3)*L(3)),:,:)

  if (k_reduced(1).lt.0) then
     eigenvectors(:,:) = conjg( eigenvectors(:,:) )
  end if


!  write(UNIT_LOG,'(3f7.3)') k
!  write(UNIT_LOG,'(9(f7.3f6.3))') eigenvectors(:,:)

#if defined(__xlc__)
  stop 'ESSL does not have zheev().'
#else
!LAPACK zheev(JOBZ,UPLO,N,A,           LXL, W,       WORK,LWORK,RWORK,INFO)
  call  zheev('V','L',  n,eigenvectors,n,eigenvalues,work,LWORK,rwork,info)
#endif
end subroutine k2eigenvalues
