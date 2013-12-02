! quadratic_module.F -*-f90-*-
! Time-stamp: <2013-02-14 18:04:14 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module quadratic_module
  implicit none
  type quadratic_type
     integer n_max(3)
     integer G_max(3)
     real*8  elevation   ! 2*elevation will be added to diagonal terms
     complex*16, allocatable ::      C_bar_DD(:,:,:, :,:, :,:)
     complex*16, allocatable :: C_tilde_short(:,:,:, :,:, :,:)
     complex*16, allocatable ::       C_tilde(:,:,:, :,:, :,:)
  end type quadratic_type

contains
  subroutine dynamical_matrix_allocate(n_atoms, L, qd)
    implicit none
    integer, intent(in) :: n_atoms
    integer, intent(in) :: L(3)
    type(quadratic_type), intent(inout) :: qd

    if (allocated(qd%C_bar_DD))      deallocate(qd%C_bar_DD)
    if (allocated(qd%C_tilde_short)) deallocate(qd%C_tilde_short)
    if (allocated(qd%C_tilde))       deallocate(qd%C_tilde)
    allocate(     qd%C_bar_DD(0:L(1)/2, &
         &                    0:L(2)-1, &
         &                    0:L(3)-1, &
         &                        3,n_atoms, &
         &                        3,n_atoms) )
    allocate(     qd%C_tilde_short(0:L(1)/2, &
         &                         0:L(2)-1, &
         &                         0:L(3)-1, &
         &                        3,n_atoms, &
         &                        3,n_atoms) )
    allocate(     qd%C_tilde(0:L(1)/2, &
         &                   0:L(2)-1, &
         &                   0:L(3)-1, &
         &                        3,n_atoms, &
         &                        3,n_atoms) )
  end subroutine dynamical_matrix_allocate
end module quadratic_module
