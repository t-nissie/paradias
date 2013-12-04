! quadratic_module.F -*-f90-*-
! Time-stamp: <2013-12-04 22:55:39 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module quadratic_module
  use, intrinsic :: iso_c_binding
  implicit none
  type quadratic_type
     integer n_max(3)
     integer G_max(3)
     real*8  elevation   ! 2*elevation will be added to diagonal terms
     type(C_PTR)                        ::      C_bar_DD_ptr
     complex(C_DOUBLE_COMPLEX), pointer ::      C_bar_DD(:,:, :,:, :,:,:), C_bar_DD_tmp(:)

     type(C_PTR)                        :: C_tilde_short_ptr
     complex(C_DOUBLE_COMPLEX), pointer :: C_tilde_short(:,:, :,:, :,:,:), C_tilde_short_tmp(:)

     type(C_PTR)                        ::       C_tilde_ptr
     complex(C_DOUBLE_COMPLEX), pointer ::       C_tilde(:,:, :,:, :,:,:), C_tilde_tmp(:)
  end type quadratic_type

contains
  subroutine dynamical_matrix_allocate(n_atoms, Lx, Ly, Lz, qd)
    use, intrinsic :: iso_c_binding
    implicit none
    integer(C_INTPTR_T), intent(in)     :: n_atoms
    integer(C_INTPTR_T), intent(in)     :: Lx, Ly, Lz
    type(quadratic_type), intent(inout) :: qd
    integer(C_INTPTR_T)                 :: alloc_local, i, L(3) !i_atom, alpha
  end subroutine dynamical_matrix_allocate
end module quadratic_module
!Local variables:
!  compile-command: "gfortran -Wall -ffree-form -c quadratic_module.F"
!End:
