! xtal_calc_other_parameters.F -*-f90-*-
! Time-stamp: <2013-03-16 14:02:28 t-nissie>
! Author: Takeshi NISHIMATSU
! Reference: International Tables for Crystallography Vol.A section 5.1
!!
#include "define.h"
subroutine xtal_calc_other_parameters(xl)
  use Xtal_module
  implicit none
  type(Xtal_type), intent(inout) :: xl
  real*8 check(3,3)
  integer j

  xl%r_prim(1,:) = cross_product(xl%prim(:,2),xl%prim(:,3))
  xl%volume      = dot_product(xl%prim(:,1),xl%r_prim(1,:))
  xl%r_prim(1,:) = xl%r_prim(1,:) / xl%volume

  xl%r_prim(2,:) = cross_product(xl%prim(:,3),xl%prim(:,1))
  xl%r_prim(2,:) = xl%r_prim(2,:) / xl%volume

  xl%r_prim(3,:) = cross_product(xl%prim(:,1),xl%prim(:,2))
  xl%r_prim(3,:) = xl%r_prim(3,:) / xl%volume

  write(UNIT_LOG,'(a)') 'r_prim = '
  write(UNIT_LOG,'(3f22.16)') (xl%r_prim(j,:),j=1,3)
  write(UNIT_LOG,'(a, f22.16)') 'volume =', xl%volume

  check = matmul(xl%r_prim, xl%prim)
  write(UNIT_LOG,'(a)') 'check = '
  write(UNIT_LOG,'(3f22.16)') (check(j,:),j=1,3)

  xl%N_in_super_cell     = xl%L(1) * xl%L(2) * xl%L(3) * xl%n_atoms
  write(UNIT_LOG,'(a,i10)') 'N_in_super_cell = ', xl%N_in_super_cell

  xl%N_in_super_cell_inv = 1.0d0 / xl%N_in_super_cell
  write(UNIT_LOG,'(a,e22.16)') 'N_in_super_cell_inv = ', xl%N_in_super_cell_inv
end subroutine xtal_calc_other_parameters
