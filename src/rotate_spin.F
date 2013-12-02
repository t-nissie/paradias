! rotate_spin.F   -*-f90-*-
! Time-stamp: <2013-03-31 19:27:22 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
subroutine rotate_spin(c,n_atoms,L,dt)
  use coordinates_module
  implicit none
  type(coordinates_type), intent(inout) :: c
  integer,                intent(in)    :: n_atoms
  integer,                intent(in)    :: L(3)
  real*8,                 intent(in)    :: dt
  integer I_atom
  integer ix,iy,iz
  real*8 scale
!$omp parallel do
  do iz=0, L(3)-1
     do iy=0, L(2)-1
        do ix=0, L(1)-1
           do I_atom = 1, n_atoms
              c%displacements_R(ix,iy,iz, :,I_atom) = &
                   & c%displacements_R(ix,iy,iz, :,I_atom) + dt * c%displacements_P(ix,iy,iz, :,I_atom)
              scale = dot_product( c%displacements_R(ix,iy,iz, :,I_atom), &
                   &               c%displacements_R(ix,iy,iz, :,I_atom) )
              c%displacements_R(ix,iy,iz, :,I_atom) = &
                   & c%displacements_R(ix,iy,iz, :,I_atom) / sqrt(scale)
           end do
        end do
     end do
  end do
!$omp end parallel do
end subroutine rotate_spin
