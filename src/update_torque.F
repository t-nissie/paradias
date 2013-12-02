! update_torque.F   -*-f90-*-
! Time-stamp: <2013-03-31 22:37:13 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
subroutine update_torque(c,n_atoms,L,dt)
  use coordinates_module
  implicit none
  type(coordinates_type), intent(inout) :: c
  integer,                intent(in)    :: n_atoms
  integer,                intent(in)    :: L(3)
  real*8,                 intent(in)    :: dt
  integer I_atom
  integer ix,iy,iz
!$omp parallel do
  do iz=0, L(3)-1
     do iy=0, L(2)-1
        do ix=0, L(1)-1
           do I_atom = 1, n_atoms
              c%displacements_P(ix,iy,iz, :,I_atom) = &
                   & c%displacements_P(ix,iy,iz, :,I_atom) - &
                   & dt * (  c%d_V_d_displacements_R(ix,iy,iz, :,I_atom) &
                   &             - c%displacements_R(ix,iy,iz, :,I_atom) &
                   & * dot_product(c%displacements_R(ix,iy,iz, :,I_atom),&
                   &         c%d_V_d_displacements_R(ix,iy,iz, :,I_atom))  )
              c%displacements_P(ix,iy,iz, :,I_atom) = &
                   & c%displacements_P(ix,iy,iz, :,I_atom) -  c%displacements_R(ix,iy,iz, :,I_atom) &
                   &                            * dot_product(c%displacements_R(ix,iy,iz, :,I_atom),&
                   &                                          c%displacements_P(ix,iy,iz, :,I_atom))
           end do
        end do
     end do
  end do
!$omp end parallel do
end subroutine update_torque
