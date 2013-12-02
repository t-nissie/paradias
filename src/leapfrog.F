! leapfrog.F   -*-f90-*-
! Time-stamp: <2013-03-31 22:36:24 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
subroutine leapfrog(xl)
  use Xtal_module
  use average_module
  use energy_module
  implicit none
  type(Xtal_type), intent(inout) :: xl
  type(average_type)             :: avg
  type(energy_type)              :: e
  integer i_step
  integer ix,iy,iz
  real*8 kinetic_energy, u_tmp(3)

  call update_torque(   xl%c,xl%n_atoms,xl%L,xl%dt_2)
  do i_step = 1, xl%n_thermalize + xl%n_average
     call rotate_spin(  xl%c,xl%n_atoms,xl%L,xl%dt)
     call d_V(xl,e%quadratic)
     call update_torque(xl%c,xl%n_atoms,xl%L,xl%dt)
     e%kinetic = kinetic_energy(xl) * xl%N_in_super_cell_inv

     u_tmp(:) = 0.0d0
!$omp parallel do private(iz, iy, ix) reduction(+:u_tmp)
     do iz = 0, xl%L(3)-1
        do iy = 0, xl%L(2)-1
           do ix = 0, xl%L(1)-1
              u_tmp(:) =  u_tmp(:) + xl%c%displacements_R(ix,iy,iz,:,1)  !!! SINGLE ATOM !!!
           end do
        end do
     end do
!$omp end parallel do
     e%externalH = - dot_product(xl%external_H_field(:),u_tmp(:)) * xl%N_in_super_cell_inv

     write(UNIT_LOG,'(a,3f13.8)') 'energy:', e%kinetic, e%quadratic, e%externalH

     if (i_step>xl%n_thermalize) call add_to_average(u_tmp,e,avg,xl%L)
  end do
  call divide_and_write_average(xl,avg)
end subroutine leapfrog
