! velocity_scaling.F   -*-f90-*-
! Time-stamp: <2013-03-31 20:39:52 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
subroutine velocity_scaling(xl)
  use Xtal_module
  use average_module
  use energy_module
  implicit none
  type(Xtal_type), intent(inout) :: xl
  type(average_type)             :: avg
  type(energy_type)              :: e
  integer i_step
  integer ix,iy,iz
  real*8 scale, kinetic_energy, u_tmp(3)

  do i_step = 1, xl%n_thermalize + xl%n_average

     call rotate_spin(xl%c,xl%n_atoms,xl%L,xl%dt)

     call d_V(xl,e%quadratic)

     !write(UNIT_LOG,'(a,3f18.8)') '        d_V_d_u(1,1,1, :,1) = ', xl%c%d_V_d_displacements_R(1,1,1, :,1)

     call update_torque(xl%c,xl%n_atoms,xl%L,xl%dt)

     ! write(UNIT_LOG,'(a,3f18.8)') 'mdified d_V_d_u(1,1,1, :,1) = ', xl%c%d_V_d_displacements_R(1,1,1, :,1) &
     !           &               - xl%c%displacements_R(1,1,1, :,1) &
     !           &   * dot_product(xl%c%displacements_R(1,1,1, :,1),&
     !           &           xl%c%d_V_d_displacements_R(1,1,1, :,1))
     ! write(UNIT_LOG,'(a,3f18.8)') 'displacements_R(1,1,1, :,1) = ', xl%c%displacements_R(1,1,1, :,1)
     ! write(UNIT_LOG,'(a,3f18.8)') 'displacements_P(1,1,1, :,1) = ', xl%c%displacements_P(1,1,1, :,1)

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

     scale = sqrt( xl%temperature / e%kinetic )   !!! Only for the case of N_free = 2 !!!
!$omp parallel do
     do iz=0, xl%L(3)-1
        xl%c%displacements_P(:,:,iz, :,:) = xl%c%displacements_P(:,:,iz, :,:) * scale
     end do
!$omp end parallel do

  end do ! i_step

  call divide_and_write_average(xl,avg)
end subroutine velocity_scaling
