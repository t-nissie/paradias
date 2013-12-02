! average_module.F -*-f90-*-
! Time-stamp: <2013-03-17 18:31:17 t-nissie>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module average_module
  use energy_module
  implicit none

  type average_type
     type(energy_type) energy
     real*8 u(3)
     real*8 uu(6)
  end type average_type

contains
  subroutine reset_average(avg)
    implicit none
    !--- arguments -------------------------------------
    type(average_type), intent(inout) :: avg
    !--- END of variable definitions -------------------
    avg%energy            = energy_type(0.0d0, 0.0d0, 0.0d0)
    avg%u(:)              = 0.0d0
    avg%uu(:)             = 0.0d0
  end subroutine reset_average

  subroutine add_to_average(u_tmp,e,avg,L)
    use xtal_module
    use coordinates_module
    implicit none
    !--- arguments -------------------------------------
    real*8,                  intent(in)    :: u_tmp(3)   !!! SINGLE ATOM !!!
    type(energy_type),       intent(in)    :: e
    type(average_type),      intent(inout) :: avg
    integer,                 intent(in)    :: L(3)
    !--- END of arguments ------------------------------

    avg%energy = avg%energy + e

    avg%u(:)  = avg%u(:) + u_tmp(:)
    avg%uu(1) = avg%uu(1) + u_tmp(1)**2
    avg%uu(2) = avg%uu(2) + u_tmp(2)**2
    avg%uu(3) = avg%uu(3) + u_tmp(3)**2
    avg%uu(4) = avg%uu(4) + u_tmp(2)*u_tmp(3)
    avg%uu(5) = avg%uu(5) + u_tmp(3)*u_tmp(1)
    avg%uu(6) = avg%uu(6) + u_tmp(1)*u_tmp(2)

    ! if (detail) then
    !    u_tmp(:) = u_tmp(:) / c%n_sum_p
    !    write(message,'("<u>   =",3f10.6)') u_tmp(1:3)
    !    call msg_ez(UNIT_LOG,message)
    !    write(message,'("sigma =",3f10.6)') sqrt(uu_tmp(1:3)/c%n_sum_p - u_tmp(1:3)**2)
    !    call msg_ez(UNIT_LOG,message)
    ! end if
  end subroutine add_to_average

  subroutine divide_and_write_average(xl,avg)
    use xtal_module
    implicit none
    type(   xtal_type), intent(in)    :: xl
    type(average_type), intent(inout) :: avg

    avg%energy = div_energy(avg%energy, xl%n_average)

    avg%u(:)  = avg%u(:)  / xl%n_average * xl%N_in_super_cell_inv
    avg%uu(:) = avg%uu(:) / xl%n_average * xl%N_in_super_cell_inv

    open(unit=UNIT_AVERAGE,file=trim(xl%filename)//'.avg', status='REPLACE')
    write(UNIT_AVERAGE,'(f5.3,3f9.5,1(e15.5,2e13.5),1(e15.5,5e13.5),e15.5,2e13.5)') &
         & xl%temperature, xl%external_H_field(:), avg%u(1:3), avg%uu(1:6), avg%energy%kinetic, avg%energy%quadratic, avg%energy%externalH
    close(UNIT_AVERAGE)
  end subroutine divide_and_write_average
end module average_module
