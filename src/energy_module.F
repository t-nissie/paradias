! energy_module.F -*-f90-*-
! Time-stamp: <2013-03-17 18:22:57 t-nissie>
! Author: Takeshi NISHIMATSU
!!
module energy_module
  implicit none
  type energy_type
     real*8 kinetic
     real*8 quadratic
     real*8 externalH
  end type energy_type

interface operator(+)
   module procedure add_energy
end interface

contains
  function add_energy(a,b)
    implicit none
    type(energy_type)             :: add_energy
    type(energy_type), intent(in) :: a, b
    add_energy%kinetic   = a%kinetic   + b%kinetic
    add_energy%quadratic = a%quadratic + b%quadratic
    add_energy%externalH = a%externalH + b%externalH
  end function add_energy

  function div_energy(a,n)
    implicit none
    type(energy_type)             :: div_energy
    type(energy_type), intent(in) :: a
    integer, intent(in)           :: n
    div_energy%kinetic    = a%kinetic    / n
    div_energy%quadratic  = a%quadratic  / n
    div_energy%externalH  = a%externalH  / n
  end function div_energy

  ! function sum_potential_energies(e)
  !   implicit none
  !   real*8 :: sum_potential_energies
  !   type(energy_type), intent(in) :: e
  !   sum_potential_energies &
  !        = e%dipole_dipole &
  !        + e%dipole_E_field&
  !        + e%unharmonic    &
  !        + e%homo_strain   &
  !        + e%homo_coupling &
  !        + e%inho_strain   &
  !        + e%inho_coupling
  ! end function sum_potential_energies

  ! subroutine write_energy(p,c,e,fd)
  !   use Param_module
  !   use Coord_module
  !   implicit none
  !   !--- arguments -------------------------------------
  !   type(Param_type),  intent(in)    :: p
  !   type(Coord_type),  intent(inout) :: c
  !   type(energy_type), intent(inout) :: e
  !   integer fd
  !   !--- local -----------------------------------------
  !   character(30), parameter :: fmt = "(a15,f18.12,' [eV/u.c.]')"
  !   character(45) message
  !   !--- END of variable definitions -------------------
  !   write(message,fmt) 'kinetic   ', e%kinetic*p%Ninv;    call msg_ez(fd,message)
  !   write(message,fmt) 'quadratic   ', e%quadratic/c%s_Nose**2*p%Ninv
  !                                                                   call msg_ez(fd,message)
  !   call msg_ez(fd,'--------------------------------------------')
  !   write(message,fmt) 'total_energy   ', e%total*p%Ninv;           call msg_ez(fd,message)
  ! end subroutine write_energy
end module energy_module
