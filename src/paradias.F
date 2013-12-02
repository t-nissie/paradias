! paradias.F -*-f90-*-
! Time-stamp: <2013-03-31 17:35:24 takeshi>
! Author: Takeshi NISHIMATSU
!!
#include "define.h"
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

program paradias
  use Xtal_module
  use quadratic_module
  implicit none
  type(Xtal_type) :: xl
  integer argn,i
  character (len=100) :: package_string
  write(package_string,'(a)') PACKAGE_STRING

  argn = command_argument_count()
  if (argn.eq.0) stop 'An input file required.'

  do i=1,argn
     call get_command_argument(i,xl%filename)
     call xtal_construct(xl)
     open(UNIT_LOG,FILE=trim(xl%filename)//'.log', STATUS='REPLACE')
     open(UNIT_INPUT, FILE=xl%filename, STATUS='OLD')
       call xtal_read(xl,UNIT_INPUT)
     close(UNIT_INPUT)
     call xtal_calc_other_parameters(xl)

     call dynamical_matrix_allocate(xl%n_atoms, xl%L, xl%qd)

     call dynamical_matrix(xl)

     if (xl%axes(1)%label.ne.'END') then
        call dispersion(xl)
     end if

     call coordinates_allocate(xl%n_atoms, xl%L, xl%c)

     call displacements_initialize(xl%n_atoms, xl%L, xl%init_avg, xl%init_dev, xl%c, trim(xl%filename) // '.restart')

     if (xl%method.eq.'vs') then
        call velocity_scaling(xl)
     else if (xl%method.eq.'lf') then
        call leapfrog(xl)
     else
        stop 'paradias: Illegal method.'
     end if

     call coordinates_write(xl%n_atoms, xl%L, xl%c, trim(xl%filename) // '.coord')

     call xtal_destruct(xl)
     close(UNIT_LOG)
  end do
end program paradias
