! atom_module.F -*-f90-*-
! Time-stamp: <2013-03-16 12:40:14 t-nissie>
! Author: Takeshi NISHIMATSU
!!
#include "define.h"
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

module Atom_module
  implicit none
  type Atom_type
     real*8 :: r(3)       = (/0.0d0, 0.0d0, 0.0d0/)
     real*8 :: Z_tensor(3,3) = reshape((/1.0d0, 0.0d0, 0.0d0, &
                                         0.0d0, 1.0d0, 0.0d0, &
                                         0.0d0, 0.0d0, 1.0d0/),(/3,3/))
     real*8 :: mass     = 1.0d0
     real*8 :: mass_inv = 1.0d0
     type(Atom_type), pointer :: next => null()
  end type Atom_type

contains
  subroutine atom_list_destruct(list)
    implicit none
    type(Atom_type), intent(inout) :: list
    type(Atom_type), pointer :: current_atom, next_atom
    next_atom => list%next
    do while (associated(next_atom))
       current_atom => next_atom
       next_atom => next_atom%next
       deallocate(current_atom)
    end do
    nullify(list%next)
  end subroutine atom_list_destruct

  subroutine atom_list_push(list, name, str)
    implicit none
    type(Atom_type), target, intent(in) :: list
    character(len=*), intent(in) :: name, str
    type(Atom_type), pointer :: p_atom
    character (len=30) dmy1, dmy2
    integer io

    p_atom => list
    do while (associated(p_atom%next))
       p_atom => p_atom%next
    end do

    allocate(p_atom%next)   !PGI Fortran Ver.6.0-5 cannot compile it!

    if (name.eq.'atom') then
       read(str, *, IOSTAT=io) dmy1, dmy2, &
            & p_atom%next%r(:),   &
            & p_atom%next%mass,   &
            & p_atom%next%Z_tensor(1,1), &
            & p_atom%next%Z_tensor(2,2), &
            & p_atom%next%Z_tensor(3,3), &
            & p_atom%next%Z_tensor(2,3), &
            & p_atom%next%Z_tensor(3,1), &
            & p_atom%next%Z_tensor(1,2)
       p_atom%next%Z_tensor(3,2) = p_atom%next%Z_tensor(2,3)
       p_atom%next%Z_tensor(1,3) = p_atom%next%Z_tensor(3,1)
       p_atom%next%Z_tensor(2,1) = p_atom%next%Z_tensor(1,2)
       if (io.ne.0) then
          p_atom%next%Z_tensor(:,:) = 0.0d0
          read(str, *) dmy1, dmy2, &
               & p_atom%next%r(:),   &
               & p_atom%next%mass,   &
               & p_atom%next%Z_tensor(1,1)
          p_atom%next%Z_tensor(2,2) = p_atom%next%Z_tensor(1,1)
          p_atom%next%Z_tensor(3,3) = p_atom%next%Z_tensor(1,1)
       end if
       p_atom%next%mass_inv = 1.0d0 / p_atom%next%mass
       write(UNIT_LOG,'(a,3f22.16,7f6.2)')  'atom =', &
            &                  p_atom%next%r(:),   &
            &                  p_atom%next%mass,   &
            &                  p_atom%next%Z_tensor(1,1), &
            &                  p_atom%next%Z_tensor(2,2), &
            &                  p_atom%next%Z_tensor(3,3), &
            &                  p_atom%next%Z_tensor(2,3), &
            &                  p_atom%next%Z_tensor(3,1), &
            &                  p_atom%next%Z_tensor(1,2)
    else
       stop '"atom" is requred.'
    end if
    nullify(p_atom%next%next)  ! make sure
    dmy1 = dmy2; dmy2 = dmy1   ! suppress compiler warning
  end subroutine atom_list_push

  integer function atom_list_count(list)
    implicit none
    type(Atom_type), intent(in) :: list
    type(Atom_type), pointer :: p_atom
    atom_list_count = 0
    p_atom => list%next
    do while (associated(p_atom))
       atom_list_count = atom_list_count + 1
       p_atom => p_atom%next
    end do
  end function atom_list_count

  logical function atom_list_each(p_atom)
    implicit none
    type(Atom_type), pointer :: p_atom
    p_atom => p_atom%next
    atom_list_each = associated(p_atom)
  end function atom_list_each

  logical function atom_list_each_with_index(i_atom,i)
    implicit none
    type(Atom_type), pointer :: i_atom
    integer, intent(inout) :: i
    i_atom => i_atom%next
    i = i + 1
    atom_list_each_with_index = associated(i_atom)
  end function atom_list_each_with_index
end module Atom_module
