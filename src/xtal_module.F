! xtal_module.F -*-f90-*-
! Time-stamp: <2013-03-31 17:27:32 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module Xtal_module
  use Atom_module
  use quadratic_module
  use coordinates_module
  implicit none
  type k_point_type
     character (len=30) label
     real*8 k(3)
  end type k_point_type

  type Axis_type
     character (len=30) label
     integer n_divide
  end type Axis_type

  type Xtal_type
     character (len=FILENAME_LEN) filename
     character (len=300) title
     character (len=30)  method
     integer n_atoms
     real*8 prim(3,3)   !primitive vector (column vector)
     real*8 r_prim(3,3) !primitive reciprocal lattice vector (row vector)
     real*8 volume
     integer L(3)
     integer n_thermalize
     integer n_average
     integer N_in_super_cell
     real*8  N_in_super_cell_inv
     real*8 dt          ! Time step width delta t [pico second]
     real*8 dt_2        ! dt/2
     real*8 temperature ! [J]
     real*8 external_H_field(3)   ! [H/J]
     real*8 init_avg(3)
     real*8 init_dev(3)
     type(k_point_type)      :: k_points(0:MAX_AXIS)
     type(Axis_type)         :: axes(MAX_AXIS)
     type(Atom_type)         :: atom_list
     type(quadratic_type)    :: qd
     type(coordinates_type)  :: c
     real*8  epsilon_inf(3,3)
  end type Xtal_type

contains
  subroutine xtal_destruct(xl)
    implicit none
    type(Xtal_type), intent(inout) :: xl
    call atom_list_destruct(xl%atom_list)
  end subroutine xtal_destruct

  subroutine xtal_construct(xl)
    implicit none
    type(Xtal_type), intent(inout) :: xl

    !default values
    xl%title = ''
    xl%method = 'vs'
    xl%prim(:,:)        = reshape((/1.0d0, 0.0d0, 0.0d0, &
         &                          0.0d0, 1.0d0, 0.0d0, &
         &                          0.0d0, 0.0d0, 1.0d0/),(/3,3/))

    xl%epsilon_inf(:,:) = reshape((/1.0d0, 0.0d0, 0.0d0, &
         &                          0.0d0, 1.0d0, 0.0d0, &
         &                          0.0d0, 0.0d0, 1.0d0/),(/3,3/))
    xl%L(:) = (/32,32,32/)
    xl%n_thermalize = 400000
    xl%n_average    = 100000
    xl%dt   = 0.002d0         ! Time step width delta t [pico second]
    xl%dt_2 = 0.001d0        ! dt/2
    xl%temperature = 1.0d0
    xl%external_H_field(:) = 0.0d0
    xl%init_avg(:) = (/0.0d0,0.0d0,0.0d0/)
    xl%init_dev(:) = (/0.3d0,0.3d0,0.3d0/)
    xl%qd%elevation = 0.0d0
    xl%qd%n_max(:)       = (/0,0,0/)
    xl%qd%G_max(:)       = (/0,0,0/)
  end subroutine xtal_construct

  subroutine xtal_read(xl,fd)
    implicit none
    type(Xtal_type), intent(inout) :: xl
    integer, intent(in) :: fd
    character (len=300) str
    character (len=30)  name, dmy1, dmy2
    integer io, line_number, i

    ! Read parameters from the imput file
    line_number = 0
    i = 0
    do
       read(fd, '(a)', IOSTAT=io) str
       if (io.ne.0) exit
       line_number = line_number + 1
       read(str,*,IOSTAT=io) name   ! extract name of a parameter
       if (io.ne.0 .or. name(1:1).eq.'#') then  ! str is a blank line or a comment line
          write(UNIT_LOG,'(a)') TRIM(str)
          cycle
       end if
       select case(TRIM(name))
          case('title')
             read(str,*) dmy1, dmy2, xl%title
             write(UNIT_LOG,'(a,a,a)') "title = '", TRIM(xl%title), "'"
          case('method')
             read(str,*) dmy1, dmy2, xl%method
             write(UNIT_LOG,'(a,a,a)') "method = '", TRIM(xl%method), "'"
          case('n_thermalize')
             read(str,*) dmy1, dmy2, xl%n_thermalize
             write(UNIT_LOG,'(a,i10)')  'n_thermalize = ', xl%n_thermalize
          case('n_average')
             read(str,*) dmy1, dmy2, xl%n_average
             write(UNIT_LOG,'(a,i10)')  'n_average = ', xl%n_average
          case('L')
             read(str,*) dmy1, dmy2, xl%L(:)
             write(UNIT_LOG,'(a,3i4)')  'L = ', xl%L(:)
          case('init_avg')
             read(str,*) dmy1, dmy2, xl%init_avg
             write(UNIT_LOG,'(a,3f10.5)')  'init_avg = ', xl%init_avg
          case('init_dev')
             read(str,*) dmy1, dmy2, xl%init_dev
             write(UNIT_LOG,'(a,3f10.5)')  'init_dev = ', xl%init_dev
          case('dt')
             read(str,*) dmy1, dmy2, xl%dt
             write(UNIT_LOG,'(a,f10.5)')  'dt = ', xl%dt
             xl%dt_2 = xl%dt / 2
          case('temperature')
             read(str,*) dmy1, dmy2, xl%temperature
             write(UNIT_LOG,'(a,f10.5)')  'temperature = ', xl%temperature
          case('external_H_field')
             read(str,*) dmy1, dmy2, xl%external_H_field(:)
             write(UNIT_LOG,'(a,3f10.5)')  'external_H_field = ', xl%external_H_field
          case('n_max')
             read(str,*) dmy1, dmy2, xl%qd%n_max(:)
             write(UNIT_LOG,'(a,3i4)')  'n_max = ', xl%qd%n_max(:)
          case('G_max')
             read(str,*) dmy1, dmy2, xl%qd%G_max(:)
             write(UNIT_LOG,'(a,3i4)')  'G_max = ', xl%qd%G_max(:)
          case('elevation')
             read(str,*) dmy1, dmy2, xl%qd%elevation
          case('prim1')
             read(str,*) dmy1, dmy2, xl%prim(:,1)
             write(UNIT_LOG,'(a,3f22.16)')  'prim1 =', xl%prim(:,1)
          case('prim2')
             read(str,*) dmy1, dmy2, xl%prim(:,2)
             write(UNIT_LOG,'(a,3f22.16)')  'prim2 =', xl%prim(:,2)
          case('prim3')
             read(str,*) dmy1, dmy2, xl%prim(:,3)
             write(UNIT_LOG,'(a,3f22.16)')  'prim3 =', xl%prim(:,3)
          case('epsilon')
             read(str,*) dmy1, dmy2, &
                  & xl%epsilon_inf(1,1), &
                  & xl%epsilon_inf(2,2), &
                  & xl%epsilon_inf(3,3), &
                  & xl%epsilon_inf(2,3), &
                  & xl%epsilon_inf(3,1), &
                  & xl%epsilon_inf(1,2)
             xl%epsilon_inf(3,2) = xl%epsilon_inf(2,3)
             xl%epsilon_inf(1,3) = xl%epsilon_inf(3,1)
             xl%epsilon_inf(2,1) = xl%epsilon_inf(1,2)
             write(UNIT_LOG,'(a)')            'epsilon_inf ='
             write(UNIT_LOG,'(3f22.16)') xl%epsilon_inf(:,:)
          case('atom')
             call atom_list_push(xl%atom_list, name, str)
          case('k-point')
             read(str,*) dmy1, dmy2, xl%k_points(i)%label, xl%k_points(i)%k(:)
             write(UNIT_LOG,'(a,a,a,3f22.16)')  "k-point = '", TRIM(xl%k_points(i)%label), "'",&
                  &                                          xl%k_points(i)%k(:)
             i = i + 1
          case('k-point-Cartesian')
             read(str,*) dmy1, dmy2, xl%k_points(i)%label, xl%k_points(i)%k(:)
             xl%k_points(i)%k(:) = matmul(xl%k_points(i)%k,xl%prim)
             write(UNIT_LOG,'(a,a,a,3f22.16)')  "k-point = '", TRIM(xl%k_points(i)%label), "'",&
                  &                                          xl%k_points(i)%k(:)
             i = i + 1
          case('axis')
             if (i.eq.0) stop ': Input error: axis befor k-point'
             read(str,*) dmy1, dmy2, xl%axes(i)%label, xl%axes(i)%n_divide
             write(UNIT_LOG,'(a,a,a,i3)')  "axis = '", TRIM(xl%axes(i)%label), "'",&
                  &                                  xl%axes(i)%n_divide
          case default
             write(0,'(a,a,i3,a,a)') trim(xl%filename),':', line_number, ':', trim(str)
             write(0,'(a,a,i3,a)') __FILE__, ':', __LINE__, ': Input error: No such parameter'
             stop 1
       end select
    end do

    if (i.eq.0) i=1   ! There was no k-point and axis input.
    xl%axes(i)%label = 'END'

    xl%n_atoms = atom_list_count(xl%atom_list)
    write(UNIT_LOG,'(a,i3)')  'n_atoms = ', xl%n_atoms
    if (xl%n_atoms.eq.0) then
       write(UNIT_LOG,'(a,a,i3,a,a)') trim(xl%filename),':', line_number, ':', trim(str)
       stop ': Input error: No atom or dipole'
    end if

    dmy1 = dmy2; dmy2 = dmy1   ! suppress compiler warning
  end subroutine xtal_read

  function cross_product(a,b)
    implicit none
    real*8 cross_product(3), a(3), b(3)
    cross_product(1)=a(2)*b(3)-a(3)*b(2)
    cross_product(2)=a(3)*b(1)-a(1)*b(3)
    cross_product(3)=a(1)*b(2)-a(2)*b(1)
  end function cross_product
end module Xtal_module
