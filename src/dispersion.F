! dispersion.F -*-f90-*-
! Time-stamp: <2013-03-16 12:29:47 t-nissie>
! Author: Takeshi NISHIMATSU
!!
#include "define.h"
subroutine dispersion(xl)
  use Xtal_module
  implicit none
  type(Xtal_type), intent(in) :: xl
  character (len=FILENAME_LEN) :: dispersion_filename
  real*8      :: last_eigenvalues(3*xl%n_atoms)
  real*8      ::  new_eigenvalues(3*xl%n_atoms)
  complex*16  :: last_eigenvectors(3*xl%n_atoms,3*xl%n_atoms)
  complex*16  ::  new_eigenvectors(3*xl%n_atoms,3*xl%n_atoms)
  real*8 k2k(3), k2k_Cartesian(3)
  real*8   k(3),   k_Cartesian(3)
  integer :: i, j, l
  real*8 :: stacked_length, axis_length, plot_length
  character (len=50) :: fmtstr
  write(fmtstr,'(a,i2,a)') "('#',", 3*xl%n_atoms, '(f12.5,f10.5))'

  write(dispersion_filename,'(2a)') TRIM(xl%filename), '.dispersion'
  open(UNIT_DISPERSION, FILE=dispersion_filename, STATUS='REPLACE')

  i = 1
  stacked_length = 0.0d0
  call dispersion_gnuplot_open(xl)
  do
     if (i.eq.1) then
        write(UNIT_DISPERSION,'(a,a)') '# k-point: ', xl%k_points(i-1)%label
     else
        write(UNIT_DISPERSION,'(a,a)') '# axis: ', xl%axes(i-1)%label
        write(UNIT_DISPERSION,'(a,a)') '# k-point: ', xl%k_points(i-1)%label
        write(UNIT_DISPERSION,'(a)'  ) ''
        write(UNIT_DISPERSION,'(a,a)') '# k-point: ', xl%k_points(i-1)%label
     end if
     if (xl%axes(i)%label.eq.'END') exit
     write(UNIT_DISPERSION,'(a,a,/,a)') '# axis: ', xl%axes(i)%label, '#'
     k2k(:) = xl%k_points(i)%k(:) - xl%k_points(i-1)%k(:)
     k2k_Cartesian(:) = matmul(k2k,xl%r_prim)
     axis_length = sqrt(dot_product(k2k_Cartesian, k2k_Cartesian))
     call dispersion_gnuplot_axis_label(xl, i, stacked_length+axis_length/2)
     do j = 0, xl%axes(i)%n_divide
        plot_length = dble(j)/xl%axes(i)%n_divide*axis_length + stacked_length
        k(:)        = dble(j)/xl%axes(i)%n_divide*k2k(:) + xl%k_points(i-1)%k(:)
        k_Cartesian(:) = matmul(k,xl%r_prim)
        call k2eigenvalues(3*xl%n_atoms, k, xl%L, xl%qd%C_tilde, &
             & new_eigenvalues, new_eigenvectors)
        call sort_eigen(xl, plot_length, j, new_eigenvalues, &
             & last_eigenvalues, last_eigenvectors, new_eigenvectors)
        write(UNIT_DISPERSION,'(100f11.6)') plot_length, k_Cartesian(:), last_eigenvalues(:)
        do l = 1, 3*xl%n_atoms
           write(UNIT_DISPERSION,fmt=fmtstr) last_eigenvectors(l,:)
        end do   ! write eigenvalues as comments in the dispersion file
        write(UNIT_DISPERSION,'(a)') '#'
     end do
     stacked_length = stacked_length + axis_length
     call dispersion_gnuplot_xtics(xl,i,stacked_length)
     i = i + 1
  end do
  call dispersion_gnuplot_close(xl,last_eigenvectors)

  close(UNIT_DISPERSION)
end subroutine dispersion
