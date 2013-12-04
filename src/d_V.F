! d_V.F -*-f90-*-
! Time-stamp: <2013-12-04 21:54:40 takeshi>
! Author: Takeshi NISHIMATSU
!!
subroutine d_V(xl,e_quadratic)
  use xtal_module
  use coordinates_module
  implicit none
  type(xtal_type), intent(inout) :: xl
  real*8,          intent(out)   :: e_quadratic
  integer :: ix,iy,iz, I_atom

  !original call displacements_FFT_R2k(xl%n_atoms,xl%L,xl%c)
  call fftw_mpi_execute_dft_r2c(c%plan_displacements_R2k, &
       & c%displacements_R, c%displacements_k)

  call d_V_d_displacements(xl%n_atoms,xl%L,xl%c,xl%qd,e_quadratic)

!$omp parallel do private(iz,iy,ix,I_atom)
  do iz=0, xl%L(3)-1
     do iy=0, xl%L(2)-1
        do ix=0, xl%L(1)-1
           do I_atom = 1, xl%n_atoms
              xl%c%d_V_d_displacements_R(ix,iy,iz,:,I_atom) &
          & = xl%c%d_V_d_displacements_R(ix,iy,iz,:,I_atom) - xl%external_H_field(:)
           end do
        end do
     end do
  end do
!$omp end parallel do
end subroutine d_V
