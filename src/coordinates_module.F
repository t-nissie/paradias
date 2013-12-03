! coordinates_module.F -*-f90-*-
! Time-stamp: <2013-12-03 17:29:22 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module coordinates_module
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'

  type coordinates_type
     integer(C_INTPTR_T)                :: local_Lz, local_z_offset
     !integer(C_INTPTR_T)               :: local_Ly, local_y_offset

     real(C_DOUBLE),            pointer ::       displacements_P(:,:,:, :, :)   ! momenta

     type(C_PTR)                        :: displacements_R_ptr, displacements_k_ptr
     real(C_DOUBLE),            pointer ::       displacements_R(:,:,:, :, :), displacements_R_tmp(:)
     complex(C_DOUBLE_COMPLEX), pointer ::       displacements_k(:,:,:, :, :), displacements_k_tmp(:)
     type(C_PTR) ::                         plan_displacements_R2k

     type(C_PTR)                        :: d_V_d_displacements_k_ptr
     complex(C_DOUBLE_COMPLEX), pointer :: d_V_d_displacements_k(:,:,:, :, :), d_V_d_displacements_k_tmp(:)
     real(C_DOUBLE),            pointer :: d_V_d_displacements_R(:,:,:, :, :)  ! in_place
     type(C_PTR) ::                   plan_d_V_d_displacements_k2R
  end type coordinates_type

contains
  subroutine coordinates_allocate(n_atoms, Lx, Ly, Lz, c)
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
    integer(C_INTPTR_T), intent(in)       :: n_atoms
    integer(C_INTPTR_T), intent(in)       :: Lx, Ly, Lz
    type(coordinates_type), intent(inout) :: c
    integer(C_INTPTR_T)                   :: alloc_local, i, L(3) !i_atom, alpha

#   if defined(HAVE_LIBFFTW3)
    integer alpha, i_atom
#   include "fftw3.f"
#   endif

    ! Get data size and allocate (note dimension reversal)
    L(1) = Lz
    L(2) = Ly
    L(3) = Lx
    alloc_local = fftw_mpi_local_size_many(3, L, 3*n_atoms, &
          & FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &
          & c%local_Lz, c%local_z_offset)

    write(6,*) alloc_local

    c%displacements_R_ptr = fftw_alloc_real(2 * alloc_local)
    call c_f_pointer(c%displacements_R_ptr, c%displacements_R_tmp, [2*(Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%displacements_R(1:n_atoms, 1:3, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1) => c%displacements_R_tmp

    c%displacements_k_ptr = fftw_alloc_complex( alloc_local)
    call c_f_pointer(c%displacements_k_ptr, c%displacements_k_tmp, [  (Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%displacements_k(1:n_atoms, 1:3, 0:Lx/2,         0:Ly-1, 0:c%local_Lz-1) => c%displacements_k_tmp

    ! First touch
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%displacements_R(:, :, :, :, i) = 0.0d0
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%displacements_k(:, :, :, :, i) = (0.0d0, 0.0d0)
    end do
!$omp end parallel do

     c%plan_displacements_R2k = fftw_mpi_plan_many_dft_r2c(3, L, 3*n_atoms, &
          & FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
          & c%displacements_R, c%displacements_k, &
          & MPI_COMM_WORLD, FFTW_MEASURE)
  end subroutine coordinates_allocate
end module coordinates_module
