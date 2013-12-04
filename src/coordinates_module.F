! coordinates_module.F -*-f90-*-
! Time-stamp: <2013-12-05 07:00:06 takeshi>
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
     integer(C_INTPTR_T)                :: vector_local_size, local_Lz, local_z_offset

     real(C_DOUBLE),            pointer ::       displacements_P(:,: ,:,:,:)   ! momenta

     type(C_PTR)                        :: displacements_R_ptr, displacements_k_ptr
     real(C_DOUBLE),            pointer ::       displacements_R(:,: ,:,:,:), displacements_R_tmp(:)
     complex(C_DOUBLE_COMPLEX), pointer ::       displacements_k(:,: ,:,:,:), displacements_k_tmp(:)
     type(C_PTR) ::                         plan_displacements_R2k

     type(C_PTR)                        :: d_V_d_displacements_k_ptr   ! in-place, no pointer for R
     complex(C_DOUBLE_COMPLEX), pointer :: d_V_d_displacements_k(:,: ,:,:,:), d_V_d_displacements_k_tmp(:)
     real(C_DOUBLE),            pointer :: d_V_d_displacements_R(:,: ,:,:,:), d_V_d_displacements_R_tmp(:)
     type(C_PTR) ::                   plan_d_V_d_displacements_k2R

     integer n_max(3)
     integer G_max(3)
     real*8  elevation   ! 2*elevation will be added to diagonal terms
     type(C_PTR)                        ::      C_bar_DD_k_ptr   ! in-place, no pointer for R
     complex(C_DOUBLE_COMPLEX), pointer ::      C_bar_DD_k(:,:, :,:, :,:,:), C_bar_DD_k_tmp(:)
     complex(C_DOUBLE_COMPLEX), pointer ::      C_bar_DD_R(:,:, :,:, :,:,:), C_bar_DD_R_tmp(:)

     type(C_PTR)                        :: C_tilde_short_k_ptr   ! in-place, no pointer for R
     complex(C_DOUBLE_COMPLEX), pointer :: C_tilde_short_k(:,:, :,:, :,:,:), C_tilde_short_k_tmp(:)
     complex(C_DOUBLE_COMPLEX), pointer :: C_tilde_short_R(:,:, :,:, :,:,:), C_tilde_short_R_tmp(:)

     complex(C_DOUBLE_COMPLEX), pointer ::       C_tilde_k(:,: ,:,: ,:,:,:)
  end type coordinates_type

contains
  subroutine coordinates_allocate(n_atoms, Lx, Ly, Lz, c)
    use, intrinsic :: iso_c_binding
    implicit none
    include 'mpif.h'
    integer(C_INTPTR_T), intent(in)       :: n_atoms
    integer(C_INTPTR_T), intent(in)       :: Lx, Ly, Lz
    type(coordinates_type), intent(inout) :: c
    integer(C_INTPTR_T)                   :: matrix_local_size, local_Lz_matrix, local_z_offset_matrix, i, L(3) !i_atom, alpha
    integer                               :: i_error
#   if defined(HAVE_LIBFFTW3)
    integer alpha, i_atom
#   include "fftw3.f"
#   endif

! Get data size
    L(1) = Lz
    L(2) = Ly   ! note dimension reversal
    L(3) = Lx
    c%vector_local_size = fftw_mpi_local_size_many(3, L,  3*n_atoms, &
          & FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &
          & c%local_Lz, c%local_z_offset)

! Make sure
    matrix_local_size   = fftw_mpi_local_size_many(3, L, (3*n_atoms)**2, &
          & FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &
          & local_Lz_matrix, local_z_offset_matrix)
    if (c%vector_local_size*3*n_atoms.ne.matrix_local_size) then
       write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
            & ': Error in getting data size. c%vector_local_size*3*n_atoms != matrix_local_size.'
       call MPI_ABORT(MPI_COMM_WORLD, ERROR_COORDINATES_MODULE_SIZE_NOT_SAME, i_error)
    end if
    if (c%local_Lz.ne.local_Lz_matrix) then
       write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
            & ': Error in getting data size. c%local_Lz != local_Lz_matrix.'
       call MPI_ABORT(MPI_COMM_WORLD, ERROR_COORDINATES_MODULE_LOCAL_LZ_NOT_SAME, i_error)
    end if
    if (c%local_z_offset.ne.local_z_offset_matrix) then
       write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
            & ': Error in getting data size. c%local_z_offset != local_z_offset_matrix.'
       call MPI_ABORT(MPI_COMM_WORLD, ERROR_COORDINATES_MODULE_LOCAL_Z_OFFSET_NOT_SAME, i_error)
    end if

! Allocate
    c%displacements_R_ptr = fftw_alloc_real(2 * c%vector_local_size)
    call c_f_pointer(c%displacements_R_ptr, &
         &           c%displacements_R_tmp, &
         &           [2*(Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%displacements_R(1:3,       1:n_atoms, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1) => c%displacements_R_tmp

    c%displacements_k_ptr = fftw_alloc_complex( c%vector_local_size)
    call c_f_pointer(c%displacements_k_ptr, &
         &           c%displacements_k_tmp, &
         &           [  (Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%displacements_k(1:3,       1:n_atoms, 0:Lx/2,         0:Ly-1, 0:c%local_Lz-1) => c%displacements_k_tmp

    c%d_V_d_displacements_k_ptr = fftw_alloc_complex(c%vector_local_size)
    call c_f_pointer(c%d_V_d_displacements_k_ptr, &
         &           c%d_V_d_displacements_k_tmp, &
         &           [  (Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%d_V_d_displacements_k(1:3, 1:n_atoms, 0:Lx/2,         0:Ly-1, 0:c%local_Lz-1) => c%d_V_d_displacements_k_tmp

    ! in-place, no pointer for R
    call c_f_pointer(c%d_V_d_displacements_k_ptr, &
         &           c%d_V_d_displacements_R_tmp, &
         &           [2*(Lx/2+1)*Ly*c%local_Lz*3*n_atoms])
    c%d_V_d_displacements_R(1:3, 1:n_atoms, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1) => c%d_V_d_displacements_R_tmp

    allocate(c%displacements_P(1:3, 1:n_atoms, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1))

    c%C_bar_DD_k_ptr  = fftw_alloc_complex(c%vector_local_size*3*n_atoms)
    call c_f_pointer(c%C_bar_DD_k_ptr, &
         &           c%C_bar_DD_k_tmp, &
         &           [  (Lx/2+1)*Ly*c%local_Lz*(3*n_atoms)**2])
    c%C_bar_DD_k(1:3, 1:n_atoms, 1:3, 1:n_atoms, 0:Lx/2, 0:Ly-1, 0:c%local_Lz-1) => c%C_bar_DD_k_tmp

    ! in-place, no pointer for R
    call c_f_pointer(c%C_bar_DD_k_ptr, &
         &           c%C_bar_DD_R_tmp, &
         &           [2*(Lx/2+1)*Ly*c%local_Lz*(3*n_atoms)**2])
    c%C_bar_DD_R(1:3, 1:n_atoms, 1:3, 1:n_atoms, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1) => c%C_bar_DD_R_tmp

    c%C_tilde_short_k_ptr  = fftw_alloc_complex(c%vector_local_size*3*n_atoms)
    call c_f_pointer(c%C_tilde_short_k_ptr, &
         &           c%C_tilde_short_k_tmp, &
         &           [  (Lx/2+1)*Ly*c%local_Lz*(3*n_atoms)**2])
    c%C_tilde_short_k(1:3, 1:n_atoms, 1:3, 1:n_atoms, 0:Lx/2, 0:Ly-1, 0:c%local_Lz-1) => c%C_tilde_short_k_tmp

    ! in-place, no pointer for R
    call c_f_pointer(c%C_tilde_short_k_ptr, &
         &           c%C_tilde_short_R_tmp, &
         &           [2*(Lx/2+1)*Ly*c%local_Lz*(3*n_atoms)**2])
    c%C_tilde_short_R(1:3, 1:n_atoms, 1:3, 1:n_atoms, 0:2*(Lx/2+1)-1, 0:Ly-1, 0:c%local_Lz-1) => c%C_tilde_short_R_tmp

    allocate(c%C_tilde_k(1:3, 1:n_atoms, 1:3, 1:n_atoms, 0:Lx/2, 0:Ly-1, 0:c%local_Lz-1))

! First touch
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%displacements_P(:, :, :, :, i) = 0.0d0
    end do
!$omp end parallel do
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
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%d_V_d_displacements_k(:, :, :, :, i) = (0.0d0, 0.0d0)
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%C_bar_DD_k(:,:, :,:, :,:,i) = (0.0d0, 0.0d0)
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%C_tilde_short_k(:,:, :,:, :,:,i) = (0.0d0, 0.0d0)
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0, c%local_Lz-1
       c%C_tilde_k(:,:, :,:, :,:,i) = (0.0d0, 0.0d0)
    end do
!$omp end parallel do

! Planning
     c%plan_displacements_R2k = fftw_mpi_plan_many_dft_r2c(3, L, 3*n_atoms, &
          & FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
          & c%displacements_R, &
          & c%displacements_k, &
          & MPI_COMM_WORLD, FFTW_MEASURE)

     c%plan_d_V_d_displacements_k2R = fftw_mpi_plan_many_dft_c2r(3, L, 3*n_atoms, &
          & FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
          & c%d_V_d_displacements_k, &
          & c%d_V_d_displacements_R, &
          & MPI_COMM_WORLD, FFTW_MEASURE)
  end subroutine coordinates_allocate

  subroutine coordinates_free(c)
    implicit none
    type(coordinates_type), intent(inout) :: c
    call fftw_destroy_plan(c%plan_displacements_R2k)
    call fftw_destroy_plan(c%plan_d_V_d_displacements_k2R)
    call fftw_free(c%displacements_R_ptr)
    call fftw_free(c%displacements_k_ptr)
    call fftw_free(c%d_V_d_displacements_k_ptr)
    call fftw_free(c%C_bar_DD_k_ptr)
    call fftw_free(c%C_tilde_short_k_ptr)
    deallocate(c%displacements_P)
    deallocate(c%C_tilde_k)
  end subroutine coordinates_free
end module coordinates_module
