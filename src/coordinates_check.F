! coordinates_check.F -*-f90-*-
! Time-stamp: <2013-12-04 15:43:52 takeshi>
! Author: Takeshi NISHIMATSU
!!
#include "define.h"
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif

program coordinates_check
  use coordinates_module
  use, intrinsic :: iso_c_binding
  implicit none
  include 'mpif.h'
  real*8,              parameter :: check_double = 0.1d0, accuracy = 1.0d-14
  integer(C_INTPTR_T), parameter :: n_atoms = 2
  integer(C_INTPTR_T), parameter :: Lx = 10
  integer(C_INTPTR_T), parameter :: Ly = 20
  integer(C_INTPTR_T), parameter :: Lz = 30
  integer(C_INTPTR_T)            :: i_atom, alpha, ix, iy, iz

  type(coordinates_type) :: c
  integer  :: i_error, my_mpi_rank, n_mpi_procs, n_threads

# if defined(_OPENMP)
  integer                            :: provided, OMP_GET_MAX_THREADS
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED, provided, i_error)
  call dfftw_init_threads(i_error)
  n_threads = OMP_GET_MAX_THREADS()
  call dfftw_plan_with_nthreads(n_threads)
  write(0,'(a,i7,a,i4)') '_OPENMP =', _OPENMP, ',    OMP_GET_MAX_THREADS =', n_threads
# else
  call MPI_INIT(i_error)
  n_threads = 1
# endif
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_mpi_rank, i_error)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi_procs, i_error)

  call coordinates_allocate(n_atoms, Lx, Ly, Lz, c)

  ! Set a value for check
  c%displacements_R(:,:,:, :, :) = check_double

  ! r2c FFT
  call fftw_mpi_execute_dft_r2c(c%plan_displacements_R2k, &
       & c%displacements_R, c%displacements_k)

  ! Check results
  do i_atom = 1, n_atoms
     do alpha = 1, 3
        if (my_mpi_rank.eq.0) then
           if (abs(check_double*Lx*Ly*Lz-real(c%displacements_k(alpha,i_atom,0,0,0))) &
                & >check_double*Lx*Ly*Lz*accuracy) then
              write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
                   & ': Error in accuracy of r2c FFT.'
              call MPI_ABORT(MPI_COMM_WORLD, 1, i_error)
           end if
        else
           if (real(c%displacements_k(alpha,i_atom,0,0,0)) > accuracy) then
              write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
                   & ': Error in accuracy of r2c FFT.'
              call MPI_ABORT(MPI_COMM_WORLD, 2, i_error)
           end if
        end if
     end do
  end do

  ! Set values for check
  c%d_V_d_displacements_k(:,:,:, :, :) = (0.0d0, 0.0d0)
  if (my_mpi_rank.eq.0) then    !!! Set a value only at the Gamma point !!!
     c%d_V_d_displacements_k(:,:,0, 0, 0) = (check_double, 0.0d0)
  end if

  ! c2r FFT
  call fftw_mpi_execute_dft_c2r(c%plan_d_V_d_displacements_k2R, &
       & c%d_V_d_displacements_k, c%d_V_d_displacements_R)

  ! Check results
!$omp parallel do
  do iz = 0, c%local_Lz-1
     do iy = 0, Ly-1
        do ix = 0, Lx-1
           do i_atom = 1, n_atoms
              do alpha = 1, 3
                 if (abs(check_double-c%d_V_d_displacements_R(alpha,i_atom,ix,iy,iz)) &
                      & >check_double*accuracy) then
                    write(0,'(a,a,i3,a)')   __FILE__, ':', __LINE__ ,  &
                         & ': Error in accuracy of c2r FFT.'
                    call MPI_ABORT(MPI_COMM_WORLD, 3, i_error)
                 end if
              end do
           end do
        end do
     end do
  end do
!$omp end parallel do

  call coordinates_free(c)
  call MPI_FINALIZE(i_error)
end program coordinates_check
!Local variables:
!  compile-command: "gfortran -I/usr/local/include -Wall -ffree-form -c coordinates_module.F && gfortran -Wall -ffree-form -I/usr/local/include -c coordinates_check.F && mpif90 -o coordinates_check coordinates_check.o coordinates_module.o -lfftw3 -lfftw3_mpi && mpirun -np 2 ./coordinates_check && sleep 1"
!End:
