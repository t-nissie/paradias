! coordinates_check.F -*-f90-*-
! Time-stamp: <2013-12-03 17:29:46 takeshi>
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
  integer(C_INTPTR_T), parameter :: n_atoms = 2
  integer(C_INTPTR_T), parameter :: Lx = 10
  integer(C_INTPTR_T), parameter :: Ly = 20
  integer(C_INTPTR_T), parameter :: Lz = 30
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

  c%displacements_R(:,:,:, :, :) = 0.1d0

  call fftw_mpi_execute_dft_r2c(c%plan_displacements_R2k, &
       & c%displacements_R, c%displacements_k)

  write(6,*) my_mpi_rank, c%displacements_k(1,1,0,0,0)
  write(6,*) my_mpi_rank, c%displacements_k(1,2,0,0,0)
  write(6,*) my_mpi_rank, c%displacements_k(1,3,0,0,0)
  write(6,*) my_mpi_rank, c%displacements_k(2,1,0,0,0)
  write(6,*) my_mpi_rank, c%displacements_k(2,2,0,0,0)
  write(6,*) my_mpi_rank, c%displacements_k(2,3,0,0,0)

  call MPI_FINALIZE(i_error)
end program coordinates_check
!Local variables:
!  compile-command: "gfortran -I/usr/local/include -Wall -ffree-form -c coordinates_module.F && gfortran -Wall -ffree-form -I/usr/local/include -c coordinates_check.F && mpif90 -o coordinates_check coordinates_check.o coordinates_module.o -lfftw3 -lfftw3_mpi && mpirun -np 2 ./coordinates_check && sleep 1"
!End:
