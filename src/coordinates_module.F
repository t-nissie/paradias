! coordinates_module.F -*-f90-*-
! Time-stamp: <2013-03-31 22:02:10 takeshi>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
module coordinates_module
  implicit none
  type coordinates_type
     real*8,     allocatable ::       displacements_P(:,:,:, :, :)   ! momenta

     real*8,     allocatable ::       displacements_R(:,:,:, :, :)
     complex*16, allocatable ::       displacements_k(:,:,:, :, :)
     integer*8,  allocatable :: plans_displacements_R2k(:)

     complex*16, allocatable ::       d_V_d_displacements_k(:,:,:, :, :)   ! t for temporally
     real*8,     allocatable ::       d_V_d_displacements_R(:,:,:, :, :)
     integer*8,  allocatable :: plans_d_V_d_displacements_k2R(:)
  end type coordinates_type

contains
  subroutine coordinates_allocate(n_atoms, L, c)
    implicit none
    integer,                intent(in)    :: n_atoms
    integer,                intent(in)    :: L(3)
    type(coordinates_type), intent(inout) :: c

    integer i

#   if defined(HAVE_LIBFFTW3)
    integer alpha, i_atom
#   include "fftw3.f"
#   endif

    if (allocated(c%displacements_P)) deallocate(c%displacements_P)
    allocate(     c%displacements_P(0:L(1)-1, 0:L(2)-1, 0:L(3)-1, 3, n_atoms))
!$omp parallel do
    do i = 0,L(3)-1
       c%displacements_P(:,:,i,:,:) = 0.0d0
    end do
!$omp end parallel do

    if (allocated(c%displacements_R)) deallocate(c%displacements_R)
    allocate(     c%displacements_R(0:L(1)-1, 0:L(2)-1, 0:L(3)-1, 3, n_atoms))
    if (allocated(c%displacements_k)) deallocate(c%displacements_k)
    allocate(     c%displacements_k(0:L(1)/2, 0:L(2)-1, 0:L(3)-1, 3, n_atoms)) !"/2" See fft.F!
!$omp parallel do
    do i = 0,L(3)-1
       c%displacements_R(:,:,i,:,:) = 0.0d0
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0,L(3)-1
       c%displacements_k(:,:,i,:,:) = (0.0d0,0.0d0)
    end do
!$omp end parallel do
#  if defined(HAVE_LIBFFTW3)
    if (allocated(c%plans_displacements_R2k)) deallocate(c%plans_displacements_R2k)
    allocate(     c%plans_displacements_R2k(3*n_atoms))
    i = 0
    do i_atom = 1, n_atoms
       do alpha = 1, 3
          i = i + 1
          call dfftw_plan_dft_r2c_3d(c%plans_displacements_R2k(i), L(1), L(2), L(3), &
               & c%displacements_R(0,0,0,alpha,i_atom), &
               & c%displacements_k(0,0,0,alpha,i_atom), FFTW_MEASURE)
       end do
    end do
#  endif

    if (allocated(c%d_V_d_displacements_k)) deallocate(c%d_V_d_displacements_k)
    allocate(     c%d_V_d_displacements_k(0:L(1)/2, 0:L(2)-1, 0:L(3)-1, 3, n_atoms)) !"/2" See fft.F!
    if (allocated(c%d_V_d_displacements_R)) deallocate(c%d_V_d_displacements_R)
    allocate(     c%d_V_d_displacements_R(0:L(1)-1, 0:L(2)-1, 0:L(3)-1, 3, n_atoms))
!$omp parallel do
    do i = 0,L(3)-1
       c%d_V_d_displacements_k(:,:,i,:,:) = (0.0d0,0.0d0)
    end do
!$omp end parallel do
!$omp parallel do
    do i = 0,L(3)-1
       c%d_V_d_displacements_R(:,:,i,:,:) = 0.0d0
    end do
!$omp end parallel do
#  if defined(HAVE_LIBFFTW3)
    if (allocated(c%plans_d_V_d_displacements_k2R)) deallocate(c%plans_d_V_d_displacements_k2R)
    allocate(     c%plans_d_V_d_displacements_k2R(3*n_atoms))
    i = 0
    do i_atom = 1, n_atoms
       do alpha = 1, 3
          i = i + 1
          call dfftw_plan_dft_c2r_3d(c%plans_d_V_d_displacements_k2R(i), L(1), L(2), L(3), &
               & c%d_V_d_displacements_k(0,0,0,alpha,i_atom), &
               & c%d_V_d_displacements_R(0,0,0,alpha,i_atom), FFTW_MEASURE)
       end do
    end do
#  endif
  end subroutine coordinates_allocate

  subroutine displacements_FFT_R2k(n_atoms,L,c)
    implicit none
    integer,                 intent(in)    :: n_atoms
    integer,                 intent(in)    :: L(3)
    type(coordinates_type),  intent(inout) :: c
#  if defined(HAVE_LIBFFTW3)
    integer i
    do i = 1, n_atoms*3
       call dfftw_execute(c%plans_displacements_R2k(i))
    end do
#  else
    integer i_atom, alpha
    do i_atom = 1, n_atoms
       do alpha = 1, 3
          call fft_3d_r2c_or_c2r_out_of_place(.true., L(1), L(2), L(3), &
               & c%displacements_R(0,0,0,alpha,i_atom), &
               & c%displacements_k(0,0,0,alpha,i_atom))
       end do
    end do
#  endif
    !write(UNIT_LOG,*) c%displacements_R(0,0,0,3,1)
    !write(UNIT_LOG,*) c%displacements_k(0,0,0,3,1)
    !write(UNIT_LOG,*) c%displacements_R(0,0,0,3,2)
    !write(UNIT_LOG,*) c%displacements_k(0,0,0,3,2)
  end subroutine displacements_FFT_R2k

  subroutine displacements_initialize(n_atoms,L,init_avg,init_dev,c,restart_fn)
    use marsaglia_tsang_uni64_module
    implicit none
    integer,                 intent(in)    :: n_atoms
    integer,                 intent(in)    :: L(3)
    real*8,                  intent(in)    :: init_avg(3),init_dev(3)
    type(coordinates_type),  intent(inout) :: c
    character(len=*),        intent(in)    :: restart_fn
    integer ios, ix, iy, iz, I_atom, alpha, dmy1, dmy2, dmy3
    real*8 dmy

    open(unit=UNIT_RESTART, file=restart_fn, status='old', action='read', iostat=ios)
    if (ios.eq.0) then
       do iz = 0, L(3)-1
          do iy = 0, L(2)-1
             do ix = 0, L(1)-1
                read(UNIT_RESTART,*) dmy1, dmy2, dmy3, c%displacements_R(ix,iy,iz,:,1), c%displacements_P(ix,iy,iz,:,1)
             end do
          end do
       end do
    else
       dmy = uni64(123456789,987654321) ! (xl%seed(1), xl%seed(2))
       do I_atom = 1, n_atoms
          do alpha = 1, 3
             do iz = 0, L(3)-1
                do iy = 0, L(2)-1
                   do ix = 0, L(1)-1
                      c%displacements_R(ix,iy,iz,alpha,I_atom) = normal_dist(init_avg(alpha),init_dev(alpha))
                   end do
                end do
             end do
          end do
       end do
    end if
    close(UNIT_RESTART)
   end subroutine displacements_initialize

  subroutine coordinates_write(n_atoms,L,c,fn)
    implicit none
    integer,                intent(in) :: n_atoms
    integer,                intent(in) :: L(3)
    type(coordinates_type), intent(in) :: c
    character(len=*) :: fn
    character(len=FILENAME_LEN) :: coord_directory_fn
    integer ix,iy,iz,i !,alpha

    ! if (coordinates_directory.eq.'never') then
    !    call msg(UNIT_LOG, __FILE__, __LINE__, 'END', 'never write .coord files.')
    !    return
    ! else if (coord_directory.eq.'') then
    !    coord_directory_fn = trim(fn)
    ! else
    !    coord_directory_fn = trim(coord_directory) // '/' // trim(fn)
    ! end if

    coord_directory_fn = fn

    open(unit=UNIT_COORDINATES_WRITE, file=coord_directory_fn, status='replace')
    do iz = 0, L(3)-1
       do iy = 0, L(2)-1
          do ix = 0, L(1)-1
                write(UNIT_COORDINATES_WRITE,'(i4,2i5,3f11.6,3e14.6)') ix, iy, iz, &
                     & (c%displacements_R(ix,iy,iz,i,1),i=1,3), & !!! SINGLE ATOM !!!
                     & (c%displacements_P(ix,iy,iz,i,1),i=1,3)    !!! SINGLE ATOM !!!
          end do
       end do
    end do
    close(UNIT_COORDINATES_WRITE)
  end subroutine coordinates_write
end module coordinates_module
