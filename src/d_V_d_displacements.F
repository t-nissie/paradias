! d_V_d_displacements.F -*-f90-*-
! Time-stamp: <2013-03-16 16:48:27 t-nissie>
! Author: Takeshi NISHIMATSU
!!
#if defined HAVE_CONFIG_H
#  include "config.h"
#endif
#include "define.h"
subroutine d_V_d_displacements(n_atoms,L,c,qd,e_quadratic)
  use coordinates_module
  use quadratic_module
  implicit none
  !--- arguments -------------------------------------
  integer,                intent(in)    :: n_atoms
  integer,                intent(in)    :: L(3)
  type(coordinates_type), intent(inout) :: c
  type(quadratic_type),   intent(in)    :: qd
  real*8,                 intent(out)   :: e_quadratic
  !--- local -----------------------------------------
  integer alpha, i_atom, ix, iy, iz
  real*8  N_inv
  !--- END of variable definitions -------------------
  !call msg(UNIT_LOG, __FILE__, __LINE__, 'BEGIN', 'calc dVddi, calc e%dipo_dipo')
!$omp parallel do private(iz,iy,ix,alpha)
  do iz=0, L(3)-1
     do alpha = 1, 3
        do i_atom = 1, n_atoms
           do iy=0, L(2)-1
              do ix=0, L(1)/2
                 c%d_V_d_displacements_k(ix,iy,iz,alpha,i_atom) &
                      & = sum(         qd%C_tilde(ix,iy,iz,alpha,i_atom,:,:) &
                      &       * c%displacements_k(ix,iy,iz,             :,:) )
              end do
           end do
        end do
     end do
  end do
!$omp end parallel do

#if defined(HAVE_LIBFFTW3)
  do ix = 1, n_atoms*3
     call dfftw_execute(c%plans_d_V_d_displacements_k2R(ix))
  end do
#else
    do i_atom = 1, n_atoms
       do alpha = 1, 3
          call fft_3d_r2c_or_c2r_out_of_place(.false., L(1), L(2), L(3), &
               & c%d_V_d_displacements_R(0,0,0,alpha,i_atom), &
               & c%d_V_d_displacements_k(0,0,0,alpha,i_atom))
       end do
    end do
#endif

    N_inv = 1.0d0/(L(1)*L(2)*L(3))
!$omp parallel do
    do iz=0, L(3)-1
       do i_atom = 1, n_atoms
          do alpha = 1,3
             c%d_V_d_displacements_R(:,:,iz,alpha,i_atom) = &
           & c%d_V_d_displacements_R(:,:,iz,alpha,i_atom) * N_inv
          end do
       end do
    end do
!$omp end parallel do

    !write(UNIT_LOG,'(6f7.3)') c%d_V_d_displacements_R(0,0,0,:,:)


    !=== e%dipole_dipole = SUM( c%dipoR(:,:,:,:) * c%dVddi(:,:,:,:) ) / 2 ===========
    e_quadratic = 0.0d0
!$omp parallel do private(iz, i_atom, alpha, iy, ix) reduction(+:e_quadratic)
    do iz=0, L(3)-1
       do i_atom = 1,n_atoms
          do alpha = 1,3
             do iy=0, L(2)-1
                do ix=0, L(1)-1
                   e_quadratic = e_quadratic + c%displacements_R(ix,iy,iz,alpha,i_atom) &
                        &              * c%d_V_d_displacements_R(ix,iy,iz,alpha,i_atom)
                end do
             end do
          end do
       end do
  end do
!$omp end parallel do
  e_quadratic = e_quadratic / (2 * L(1) * L(2) * L(3) * n_atoms)
  !================================================================================

!   !call msg(UNIT_LOG, __FILE__, __LINE__, 'END', '')
end subroutine d_V_d_displacements
