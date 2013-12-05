! ifc_matrix.F -*-f90-*-
! Time-stamp: <2013-12-05 14:05:56 takeshi>
! Author: Takeshi NISHIMATSU
! Reference: Xavier Gonze and Changyol Lee:
!            PRB Vol.55, pp.10355-10368 (1997),
!            See equation (72) and (76).
!!
#include "define.h"
subroutine ifc_matrix(n_atoms, Lx, Ly, Lz, c)
  use, intrinsic :: iso_c_binding
  use :: coordinates_module
  implicit none
  include 'mpif.h'
  integer(C_INTPTR_T), intent(in)       :: n_atoms
  integer(C_INTPTR_T), intent(in)       :: Lx, Ly, Lz
  type(coordinates_type), intent(inout) :: c
  integer(C_INTPTR_T)                   :: alpha !i_atom, j_atom, alpha, beta, ix, iy, iz

  !--- BEGIN sc only ------------------------------------------
  c%C_tilde_short_R(:,:,:, :,:, :,:) = 0.0d0
  if (c%local_z_offset.eq.0) then
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1, Lx-1,0,0) = -1.0d0; end do
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1,    1,0,0) = -1.0d0; end do
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1, 0,Ly-1,0) = -1.0d0; end do
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1, 0,   1,0) = -1.0d0; end do
  end if
  if (c%local_z_offset<=Lz-1 .and. Lz-1<c%local_z_offset+c%local_Lz) then
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1, 0,0,Lz-1-c%local_z_offset) = -1.0d0; end do
  end if
  if (c%local_z_offset<=1    .and.    1<c%local_z_offset+c%local_Lz) then
     do alpha=1,3; c%C_tilde_short_R(alpha,1, alpha,1, 0,0,   1-c%local_z_offset) = -1.0d0; end do
  end if
  !--- END   sc only ------------------------------------------

  !Preparation of C_tilde_short_k
  call fftw_mpi_execute_dft_r2c(c%plan_C_tilde_short_R2k, &
       & c%C_tilde_short_R, c%C_tilde_short_k)




  ! call dipole_dipole_3d_Ewald(xl, xl%qd)

  ! !--- BEGIN of C_tilde = C_tilde_short + Z C_bar_DD Z --------
  ! p_atom_J => xl%atom_list
  ! J = 0
  ! do while(atom_list_each_with_index(p_atom_J,J))
  !    p_atom_I => xl%atom_list
  !    I = 0
  !    do while(atom_list_each_with_index(p_atom_I,I))

  !       do iz = 0, xl%L(3)-1
  !          do iy = 0, xl%L(2)-1
  !             do ix = 0, xl%L(1)/2
  !                xl%qd%C_tilde(ix,iy,iz, :,I, :,J) = &
  !                     &                xl%qd%C_tilde_short(ix,iy,iz, :,I, :,J)   &
  !                     & + matmul( matmul(p_atom_I%Z_tensor(:,:),              &
  !                     &                  xl%qd%C_bar_DD(ix,iy,iz, :,I, :,J)), &
  !                     &                  p_atom_J%Z_tensor(:,:) )
  !                if (I.eq.J) then
  !                   do alpha = 1, 3
  !                      xl%qd%C_tilde(ix,iy,iz, alpha,I, alpha,J) = &
  !                    & xl%qd%C_tilde(ix,iy,iz, alpha,I, alpha,J) + 2*xl%qd%elevation
  !                   end do
  !                end if
  !                end do
  !          end do
  !       end do

  !    end do
  ! end do
  ! !--- END of C_tilde = C_tilde_short + Z C_bar_DD Z ----------
  !  !do alpha = 1,3
  !  !   write(UNIT_LOG,'(3(f7.2,f6.2))') xl%qd%C_tilde(0,0,0, alpha,1, :,1)
  !  !end do

  ! !
  ! ! C_tilde should be weighted with effective mass around here or in the loop above.
  ! !
end subroutine ifc_matrix
