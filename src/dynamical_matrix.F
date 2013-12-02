! dynamical_matrix.F -*-f90-*-
! Time-stamp: <2013-03-16 12:39:02 t-nissie>
! Author: Takeshi NISHIMATSU
! Reference: Xavier Gonze and Changyol Lee:
!            PRB Vol.55, pp.10355-10368 (1997),
!            See equation (72) and (76).
!!
#include "define.h"
subroutine dynamical_matrix(xl)
  use Xtal_module
  implicit none
  type(Xtal_type), target, intent(inout) :: xl
  integer I, J
  integer ix, iy, iz, alpha
  type(Atom_type),pointer :: p_atom_I, p_atom_J
  
  real*8, allocatable :: C_short(:,:,:, :,:, :,:)
  allocate(C_short(0:xl%L(1)-1, &
       &           0:xl%L(2)-1, &
       &           0:xl%L(3)-1, &
       & 3,xl%n_atoms, &
       & 3,xl%n_atoms) )

  !--- begin sc only ------------------------------------------
  C_short(:,:,:, :,:, :,:) = 0.0d0
  do alpha=1,3; C_short(xl%L(1)-1,0,0, alpha,1, alpha,1) = -1.0d0; end do
  do alpha=1,3; C_short(        1,0,0, alpha,1, alpha,1) = -1.0d0; end do
  do alpha=1,3; C_short(0,xl%L(2)-1,0, alpha,1, alpha,1) = -1.0d0; end do
  do alpha=1,3; C_short(0,        1,0, alpha,1, alpha,1) = -1.0d0; end do
  do alpha=1,3; C_short(0,0,xl%L(3)-1, alpha,1, alpha,1) = -1.0d0; end do
  do alpha=1,3; C_short(0,0,        1, alpha,1, alpha,1) = -1.0d0; end do
  !--- end sc only --------------------------------------------

  !--- BEGIN of preparation of  C_tilde_short -----------------
  call fft_3d_r2c_NxN_lower_then_copy(xl%L(1),xl%L(2),xl%L(3), 3*xl%n_atoms, C_short, xl%qd%C_tilde_short)
   !do alpha = 1,3
   !   write(UNIT_LOG,'(3(f7.2,f6.2))') xl%qd%C_tilde_short(0,0,0, alpha,1, :,1)
   !end do
  !--- END of preparation of  C_tilde_short -------------------

  call dipole_dipole_3d_Ewald(xl, xl%qd)

  !--- BEGIN of C_tilde = C_tilde_short + Z C_bar_DD Z --------
  p_atom_J => xl%atom_list
  J = 0
  do while(atom_list_each_with_index(p_atom_J,J))
     p_atom_I => xl%atom_list
     I = 0
     do while(atom_list_each_with_index(p_atom_I,I))

        do iz = 0, xl%L(3)-1
           do iy = 0, xl%L(2)-1
              do ix = 0, xl%L(1)/2
                 xl%qd%C_tilde(ix,iy,iz, :,I, :,J) = &
                      &                xl%qd%C_tilde_short(ix,iy,iz, :,I, :,J)   &
                      & + matmul( matmul(p_atom_I%Z_tensor(:,:),              &
                      &                  xl%qd%C_bar_DD(ix,iy,iz, :,I, :,J)), &
                      &                  p_atom_J%Z_tensor(:,:) )
                 if (I.eq.J) then
                    do alpha = 1, 3
                       xl%qd%C_tilde(ix,iy,iz, alpha,I, alpha,J) = &
                     & xl%qd%C_tilde(ix,iy,iz, alpha,I, alpha,J) + 2*xl%qd%elevation
                    end do
                 end if
                 end do
           end do
        end do

     end do
  end do
  !--- END of C_tilde = C_tilde_short + Z C_bar_DD Z ----------
   !do alpha = 1,3
   !   write(UNIT_LOG,'(3(f7.2,f6.2))') xl%qd%C_tilde(0,0,0, alpha,1, :,1)
   !end do

  !
  ! C_tilde should be weighted with effective mass around here or in the loop above.
  !
end subroutine dynamical_matrix
