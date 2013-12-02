! dipole_dipole_3d_Ewald.F -*-f90-*-
! Time-stamp: <2013-03-16 12:39:02 t-nissie>
! Author: Takeshi NISHIMATSU
! Reference: Xavier Gonze and Changyol Lee: PRB Vol.55, pp.10355-10368 (1997)
!!
#include "define.h"
subroutine dipole_dipole_3d_Ewald(xl, qd)
  use Xtal_module
  use Dipole_Dipole_Functions_module
  implicit none
  !--- arguments -------------------------------------
  type(Xtal_type), target, intent(inout)  :: xl
  type(quadratic_type), intent(inout)  :: qd
  !--- local -----------------------------------------
  real*8 epsilon_inf_det
  real*8 epsilon_inf_inv(3,3)
  real*8 kappa   !An arbitrary parameter for Ewald sum.
  real*8 tmp(0:xl%L(1)-1, &
       &     0:xl%L(2)-1, &
       &     0:xl%L(3)-1, &
       &     3,xl%n_atoms,   &
       &     3,xl%n_atoms)     ! auto allocation
  integer I, J
  integer ix, iy, iz
  integer nx, ny, nz
  real*8  r(3), rRn(3), rRn_Cartesian(3)
  real*8  k(3),  kG(3),  kG_Cartesian(3), kGekG
  real*8  Delta(3), D2, y2, y, xx_y2(3,3), kGkG_kGekG(3,3)
  type(Atom_type),pointer :: p_atom_I, p_atom_J
  !character (len=50) fmt
  !--- END of variable definitions ------------------- 

  epsilon_inf_det = determinant3x3(xl%epsilon_inf(:,:))
  !write(UNIT_LOG,'(a, f22.16)') 'epsilon_inf_det =', epsilon_inf_det

  epsilon_inf_inv(:,:) = inv3x3(xl%epsilon_inf(:,:))

  kappa = ewald_kappa(epsilon_inf_det, xl%volume, xl%L)
  write(UNIT_LOG,'(a,f22.16)')  'kappa = ', kappa

  tmp(:,:,:, :,:, :,:) =  0.0d0

  !--- BEGIN of self-interaction correction ----------------
  do J = 1, xl%n_atoms
     tmp(0, 0, 0, :,J, :,J) = &
          & - 2 * M_2_SQRTPI / 3 * epsilon_inf_inv(:,:)
  end do
  !--- END of self-interaction correction ------------------

  !--- BEGIN of real-space summation -----------------------
  p_atom_J => xl%atom_list
  J = 0
  do while(atom_list_each_with_index(p_atom_J,J))
     p_atom_I => xl%atom_list
     I = 0
     do while(atom_list_each_with_index(p_atom_I,I))
        r(:) = p_atom_I%r(:) - p_atom_J%r(:)
        !write(UNIT_LOG,'(2i2,2(f7.3,2f8.3))') I, J, p_atom_I%r(:),  p_atom_J%r(:)

        do iz = 0, xl%L(3)-1
           do iy = 0, xl%L(2)-1
              do ix = 0, xl%L(1)-1

                 do nz = -qd%n_max(3)-1, qd%n_max(3)
                    rRn(3) = r(3) + iz + xl%L(3)*nz
                    do ny = -qd%n_max(2)-1, qd%n_max(2)
                       rRn(2) = r(2) + iy + xl%L(2)*ny
                       do nx = -qd%n_max(1)-1, qd%n_max(1)
                          if (I.eq.J .and.                             &
                               &  ix.eq.0 .and. iy.eq.0 .and. iz.eq.0 .and. &
                               &  nx.eq.0 .and. ny.eq.0 .and. nz.eq.0) cycle
                          rRn(1) = r(1) + ix + xl%L(1)*nx
                          rRn_Cartesian(:) = matmul(xl%prim(:,:),rRn(:))
                          Delta(:) = matmul(epsilon_inf_inv(:,:),rRn_Cartesian(:))
                          D2 = dot_product(rRn_Cartesian(:),Delta(:))
                          y2 = kappa**2 * D2
                          y = sqrt(y2)
                          xx_y2(:,:) = dyad(Delta(:),Delta(:)) / D2
                          tmp(ix,iy,iz, :,I, :,J) = &
                          tmp(ix,iy,iz, :,I, :,J) & !You need this line!
                               & + epsilon_inf_inv(:,:) * decay1(y,y2) &
                               & - xx_y2(:,:)           * decay2(y,y2)
                       end do
                    end do
                 end do

              end do
           end do
        end do

     end do
  end do

  tmp = kappa**3 * tmp / sqrt(epsilon_inf_det)   !Do NOT forget it!
  !--- END of real-space summation -------------------------

  !--- BEGIN of FFT ----------------------------------------
  call fft_3d_r2c_NxN_lower_then_copy(xl%L(1),xl%L(2),xl%L(3), 3*xl%n_atoms, tmp, qd%C_bar_DD)
  !--- END of FFT ------------------------------------------ 

  !--- BEGIN of reciplocal-space summation -----------------
  p_atom_J => xl%atom_list
  J = 0
  do while(atom_list_each_with_index(p_atom_J,J))
     p_atom_I => xl%atom_list
     I = 0
     do while(atom_list_each_with_index(p_atom_I,I))
        r(:) = p_atom_I%r(:) - p_atom_J%r(:)

        do iz = 0, xl%L(3)-1
           k(3) = dble(iz)/xl%L(3)
           do iy = 0, xl%L(2)-1
              k(2) = dble(iy)/xl%L(2)
              do ix = 0, xl%L(1)/2
                 k(1) = dble(ix)/xl%L(1)

                 do nz = -qd%G_max(3)-1, qd%G_max(3)
                    kG(3) = 2*M_PI*(k(3)+nz)
                    do ny = -qd%G_max(2)-1, qd%G_max(2)
                       kG(2) = 2*M_PI*(k(2)+ny)
                       do nx = -qd%G_max(1)-1, qd%G_max(1)
                          if (ix.eq.0 .and. iy.eq.0 .and. iz.eq.0 .and. &
                            & nx.eq.0 .and. ny.eq.0 .and. nz.eq.0) cycle
                          kG(1) = 2*M_PI*(k(1)+nx)
                          kG_Cartesian(:) = matmul(kG(:),xl%r_prim(:,:))
                          kGekG = dot_product(              kG_Cartesian(:), &
                               & matmul(xl%epsilon_inf(:,:),kG_Cartesian(:)))
                          kGkG_kGekG(:,:) = dyad(kG_Cartesian(:), &
                               &                 kG_Cartesian(:)) / kGekG
                          qd%C_bar_DD(ix,iy,iz, :,I, :,J) = &
                        & qd%C_bar_DD(ix,iy,iz, :,I, :,J)   &
                        & + 4 * M_PI / xl%volume * kGkG_kGekG(:,:) &
                        &   * exp( cmplx(-kGekG/4/kappa**2, dot_product(kG,r)) )
                       end do
                    end do
                 end do

              end do
           end do
        end do

     end do
  end do
  !--- END of reciplocal-space summation -------------------

end subroutine dipole_dipole_3d_Ewald
