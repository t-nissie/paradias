! sort_eigen.F -*-f90-*-
! Time-stamp: <2013-02-13 11:54:19 takeshi>
! Author: Takeshi NISHIMATSU
! Purpose: pretty plotting
!!
#include "define.h"
subroutine sort_eigen(xl, plot_length, not_reset, new_eigenvalues, &
     & last_eigenvalues, last_eigenvectors, new_eigenvectors)
  use Xtal_module
  implicit none
  type(Xtal_type), intent(in)  :: xl
  real*8,           intent(in)  :: plot_length
  integer,          intent(in)  :: not_reset
  real*8,           intent(in)  :: new_eigenvalues(3*xl%n_atoms)
  real*8,           intent(out) :: last_eigenvalues(3*xl%n_atoms)
  complex*16,     intent(inout) :: last_eigenvectors(3*xl%n_atoms, 3*xl%n_atoms)
  complex*16,     intent(in)    ::  new_eigenvectors(3*xl%n_atoms, 3*xl%n_atoms)
  real*8  abs_similarity(3*xl%n_atoms, 3*xl%n_atoms)!|<new_eigenvector_i|last_eigenvector_j>|
  logical mask(3*xl%n_atoms)   ! auto allocation
  integer i(1), j

  if ( not_reset.eq.0 .or. plot_length.eq.0.0d0 ) then
     last_eigenvalues(:)    = new_eigenvalues(:)
     last_eigenvectors(:,:) = new_eigenvectors(:,:)
  else                                                    ! Successive k-points
     abs_similarity = abs(matmul(conjg(transpose(new_eigenvectors)),last_eigenvectors))
     mask(:) = .true.
     do j = 1, 3*xl%n_atoms
        i(:) = maxloc( abs_similarity(:,j), mask(:) )
        mask(i(1)) = .false.
        last_eigenvalues(j) = new_eigenvalues(i(1))
        last_eigenvectors(:,j) = new_eigenvectors(:,i(1))
     end do
  end if
end subroutine sort_eigen
