  real*8 function kinetic_energy(xl)
    use Xtal_module
    use coordinates_module
    use Atom_module
    implicit none
    !--- arguments -------------------------------------
    type(Xtal_type),target,intent(in) :: xl
    !--- local -----------------------------------------
    real*8                  :: tmp_e
    integer                 :: alpha, ix,iy,iz
    !type(Atom_type),pointer :: p_atom_I
    integer                 :: I_atom


    tmp_e = 0.0d0
    I_atom = 1
       
!$omp parallel do private(iz, alpha, iy, ix) reduction(+:tmp_e)
    do iz=0, xl%L(3)-1

          do alpha = 1,3
             do iy=0, xl%L(2)-1
                do ix=0, xl%L(1)-1
                   tmp_e = tmp_e + xl%c%displacements_P(ix,iy,iz,alpha,I_atom)**2 !* p_atom_I%mass_inv
                end do
             end do
          end do

    end do
!$omp end parallel do
  kinetic_energy = tmp_e / 2
end function kinetic_energy
