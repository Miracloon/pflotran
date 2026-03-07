module Geomechanics_Field_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

  implicit none

  private

  type, public :: geomech_field_type
    Vec :: work
    Vec :: work_loc
    ! residual vectors
    Vec :: rhs
    Vec :: press          ! store pressure from subsurf
    Vec :: press_loc
    Vec :: press_init_loc ! store initial pressure
    Vec :: temp           ! store temperature from subsurf
    Vec :: temp_loc
    Vec :: temp_init_loc  ! store initial temperature
    Vec :: subsurf_vec_1dof ! MPI
    Vec :: imech_loc
    Vec :: strain
    Vec :: strain_loc
    Vec :: stress
    Vec :: stress_loc
    Vec :: stress_total
    Vec :: stress_total_loc
    Vec :: strain_subsurf  ! Stores strains after scattering from geomech to subsurf
    Vec :: stress_subsurf  ! Stores stresses after scattering from geomech to subsurf
    Vec :: strain_subsurf_loc
    Vec :: stress_subsurf_loc

    Vec :: porosity ! store porosity from subsurf
    Vec :: porosity_loc
    Vec :: porosity_init_loc

    ! spatially-varying material properties
    Vec :: youngs_modulus_loc
    Vec :: poissons_ratio_loc
    Vec :: density_loc
    Vec :: biot_coeff_loc
    Vec :: thermal_exp_coeff_loc
    Vec :: body_force_x_loc
    Vec :: body_force_y_loc
    Vec :: body_force_z_loc

    Vec :: fluid_density ! store density from subsurf
    Vec :: fluid_density_loc
    Vec :: fluid_density_init_loc

    ! Solution vectors (xx = current iterate)
    Vec :: disp_xx, disp_xx_loc
    Vec :: disp_xx_init_loc

    Mat :: A ! stores the matrix of the linear system

  end type geomech_field_type

  public :: GeomechFieldCreate, &
            GeomechFieldDestroy

contains

! ************************************************************************** !

function GeomechFieldCreate()
  !
  ! Allocates and initializes a new geomechanics
  ! Field object
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/05/13
  !

  implicit none

  type(geomech_field_type), pointer :: GeomechFieldCreate
  type(geomech_field_type), pointer :: geomech_field

  allocate(geomech_field)

  ! nullify PetscVecs
  PetscObjectNullify(geomech_field%work)
  PetscObjectNullify(geomech_field%work_loc)

  PetscObjectNullify(geomech_field%rhs)
  PetscObjectNullify(geomech_field%disp_xx)
  PetscObjectNullify(geomech_field%disp_xx_loc)
  PetscObjectNullify(geomech_field%disp_xx_init_loc)

  PetscObjectNullify(geomech_field%press)
  PetscObjectNullify(geomech_field%press_loc)
  PetscObjectNullify(geomech_field%press_init_loc)
  PetscObjectNullify(geomech_field%temp)
  PetscObjectNullify(geomech_field%temp_loc)
  PetscObjectNullify(geomech_field%temp_init_loc)
  PetscObjectNullify(geomech_field%subsurf_vec_1dof)
  PetscObjectNullify(geomech_field%imech_loc)

  PetscObjectNullify(geomech_field%strain)
  PetscObjectNullify(geomech_field%strain_loc)
  PetscObjectNullify(geomech_field%stress)
  PetscObjectNullify(geomech_field%stress_loc)
  PetscObjectNullify(geomech_field%stress_total)
  PetscObjectNullify(geomech_field%stress_total_loc)

  PetscObjectNullify(geomech_field%strain_subsurf)
  PetscObjectNullify(geomech_field%stress_subsurf)
  PetscObjectNullify(geomech_field%strain_subsurf_loc)
  PetscObjectNullify(geomech_field%stress_subsurf_loc)

  PetscObjectNullify(geomech_field%porosity)
  PetscObjectNullify(geomech_field%porosity_loc)
  PetscObjectNullify(geomech_field%porosity_init_loc)

  PetscObjectNullify(geomech_field%youngs_modulus_loc)
  PetscObjectNullify(geomech_field%poissons_ratio_loc)
  PetscObjectNullify(geomech_field%density_loc)
  PetscObjectNullify(geomech_field%biot_coeff_loc)
  PetscObjectNullify(geomech_field%thermal_exp_coeff_loc)
  PetscObjectNullify(geomech_field%body_force_x_loc)
  PetscObjectNullify(geomech_field%body_force_y_loc)
  PetscObjectNullify(geomech_field%body_force_z_loc)

  PetscObjectNullify(geomech_field%fluid_density)
  PetscObjectNullify(geomech_field%fluid_density_loc)
  PetscObjectNullify(geomech_field%fluid_density_init_loc)

  PetscObjectNullify(geomech_field%A)

  GeomechFieldCreate => geomech_field

end function GeomechFieldCreate

! ************************************************************************** !

subroutine GeomechFieldDestroy(geomech_field)
  !
  ! Deallocates a geomechanics field object
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/05/13

  use Petsc_Utility_module

  implicit none

  type(geomech_field_type), pointer :: geomech_field

  if (.not.associated(geomech_field)) return

  ! Destroy PetscVecs
  call PUVecDestroy(geomech_field%work)
  call PUVecDestroy(geomech_field%work_loc)
  call PUVecDestroy(geomech_field%rhs)
  call PUVecDestroy(geomech_field%disp_xx)
  call PUVecDestroy(geomech_field%disp_xx_loc)
  call PUVecDestroy(geomech_field%disp_xx_init_loc)
  call PUVecDestroy(geomech_field%press)
  call PUVecDestroy(geomech_field%press_loc)
  call PUVecDestroy(geomech_field%press_init_loc)
  call PUVecDestroy(geomech_field%temp)
  call PUVecDestroy(geomech_field%temp_loc)
  call PUVecDestroy(geomech_field%temp_init_loc)
  call PUVecDestroy(geomech_field%subsurf_vec_1dof)
  call PUVecDestroy(geomech_field%imech_loc)
  call PUVecDestroy(geomech_field%strain)
  call PUVecDestroy(geomech_field%strain_loc)
  call PUVecDestroy(geomech_field%stress)
  call PUVecDestroy(geomech_field%stress_loc)
  call PUVecDestroy(geomech_field%stress_total)
  call PUVecDestroy(geomech_field%stress_total_loc)
  call PUVecDestroy(geomech_field%strain_subsurf)
  call PUVecDestroy(geomech_field%stress_subsurf)
  call PUVecDestroy(geomech_field%strain_subsurf_loc)
  call PUVecDestroy(geomech_field%stress_subsurf_loc)
  call PUVecDestroy(geomech_field%porosity)
  call PUVecDestroy(geomech_field%porosity_loc)
  call PUVecDestroy(geomech_field%porosity_init_loc)
  call PUVecDestroy(geomech_field%youngs_modulus_loc)
  call PUVecDestroy(geomech_field%poissons_ratio_loc)
  call PUVecDestroy(geomech_field%density_loc)
  call PUVecDestroy(geomech_field%biot_coeff_loc)
  call PUVecDestroy(geomech_field%thermal_exp_coeff_loc)
  call PUVecDestroy(geomech_field%body_force_x_loc)
  call PUVecDestroy(geomech_field%body_force_y_loc)
  call PUVecDestroy(geomech_field%body_force_z_loc)
  call PUVecDestroy(geomech_field%fluid_density)
  call PUVecDestroy(geomech_field%fluid_density_loc)
  call PUVecDestroy(geomech_field%fluid_density_init_loc)

  PetscObjectNullify(geomech_field%A)

  if (associated(geomech_field)) deallocate(geomech_field)
  nullify(geomech_field)

end subroutine GeomechFieldDestroy

end module Geomechanics_Field_module
