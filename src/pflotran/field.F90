module Field_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!

#include "petsc/finclude/petscvec.h"
  use petscvec
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: field_type

    !get material id
    ! 1 degree of freedom
    Vec :: porosity0
    Vec :: porosity_base_store
    Vec :: porosity_t
    Vec :: porosity_tpdt
    Vec :: porosity_geomech_store
    Vec :: tortuosity0

    Vec :: perm0_xx, perm0_yy, perm0_zz
    Vec :: perm0_xz, perm0_xy, perm0_yz

    Vec :: work, work_loc

    Vec :: volume0
    Vec :: compressibility0

    !TODO(geh): move these Vecs into their respective pms
    ! residual vectors
    Vec :: flow_r
    Vec :: tran_r

    ! Solution vectors (yy = previous solution, xx = current iterate)
    Vec :: flow_xx, flow_xx_loc, flow_dxx, flow_yy
    Vec :: tran_xx, tran_xx_loc, tran_dxx, tran_yy
    Vec :: flow_xxdot, flow_xxdot_loc
    Vec :: flow_rhs

    ! Accumulation terms (t = old time level, tpdt = new time level)
    Vec :: flow_accum_t, flow_accum_tpdt
    Vec :: tran_accum_t

    ! vectors for advanced nonlinear solvers other than Newton - Heeho
    Vec :: flow_scaled_xx, flow_work_loc

    Vec :: tran_log_xx, tran_work_loc

    ! mass transfer
    Vec :: flow_mass_transfer
    Vec :: tran_mass_transfer

    Vec :: flow_ts_mass_balance, flow_total_mass_balance
    Vec :: tran_ts_mass_balance, tran_total_mass_balance

    ! vector that holds the second layer of ghost cells for tvd
    Vec :: tvd_ghosts

    ! vectors to save temporally average quantities
    Vec, pointer :: avg_vars_vec(:)
    PetscInt :: nvars

    ! vectors to save temporally average flowrates
    Vec :: flowrate_inst
    Vec :: flowrate_aveg

    ! vectors to save velocity at face
    Vec :: vx_face_inst
    Vec :: vy_face_inst
    Vec :: vz_face_inst

    Vec, pointer :: max_change_vecs(:)

  end type field_type

  public :: FieldCreate, &
            FieldDestroy

contains

! ************************************************************************** !

function FieldCreate()
  !
  ! Allocates and initializes a new Field object
  !
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  !

  implicit none

  type(field_type), pointer :: FieldCreate

  type(field_type), pointer :: field

  allocate(field)

  ! nullify PetscVecs
  PetscObjectNullify(field%porosity0)
  PetscObjectNullify(field%porosity_base_store)
  PetscObjectNullify(field%porosity_geomech_store)
  PetscObjectNullify(field%porosity_t)
  PetscObjectNullify(field%porosity_tpdt)
  PetscObjectNullify(field%tortuosity0)

  PetscObjectNullify(field%perm0_xx)
  PetscObjectNullify(field%perm0_yy)
  PetscObjectNullify(field%perm0_zz)
  PetscObjectNullify(field%perm0_xy)
  PetscObjectNullify(field%perm0_xz)
  PetscObjectNullify(field%perm0_yz)

  PetscObjectNullify(field%work)
  PetscObjectNullify(field%work_loc)

  PetscObjectNullify(field%volume0)
  PetscObjectNullify(field%compressibility0)

  PetscObjectNullify(field%flow_r)
  PetscObjectNullify(field%flow_xx)
  PetscObjectNullify(field%flow_xx_loc)
  PetscObjectNullify(field%flow_scaled_xx)
  PetscObjectNullify(field%flow_work_loc)
  PetscObjectNullify(field%flow_dxx)
  PetscObjectNullify(field%flow_yy)
  PetscObjectNullify(field%flow_accum_t)
  PetscObjectNullify(field%flow_accum_tpdt)
  PetscObjectNullify(field%flow_xxdot)
  PetscObjectNullify(field%flow_xxdot_loc)
  PetscObjectNullify(field%flow_rhs)

  PetscObjectNullify(field%tran_r)
  PetscObjectNullify(field%tran_log_xx)
  PetscObjectNullify(field%tran_xx)
  PetscObjectNullify(field%tran_xx_loc)
  PetscObjectNullify(field%tran_dxx)
  PetscObjectNullify(field%tran_yy)
  PetscObjectNullify(field%tran_accum_t)
  PetscObjectNullify(field%tran_work_loc)

  PetscObjectNullify(field%tvd_ghosts)

  PetscObjectNullify(field%flow_mass_transfer)
  PetscObjectNullify(field%tran_mass_transfer)

  PetscObjectNullify(field%flow_ts_mass_balance)
  PetscObjectNullify(field%flow_total_mass_balance)
  PetscObjectNullify(field%tran_ts_mass_balance)
  PetscObjectNullify(field%tran_total_mass_balance)

  nullify(field%avg_vars_vec)
  field%nvars = 0

  PetscObjectNullify(field%flowrate_inst)
  PetscObjectNullify(field%flowrate_aveg)

  PetscObjectNullify(field%vx_face_inst)
  PetscObjectNullify(field%vy_face_inst)
  PetscObjectNullify(field%vz_face_inst)

  nullify(field%max_change_vecs)

  FieldCreate => field

end function FieldCreate

! ************************************************************************** !

subroutine FieldDestroy(field)
  !
  ! Deallocates a field object
  !
  ! Author: Glenn Hammond
  ! Date: 11/15/07

  use Petsc_Utility_module

  implicit none

  type(field_type), pointer :: field

  PetscErrorCode :: ierr
  PetscInt :: ivar
  PetscInt :: num_vecs

  if (.not.associated(field)) return

  ! Destroy PetscVecs
  call PUVecDestroy(field%porosity0)
  call PUVecDestroy(field%porosity_base_store)
  call PUVecDestroy(field%porosity_geomech_store)
  call PUVecDestroy(field%porosity_t)
  call PUVecDestroy(field%porosity_tpdt)
  call PUVecDestroy(field%tortuosity0)
  call PUVecDestroy(field%perm0_xx)
  call PUVecDestroy(field%perm0_yy)
  call PUVecDestroy(field%perm0_zz)
  call PUVecDestroy(field%perm0_xy)
  call PUVecDestroy(field%perm0_xz)
  call PUVecDestroy(field%perm0_yz)
  call PUVecDestroy(field%work)
  call PUVecDestroy(field%work_loc)
  call PUVecDestroy(field%volume0)
  call PUVecDestroy(field%compressibility0)
  call PUVecDestroy(field%flow_r)
  call PUVecDestroy(field%flow_xx)
  call PUVecDestroy(field%flow_xx_loc)
  call PUVecDestroy(field%flow_scaled_xx)
  call PUVecDestroy(field%flow_work_loc)
  call PUVecDestroy(field%flow_dxx)
  call PUVecDestroy(field%flow_yy)
  call PUVecDestroy(field%flow_accum_t)
  call PUVecDestroy(field%flow_accum_tpdt)
  call PUVecDestroy(field%flow_xxdot)
  call PUVecDestroy(field%flow_xxdot_loc)
  call PUVecDestroy(field%flow_rhs)
  call PUVecDestroy(field%tran_r)
  call PUVecDestroy(field%tran_log_xx)
  call PUVecDestroy(field%tran_xx)
  call PUVecDestroy(field%tran_xx_loc)
  call PUVecDestroy(field%tran_dxx)
  call PUVecDestroy(field%tran_yy)
  call PUVecDestroy(field%tran_accum_t)
  call PUVecDestroy(field%tran_work_loc)
  call PUVecDestroy(field%flow_mass_transfer)
  call PUVecDestroy(field%tran_mass_transfer)
  call PUVecDestroy(field%flow_ts_mass_balance)
  call PUVecDestroy(field%flow_total_mass_balance)
  call PUVecDestroy(field%tran_ts_mass_balance)
  call PUVecDestroy(field%tran_total_mass_balance)
  call PUVecDestroy(field%tvd_ghosts)
  do ivar = 1,field%nvars
    call PUVecDestroy(field%avg_vars_vec(ivar))
  enddo
  call PUVecDestroy(field%flowrate_inst)
  call PUVecDestroy(field%flowrate_aveg)
  call PUVecDestroy(field%vx_face_inst)
  call PUVecDestroy(field%vy_face_inst)
  call PUVecDestroy(field%vz_face_inst)
  if (associated(field%max_change_vecs)) then
    !geh: kludge as the compiler returns i4 in 64-bit
    num_vecs = size(field%max_change_vecs)
    call VecDestroyVecs(num_vecs,field%max_change_vecs,ierr);CHKERRQ(ierr)
  endif

  deallocate(field)
  nullify(field)

end subroutine FieldDestroy

end module Field_module
