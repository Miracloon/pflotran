module Immiscible_module

#include "petsc/finclude/petscsnes.h"
  use petscsnes
  use Immiscible_Aux_module
  use Immiscible_Common_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: ImmiscibleSetup, &
            ImmiscibleInitializeTimestep, &
            ImmiscibleUpdateSolution, &
            ImmiscibleTimeCut, &
            ImmiscibleUpdateAuxVars, &
            ImmiscibleUpdateFixedAccum, &
            ImmiscibleComputeMassBalance, &
            ImmiscibleZeroMassBalanceDelta, &
            ImmiscibleResidual, &
            ImmiscibleSetPlotVariables, &
            ImmiscibleMapBCAuxVarsToGlobal, &
            ImmiscibleDestroy

contains

! ************************************************************************** !

subroutine ImmiscibleSetup(realization)
  !
  ! Creates arrays for auxiliary variables
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Fluid_module
  use Grid_module
  use Material_Aux_module
  use Output_Aux_module
  use Characteristic_Curves_module
  use Matrix_Zeroing_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type),pointer :: patch
  type(grid_type), pointer :: grid
  type(output_variable_list_type), pointer :: list
  type(material_parameter_type), pointer :: material_parameter
  type(immiscible_parameter_type), pointer :: immiscible_parameter

  PetscInt :: ghosted_id, iconn, sum_connection, local_id
  PetscBool :: error_found
  PetscInt :: flag(3)
  PetscBool, allocatable :: dof_is_active(:)
  PetscInt :: temp_int, idof
  PetscErrorCode :: ierr
                                                ! extra index for derivatives
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_bc(:)
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  option => realization%option
  patch => realization%patch
  grid => patch%grid

  patch%aux%Immiscible => ImmiscibleAuxCreate(option)
  immiscible_parameter => patch%aux%Immiscible%immiscible_parameter

  immis_numerical_derivatives = option%flow%numerical_derivatives
  if (immis_numerical_derivatives) then
    allocate(immis_min_pert(IMMIS_MAX_DOF))
    immis_min_pert(IMMIS_GAS_PRESSURE_DOF) = immis_pres_min_pert
    immis_min_pert(IMMIS_GAS_SATURATION_DOF) = immis_sat_min_pert
  endif

  ! ensure that material properties specific to this module are properly
  ! initialized
  material_parameter => patch%aux%Material%material_parameter
  error_found = PETSC_FALSE

  material_auxvars => patch%aux%Material%auxvars
  flag = 0
  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    if (patch%imat(ghosted_id) <= 0) cycle
    if (material_auxvars(ghosted_id)%volume < 0.d0 .and. flag(1) == 0) then
      flag(1) = 1
      option%io_buffer = 'ERROR: Non-initialized cell volume.'
      call PrintMsgByRank(option)
    endif
    if (material_auxvars(ghosted_id)%porosity_base < 0.d0 .and. &
        flag(2) == 0) then
      flag(2) = 1
      option%io_buffer = 'ERROR: Non-initialized porosity.'
      call PrintMsgByRank(option)
    endif
    if (minval(material_auxvars(ghosted_id)%permeability) < 0.d0 .and. &
        flag(3) == 0) then
      option%io_buffer = 'ERROR: Non-initialized permeability.'
      call PrintMsgByRank(option)
      flag(3) = 1
    endif
  enddo

  error_found = error_found .or. (maxval(flag) > 0)
  call MPI_Allreduce(MPI_IN_PLACE,error_found,ONE_INTEGER_MPI,MPI_C_BOOL, &
                     MPI_LOR,option%mycomm,ierr);CHKERRQ(ierr)
  if (error_found) then
    option%io_buffer = 'Material property errors found in ImmiscibleSetup.'
    call PrintErrMsg(option)
  endif

  temp_int = 0
  if (immis_numerical_derivatives) then
    temp_int = option%nflowdof
  endif
  allocate(immiscible_auxvars(0:temp_int,grid%ngmax))
  do ghosted_id = 1, grid%ngmax
    do idof = 0, temp_int
      call ImmiscibleAuxVarInit(immiscible_auxvars(idof,ghosted_id),option)
    enddo
  enddo
  patch%aux%Immiscible%auxvars => immiscible_auxvars
  patch%aux%Immiscible%num_aux = grid%ngmax

  ! count the number of boundary connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%boundary_condition_list)
  if (sum_connection > 0) then
    allocate(immiscible_auxvars_bc(sum_connection))
    do iconn = 1, sum_connection
      call ImmiscibleAuxVarInit(immiscible_auxvars_bc(iconn),option)
    enddo
    patch%aux%Immiscible%auxvars_bc => immiscible_auxvars_bc
  endif
  patch%aux%Immiscible%num_aux_bc = sum_connection

  ! count the number of source/sink connections and allocate
  ! auxvar data structures for them
  sum_connection = CouplerGetNumConnectionsInList(patch%source_sink_list)
  if (sum_connection > 0) then
    allocate(immiscible_auxvars_ss(sum_connection))
    do iconn = 1, sum_connection
      call ImmiscibleAuxVarInit(immiscible_auxvars_ss(iconn),option)
    enddo
    patch%aux%Immiscible%auxvars_ss => immiscible_auxvars_ss
  endif
  patch%aux%Immiscible%num_aux_ss = sum_connection

  list => realization%output_option%output_snap_variable_list
  call ImmiscibleSetPlotVariables(realization,list)
  list => realization%output_option%output_obs_variable_list
  call ImmiscibleSetPlotVariables(realization,list)

  XXFlux => ImmiscibleFlux
  XXBCFlux => ImmiscibleBCFlux

  if (Initialized(immis_debug_cell_id) .and. &
      option%comm%size > 1) then
    option%io_buffer = 'Cannot debug cells in parallel.'
    call PrintErrMsg(option)
  endif

  allocate(dof_is_active(option%nflowdof))
  dof_is_active = PETSC_TRUE
  call PatchCreateZeroArray(patch,dof_is_active, &
                            patch%aux%Immiscible%matrix_zeroing,option)
  deallocate(dof_is_active)

  immis_ts_count = 0
  immis_ts_cut_count = 0
  immis_ni_count = 0

end subroutine ImmiscibleSetup

! ************************************************************************** !

subroutine ImmiscibleInitializeTimestep(realization)
  !
  ! Update data in module prior to time step
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Realization_Subsurface_class
  use Upwind_Direction_module

  implicit none

  class(realization_subsurface_type) :: realization

  call ImmiscibleUpdateFixedAccum(realization)
  immis_ni_count = 0

end subroutine ImmiscibleInitializeTimestep

! ************************************************************************** !

subroutine ImmiscibleUpdateSolution(realization)
  !
  ! Updates data in module after a successful time
  ! step
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  if (realization%option%compute_mass_balance_new) then
    call ImmiscibleUpdateMassBalance(realization)
  endif

  immis_ts_count = immis_ts_count + 1
  immis_ts_cut_count = 0
  immis_ni_count = 0

end subroutine ImmiscibleUpdateSolution

! ************************************************************************** !

subroutine ImmiscibleTimeCut(realization)
  !
  ! Resets arrays for time step cut
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Realization_Subsurface_class

  implicit none

  class(realization_subsurface_type) :: realization

  immis_ts_cut_count = immis_ts_cut_count + 1

  call ImmiscibleInitializeTimestep(realization)

end subroutine ImmiscibleTimeCut

! ************************************************************************** !

subroutine ImmiscibleComputeMassBalance(realization,mass_balance)
  !
  ! Initializes mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Field_module
  use Grid_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization
  PetscReal :: mass_balance(2)

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: iphase

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  immiscible_auxvars => patch%aux%Immiscible%auxvars
  material_auxvars => patch%aux%Material%auxvars

  mass_balance = 0.d0

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    do iphase = 1, option%nflowdof
      ! volume_phase = saturation*porosity*volume
      mass_balance(iphase) = mass_balance(iphase) + &
          immiscible_auxvars(ZERO_INTEGER,ghosted_id)%den_kg(iphase) * &
          immiscible_auxvars(ZERO_INTEGER,ghosted_id)%sat(iphase) * &
          material_auxvars(ghosted_id)%porosity * &
          material_auxvars(ghosted_id)%volume
    enddo
  enddo

end subroutine ImmiscibleComputeMassBalance

! ************************************************************************** !

subroutine ImmiscibleZeroMassBalanceDelta(realization)
  !
  ! Zeros mass balance delta array
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%Immiscible%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance_delta = 0.d0
  enddo
  do iconn = 1, patch%aux%Immiscible%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance_delta = 0.d0
  enddo

end subroutine ImmiscibleZeroMassBalanceDelta

! ************************************************************************** !

subroutine ImmiscibleUpdateMassBalance(realization)
  !
  ! Updates mass balance
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)

  PetscInt :: iconn

  option => realization%option
  patch => realization%patch

  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss

  do iconn = 1, patch%aux%Immiscible%num_aux_bc
    global_auxvars_bc(iconn)%mass_balance(1,:) = &
      global_auxvars_bc(iconn)%mass_balance(1,:) + &
      global_auxvars_bc(iconn)%mass_balance_delta(1,:)*option%flow_dt
  enddo
  do iconn = 1, patch%aux%Immiscible%num_aux_ss
    global_auxvars_ss(iconn)%mass_balance(1,:) = &
      global_auxvars_ss(iconn)%mass_balance(1,:) + &
      global_auxvars_ss(iconn)%mass_balance_delta(1,:)*option%flow_dt
  enddo

end subroutine ImmiscibleUpdateMassBalance

! ************************************************************************** !

subroutine ImmiscibleUpdateAuxVars(realization)
  !
  ! Updates the auxiliary variables associated with the Immiscible problem
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Coupler_module
  use Connection_module
  use Material_module
  use Material_Aux_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(connection_set_type), pointer :: cur_connection_set
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_bc(:)
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_ss(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)

  PetscInt :: ghosted_id, local_id, sum_connection, iconn, natural_id
  PetscInt :: ghosted_start, ghosted_end, ghosted_offset
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal :: xxbc(realization%option%nflowdof)
  PetscInt :: gas_pressure_index, gas_saturation_index
  PetscErrorCode :: ierr

  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field

  immiscible_auxvars => patch%aux%Immiscible%auxvars
  immiscible_auxvars_bc => patch%aux%Immiscible%auxvars_bc
  immiscible_auxvars_ss => patch%aux%Immiscible%auxvars_ss
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  call VecGetArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  do ghosted_id = 1, grid%ngmax
    if (grid%nG2L(ghosted_id) < 0) cycle ! bypass ghosted corner cells
    !geh - Ignore inactive cells with inactive materials
    if (patch%imat(ghosted_id) <= 0) cycle
    ! IMMIS_UPDATE_FOR_ACCUM indicates call from non-perturbation
    option%iflag = IMMIS_UPDATE_FOR_ACCUM
    natural_id = grid%nG2A(ghosted_id)
    ghosted_end = ghosted_id * option%nflowdof
    ghosted_start = ghosted_end - option%nflowdof + 1
    if (grid%nG2L(ghosted_id) == 0) natural_id = -natural_id
    call ImmiscibleAuxVarCompute(xx_loc_p(ghosted_start:ghosted_end), &
                            immiscible_auxvars(ZERO_INTEGER,ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            patch%characteristic_curves_array( &
                              patch%cc_id(ghosted_id))%ptr, &
                            natural_id, &
                            PETSC_TRUE,option)
  enddo

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    gas_pressure_index = &
      boundary_condition%flow_aux_mapping(IMMIS_GAS_PRESSURE_INDEX)
    gas_saturation_index = &
      boundary_condition%flow_aux_mapping(IMMIS_GAS_SATURATION_INDEX)
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      !geh: negate to indicate boundary connection, not actual cell
      natural_id = -grid%nG2A(ghosted_id)
      ghosted_offset = (ghosted_id-1)*option%nflowdof
      if (gas_pressure_index > 0 .and. gas_saturation_index > 0) then
        xxbc(IMMIS_GAS_PRESSURE_DOF) = &
            boundary_condition%flow_aux_real_var(gas_pressure_index,iconn)
        xxbc(IMMIS_GAS_SATURATION_DOF) = &
            boundary_condition%flow_aux_real_var(gas_saturation_index,iconn)
      else
        xxbc(:) = xx_loc_p(ghosted_offset+1:ghosted_offset+2)
      endif
      ! IMMIS_UPDATE_FOR_BOUNDARY indicates call from non-perturbation
      option%iflag = IMMIS_UPDATE_FOR_BOUNDARY
      call ImmiscibleAuxVarCompute(xxbc,immiscible_auxvars_bc(sum_connection), &
                              global_auxvars_bc(sum_connection), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              natural_id, &
                              PETSC_FALSE,option)
    enddo
    boundary_condition => boundary_condition%next
  enddo

  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit
    cur_connection_set => source_sink%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      call ImmiscibleAuxVarCopy(immiscible_auxvars(ZERO_INTEGER,ghosted_id), &
                                immiscible_auxvars_ss(sum_connection),option)
      call GlobalAuxVarCopy(global_auxvars(ghosted_id), &
                            global_auxvars_ss(sum_connection),option)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArray(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)

  patch%aux%Immiscible%auxvars_up_to_date = PETSC_TRUE

end subroutine ImmiscibleUpdateAuxVars

! ************************************************************************** !

subroutine ImmiscibleUpdateFixedAccum(realization)
  !
  ! Updates the fixed portion of the
  ! accumulation term
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Patch_module
  use Option_module
  use Field_module
  use Grid_module
  use Material_Aux_module
  use Petsc_Utility_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(field_type), pointer :: field
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(material_parameter_type), pointer :: material_parameter

  PetscInt :: ghosted_id, local_id, local_start, local_end, natural_id
  PetscInt :: imat
  PetscReal, pointer :: xx_p(:)
  PetscReal, pointer :: accum_t_p(:)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: ndof

  PetscErrorCode :: ierr

  if (.not.immis_calc_accum) return

  option => realization%option
  field => realization%field
  patch => realization%patch
  grid => patch%grid

  ndof = option%nflowdof

  immiscible_auxvars => patch%aux%Immiscible%auxvars
  global_auxvars => patch%aux%Global%auxvars
  material_auxvars => patch%aux%Material%auxvars
  material_parameter => patch%aux%Material%material_parameter

  call VecGetArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecGetArray(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

  do local_id = 1, grid%nlmax
    ghosted_id = grid%nL2G(local_id)
    !geh - Ignore inactive cells with inactive materials
    imat = patch%imat(ghosted_id)
    if (imat <= 0) cycle
    local_end = local_id * ndof
    local_start = local_end - ndof + 1
    natural_id = grid%nG2A(ghosted_id)
    ! IMMIS_UPDATE_FOR_FIXED_ACCUM indicates call from non-perturbation
    option%iflag = IMMIS_UPDATE_FOR_FIXED_ACCUM
    call ImmiscibleAuxVarCompute(xx_p(local_start:local_end), &
                            immiscible_auxvars(ZERO_INTEGER,ghosted_id), &
                            global_auxvars(ghosted_id), &
                            material_auxvars(ghosted_id), &
                            patch%characteristic_curves_array( &
                              patch%cc_id(ghosted_id))%ptr, &
                            natural_id, &
                            PETSC_TRUE,option)
    call ImmiscibleAccumulation(immiscible_auxvars(ZERO_INTEGER,ghosted_id), &
                           global_auxvars(ghosted_id), &
                           material_auxvars(ghosted_id), &
                           option,Res,Jdum,.not.immis_numerical_derivatives)
    call PetUtilVecSVBL(accum_t_p,local_id,Res,ndof,PETSC_TRUE)
  enddo

  call VecRestoreArrayRead(field%flow_xx,xx_p,ierr);CHKERRQ(ierr)
  call VecRestoreArray(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

end subroutine ImmiscibleUpdateFixedAccum

! ************************************************************************** !

subroutine ImmiscibleResidual(snes,xx,r,A,realization,debug,ierr)
  !
  ! Computes the residual equation
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Field_module
  use Patch_module
  use Discretization_module
  use Option_module
  use Debug_module

  use Connection_module
  use Coupler_module
  use Debug_module
  use Grid_module
  use Material_Aux_module
  use Upwind_Direction_module
  use Petsc_Utility_module
  use Matrix_Zeroing_module

  implicit none

  SNES :: snes
  Vec :: xx
  Vec :: r
  Mat :: A
  class(realization_subsurface_type) :: realization
  type(debug_type), pointer :: debug
  PetscErrorCode :: ierr

!  PetscViewer :: viewer
  type(discretization_type), pointer :: discretization
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch
  type(option_type), pointer :: option
  type(field_type), pointer :: field
  type(coupler_type), pointer :: boundary_condition
  type(coupler_type), pointer :: source_sink
  type(material_parameter_type), pointer :: material_parameter
  type(immiscible_parameter_type), pointer :: immiscible_parameter
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars(:,:)
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_ss(:)
  type(material_auxvar_type), pointer :: material_auxvars(:)
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set

  PetscInt :: iconn
  PetscReal :: scale
  PetscReal :: ss_flow_vol_flux(IMMIS_MAX_DOF)
  PetscInt :: sum_connection
  PetscInt :: local_id, ghosted_id
  PetscInt :: local_id_up, local_id_dn, ghosted_id_up, ghosted_id_dn
  PetscInt :: imat, imat_up, imat_dn

  PetscReal, pointer :: r_p(:)
  PetscReal, pointer :: accum_t_p(:), accum_tpdt_p(:)
  PetscReal, pointer :: xx_loc_p(:)
  PetscReal, pointer :: vec_p(:)

!  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: icc_up, icc_dn
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: v_darcy(IMMIS_MAX_DOF)

  PetscInt :: ndof
  PetscInt :: istart, iend

  ndof = realization%option%nflowdof

  discretization => realization%discretization
  option => realization%option
  patch => realization%patch
  grid => patch%grid
  field => realization%field
  material_parameter => patch%aux%Material%material_parameter
  immiscible_auxvars => patch%aux%Immiscible%auxvars
  immiscible_auxvars_bc => patch%aux%Immiscible%auxvars_bc
  immiscible_parameter => patch%aux%Immiscible%immiscible_parameter
  global_auxvars => patch%aux%Global%auxvars
  global_auxvars_bc => patch%aux%Global%auxvars_bc
  global_auxvars_ss => patch%aux%Global%auxvars_ss
  material_auxvars => patch%aux%Material%auxvars

  call MatZeroEntries(A,ierr);CHKERRQ(ierr)

  ! Communication -----------------------------------------
  ! must be called before ImmiscibleUpdateAuxVars()
  call DiscretizationGlobalToLocal(discretization,xx,field%flow_xx_loc,NFLOWDOF)
  call ImmiscibleUpdateAuxVars(realization)

  ! override flags since they will soon be out of date
  patch%aux%Immiscible%auxvars_up_to_date = PETSC_FALSE

  if (immis_numerical_derivatives) then
    ! Perturb aux vars
    call VecGetArrayRead(field%flow_xx_loc,xx_loc_p,ierr);CHKERRQ(ierr)
    do ghosted_id = 1, grid%ngmax  ! For each local node do...
      if (patch%imat(ghosted_id) <= 0) cycle
      iend = ghosted_id*ndof
      istart = iend-ndof+1
      call ImmiscibleAuxVarPerturb(xx_loc_p(istart:iend), &
                              immiscible_auxvars(:,ghosted_id), &
                              global_auxvars(ghosted_id), &
                              material_auxvars(ghosted_id), &
                              patch%characteristic_curves_array( &
                                patch%cc_id(ghosted_id))%ptr, &
                              grid%nG2A(ghosted_id),option)
    enddo
    call VecRestoreArrayRead(field%flow_xx_loc,xx_loc_p, &
                                ierr);CHKERRQ(ierr)
  endif

  if (option%compute_mass_balance_new) then
    call ImmiscibleZeroMassBalanceDelta(realization)
  endif

  option%iflag = 1
  ! now assign access pointer to local variables
  call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)

  ! Accumulation terms ------------------------------------
  ! accumulation at t(k) (doesn't change during Newton iteration)
  if (immis_calc_accum) then
    call VecGetArrayRead(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)
    r_p = -accum_t_p
    call VecRestoreArrayRead(field%flow_accum_t,accum_t_p,ierr);CHKERRQ(ierr)

    ! accumulation at t(k+1)
    call VecGetArray(field%flow_accum_tpdt,accum_tpdt_p,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      !geh - Ignore inactive cells with inactive materials
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      call ImmiscibleAccumDerivative(immiscible_auxvars(:,ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                option,Res,Jup)
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
      call PetUtilVecSVBL(accum_tpdt_p,local_id,Res,ndof,PETSC_TRUE)
      call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jup,ndof)
    enddo
    call VecRestoreArray(field%flow_accum_tpdt,accum_tpdt_p,ierr);CHKERRQ(ierr)
  else
    r_p = 0.d0
  endif

  if (immis_calc_flux) then
    ! Interior Flux Terms -----------------------------------
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0
    do
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        ghosted_id_up = cur_connection_set%id_up(iconn)
        ghosted_id_dn = cur_connection_set%id_dn(iconn)

        local_id_up = grid%nG2L(ghosted_id_up) ! = zero for ghost nodes
        local_id_dn = grid%nG2L(ghosted_id_dn) ! Ghost to local mapping

        imat_up = patch%imat(ghosted_id_up)
        imat_dn = patch%imat(ghosted_id_dn)
        if (imat_up <= 0 .or. imat_dn <= 0) cycle

        icc_up = patch%cc_id(ghosted_id_up)
        icc_dn = patch%cc_id(ghosted_id_dn)

        call XXFluxDerivative(immiscible_auxvars(:,ghosted_id_up), &
                              global_auxvars(ghosted_id_up), &
                              material_auxvars(ghosted_id_up), &
                              immiscible_auxvars(:,ghosted_id_dn), &
                              global_auxvars(ghosted_id_dn), &
                              material_auxvars(ghosted_id_dn), &
                              cur_connection_set%area(iconn), &
                              cur_connection_set%dist(:,iconn), &
                              immiscible_parameter,option,v_darcy, &
                              Res,Jup,Jdn, &
                              PUCast(local_id_up == immis_debug_cell_id .or. &
                                     local_id_dn == immis_debug_cell_id))
        patch%internal_velocities(:,sum_connection) = v_darcy(:)
        if (associated(patch%internal_flow_fluxes)) then
          patch%internal_flow_fluxes(:,sum_connection) = Res(:)
        endif

        if (local_id_up > 0) then
          call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
          call PetUtilMatSVBL(A,ghosted_id_up,ghosted_id_up,Jup,ndof)
          call PetUtilMatSVBL(A,ghosted_id_up,ghosted_id_dn,Jdn,ndof)
        endif

        if (local_id_dn > 0) then
          Res = -Res
          call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
          Jup = -Jup
          Jdn = -Jdn
          call PetUtilMatSVBL(A,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
          call PetUtilMatSVBL(A,ghosted_id_dn,ghosted_id_up,Jup,ndof)
        endif
      enddo

      cur_connection_set => cur_connection_set%next
    enddo
  endif

  if (immis_calc_bcflux) then
    ! Boundary Flux Terms -----------------------------------
    boundary_condition => patch%boundary_condition_list%first
    sum_connection = 0
    do
      if (.not.associated(boundary_condition)) exit

      cur_connection_set => boundary_condition%connection_set

      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1

        local_id = cur_connection_set%id_dn(iconn)
        ghosted_id = grid%nL2G(local_id)

        imat_dn = patch%imat(ghosted_id)
        if (imat_dn <= 0) cycle

        icc_dn = patch%cc_id(ghosted_id)

        call XXBCFluxDerivative(boundary_condition%flow_bc_type, &
                                boundary_condition%flow_aux_mapping, &
                                boundary_condition% &
                                  flow_aux_real_var(:,iconn), &
                                immiscible_auxvars_bc(sum_connection), &
                                global_auxvars_bc(sum_connection), &
                                immiscible_auxvars(:,ghosted_id), &
                                global_auxvars(ghosted_id), &
                                material_auxvars(ghosted_id), &
                                cur_connection_set%area(iconn), &
                                cur_connection_set%dist(:,iconn), &
                                immiscible_parameter,option, &
                                v_darcy,Res,Jdn, &
                                PUCast(local_id == immis_debug_cell_id))
        patch%boundary_velocities(:,sum_connection) = v_darcy(:)
        if (associated(patch%boundary_flow_fluxes)) then
          patch%boundary_flow_fluxes(:,sum_connection) = Res(:)
        endif
        if (option%compute_mass_balance_new) then
          ! contribution to boundary
          global_auxvars_bc(sum_connection)%mass_balance_delta(1,:) = &
            global_auxvars_bc(sum_connection)%mass_balance_delta(1,:) - Res(:)
        endif
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
        Jdn = -Jdn
        call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jdn,ndof)
      enddo
      boundary_condition => boundary_condition%next
    enddo
  endif

  ! Source/sink terms -------------------------------------
  source_sink => patch%source_sink_list%first
  sum_connection = 0
  do
    if (.not.associated(source_sink)) exit

    cur_connection_set => source_sink%connection_set

    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      ghosted_id = grid%nL2G(local_id)
      if (patch%imat(ghosted_id) <= 0) cycle
      if (associated(source_sink%flow_aux_real_var)) then
        scale = source_sink%flow_aux_real_var( &
                  source_sink%flow_aux_mapping(IMMIS_RATE_SCALE_INDEX),iconn)
      else
        scale = 1.d0
      endif
      call ImmiscibleSrcSinkDerivative(option,source_sink,scale, &
                                       immiscible_auxvars(:,ghosted_id), &
                                       global_auxvars(ghosted_id), &
                                       material_auxvars(ghosted_id), &
                                       ss_flow_vol_flux,Res,Jdn)
      if (associated(patch%ss_flow_vol_fluxes)) then
        patch%ss_flow_vol_fluxes(:,sum_connection) = ss_flow_vol_flux
      endif
      if (associated(patch%ss_flow_fluxes)) then
        patch%ss_flow_fluxes(:,sum_connection) = Res(:)
      endif
      if (option%compute_mass_balance_new) then
        ! contribution to boundary
        global_auxvars_ss(sum_connection)%mass_balance_delta(1,:) = &
          global_auxvars_ss(sum_connection)%mass_balance_delta(1,:) - Res(:)
      endif
      Res = -Res
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
      call PetUtilMatSVBL(A,ghosted_id,ghosted_id,Jdn,ndof)
    enddo
    source_sink => source_sink%next
  enddo

  call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)

  call MatrixZeroingZeroVecEntries(realization%patch%aux%Immiscible% &
                                     matrix_zeroing,r)

  if (immis_simultaneous_res_jac_calc) then

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      ! zero out inactive cells

    call MatrixZeroingZeroMatEntries(realization%patch%aux%Immiscible% &
                                       matrix_zeroing,A)
  endif

  ! Mass Transfer
  if (.not.PetscObjectIsNull(field%flow_mass_transfer)) then
    ! scale by -1.d0 for contribution to residual.  A negative contribution
    ! indicates mass being added to system.
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
    call VecGetArray(field%flow_mass_transfer,vec_p,ierr);CHKERRQ(ierr)
    ! geh: leave in expanded do loop form instead of VecAXPY for flexibility
    !      in the future
    do local_id = 1, grid%nlmax  ! For each local node do...
      ghosted_id = grid%nL2G(local_id)
      imat = patch%imat(ghosted_id)
      if (imat <= 0) cycle
      r_p(local_id) = r_p(local_id) - vec_p(local_id)
    enddo
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
    call VecRestoreArray(field%flow_mass_transfer,vec_p, &
                            ierr);CHKERRQ(ierr)
  endif

  if (debug%vecview_residual) then
    call DebugVecView(debug,r,'Immis_residual','', &
                      immis_ts_count,immis_ts_cut_count, &
                      immis_ni_count,option)
  endif
  if (debug%vecview_solution) then
    call DebugVecView(debug,xx,'Immis_xx','', &
                      immis_ts_count,immis_ts_cut_count, &
                      immis_ni_count,option)
  endif

end subroutine ImmiscibleResidual

! ************************************************************************** !

subroutine ImmiscibleSetPlotVariables(realization,list)
  !
  ! Adds variables to be printed to list
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Output_Aux_module
  use Variables_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(output_variable_list_type), pointer :: list

  character(len=MAXWORDLENGTH) :: name, units

  if (associated(list%first)) then
    return
  endif

  if (list%flow_vars) then

    name = 'Gas Pressure'
    units = 'Pa'
    call OutputVariableAddToList(list,name,OUTPUT_PRESSURE,units, &
                                 GAS_PRESSURE)

    name = 'Gas Saturation'
    units = ''
    call OutputVariableAddToList(list,name,OUTPUT_SATURATION,units, &
                                 Gas_SATURATION)

  endif

end subroutine ImmiscibleSetPlotVariables

! ************************************************************************** !

subroutine ImmiscibleMapBCAuxVarsToGlobal(realization)
  !
  ! Maps variables in immiscible auxvar to global equivalent.
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Option_module
  use Patch_module
  use Coupler_module
  use Connection_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(patch_type), pointer :: patch
  type(coupler_type), pointer :: boundary_condition
  type(connection_set_type), pointer :: cur_connection_set
  type(immiscible_auxvar_type), pointer :: immiscible_auxvars_bc(:)
  type(global_auxvar_type), pointer :: global_auxvars_bc(:)

  PetscInt :: sum_connection, iconn

  option => realization%option
  patch => realization%patch

  if (option%ntrandof == 0) return ! no need to update

  immiscible_auxvars_bc => patch%aux%Immiscible%auxvars_bc
  global_auxvars_bc => patch%aux%Global%auxvars_bc

  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0
  do
    if (.not.associated(boundary_condition)) exit
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      global_auxvars_bc(sum_connection)%sat = &
        immiscible_auxvars_bc(sum_connection)%sat
      global_auxvars_bc(sum_connection)%den_kg = &
        immiscible_auxvars_bc(sum_connection)%den_kg
    enddo
    boundary_condition => boundary_condition%next
  enddo

end subroutine ImmiscibleMapBCAuxVarsToGlobal

! ************************************************************************** !

subroutine ImmiscibleDestroy(realization)
  !
  ! Deallocates variables associated with Richard
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Realization_Subsurface_class
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization

  ! place anything that needs to be freed here.
  ! auxvars are deallocated in auxiliary..

end subroutine ImmiscibleDestroy

end module Immiscible_module
