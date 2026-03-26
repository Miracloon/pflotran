module TH_Well_module

#include "petsc/finclude/petscmat.h"
  use petscmat

  use PFLOTRAN_Constants_module
  use TH_Aux_module
  use PM_Well_class

  implicit none

  private

  PetscReal :: liquid_pressure = 1.d6
  PetscReal :: thermal_diffusivity = 0.143d-6

  PetscInt, parameter :: INLET_BC = 1
  PetscInt, parameter :: OUTLET_BC = 2

  public :: THWellSetup, &
            THWellUpdateAuxVars, &
            THWellPerturb, &
            THWell, &
            THWellAccumulationTerms, &
            THWellMatrixZeroing, &
            THWellUpdateHeatFlux

contains

! ************************************************************************** !

subroutine THWellSetup(pm_well_base,realization)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Connection_module
  use Realization_Subsurface_class
  use Patch_module
  use Grid_module
  use Matrix_Zeroing_module
  use String_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base
  class(realization_subsurface_type), pointer :: realization

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(connection_set_type), pointer :: connection_set
  type(connection_set_type), pointer :: bc_connection_set
  type(th_parameter_type), pointer :: th_parameter
  type(th_auxvar_type), pointer :: auxvars(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscReal :: dist_up, dist_dn
  PetscReal :: surface_area
  PetscInt :: local_id
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, j, k, iconn, idof
  PetscInt :: nx, ny, nz
  PetscInt :: offset_i, offset_j, offset_k
  PetscInt :: istart, iend, kstart, kend, dk
  PetscInt :: num_connections
  PetscInt :: num_well_cells
  PetscReal :: pipe_cross_sectional_area
  PetscReal :: rock_thermal_conductivity
  PetscInt, allocatable :: inactive_cells(:)
  PetscReal, allocatable :: len_(:,:)
  PetscReal :: segment_length
  PetscReal :: dx, dy, dz
  PetscReal :: tempreal

  if (.not.associated(pm_well_base)) return

  th_scale_by_volume = PETSC_FALSE
  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  select case(grid%itype)
    case(STRUCTURED_GRID)
    case default
      if (pm_well%iscenario /= 4) then
        option%io_buffer = 'THWell currently supports only structured grids.'
        call PrintErrMsg(option)
      endif
  end select

  pm_well%realization_base => realization
  call pm_well%Setup()

  th_parameter => patch%aux%TH%th_parameter
  auxvars => patch%aux%TH%auxvars

  if (Uninitialized(pm_well%inlet_temperature)) then
    option%io_buffer = 'Well INLET_TEMPERATURE must be specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%pipe_inner_diameter)) then
    option%io_buffer = 'Well PIPE_INNER_DIAMETER must be specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%pipe_wall_thickness) .and. &
      Initialized(pm_well%pipe_wall_thermal_conductivity)) then
    option%io_buffer = 'Well PIPE_WALL_THICKNESS must be specified when &
      &PIPE_WALL_THERMAL_CONDUCTIVITY is specified.'
    call PrintErrMsg(option)
  endif
  if (Initialized(pm_well%pipe_wall_thickness) .and. &
      Uninitialized(pm_well%pipe_wall_thermal_conductivity)) then
    option%io_buffer = 'Well PIPE_WALL_THERMAL_CONDUCTIVITY must be &
      &specified when PIPE_WALL_THICKNESS is specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%insulator_thickness) .and. &
      Initialized(pm_well%insulator_thermal_conductivity)) then
    option%io_buffer = 'Well INSULATOR_THICKNESS must be specified when &
      &INSULATOR_THERMAL_CONDUCTIVITY is specified.'
    call PrintErrMsg(option)
  endif
  if (Initialized(pm_well%insulator_thickness) .and. &
      Uninitialized(pm_well%insulator_thermal_conductivity)) then
    option%io_buffer = 'Well INSULATOR_THERMAL_CONDUCTIVITY must be &
      &specified when INSULATOR_THICKNESS is specified.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(pm_well%flow_velocity)) then
    option%io_buffer = 'Well PIPE_FLOW_VELOCITY must be specified.'
    call PrintErrMsg(option)
  elseif (pm_well%flow_velocity < 0.d0) then
    option%io_buffer = 'Well PIPE_FLOW_VELOCITY must be non-negative.'
    call PrintErrMsg(option)
  endif

  pipe_cross_sectional_area = PI*(0.5d0*pm_well%pipe_inner_diameter)**2

  nx = grid%structured_grid%nx
  ny = grid%structured_grid%ny
  nz = grid%structured_grid%nz
  select case(pm_well%iscenario)
    case(1)
      num_well_cells = nx
    case(2)
      num_well_cells = 380
    case(3)
      num_well_cells = 2*(int(0.75001d0*grid%structured_grid%nz)) + &
                       int(0.5001d0*grid%structured_grid%nx) - 2
    case(4)
      num_well_cells = pm_well%well_grid%nsegments
    case default
      option%io_buffer = 'Unrecognized scenario in THWellSetup()'
      call PrintErrMsg(option)
  end select
  num_connections = num_well_cells-1

  connection_set => &
    ConnectionCreate(num_connections,INTERNAL_FACE_CONNECTION_TYPE,NULL_GRID)

  bc_connection_set => &
    ConnectionCreate(2,BOUNDARY_FACE_CONNECTION_TYPE,NULL_GRID)

  allocate(len_(3,grid%ngmax))
  len_ = 0.d0
  select case(pm_well%iscenario)
    case(1)
      offset_i = 0
      offset_j = int(ny/2)
      offset_k = int(nz/2)
      do iconn = 1, num_connections
        ghosted_id_up = nx*ny*offset_k + nx*offset_j + offset_i + iconn
        ghosted_id_dn = ghosted_id_up + 1
        connection_set%id_up(iconn) = ghosted_id_up
        connection_set%id_dn(iconn) = ghosted_id_dn
        connection_set%dist(:,iconn) = 0.d0
        dist_up = 0.5d0*grid%structured_grid%dx(ghosted_id_up)
        dist_dn = 0.5d0*grid%structured_grid%dx(ghosted_id_dn)
        connection_set%dist(0,iconn) = dist_up + dist_dn
        connection_set%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connection_set%dist(1,iconn) = 1.d0
        connection_set%area(iconn) = pipe_cross_sectional_area
        len_(1,ghosted_id_up) = len_(1,ghosted_id_up) + dist_up
        len_(1,ghosted_id_dn) = len_(1,ghosted_id_dn) + dist_dn
        if (iconn == 1) then
          bc_connection_set%id_dn(INLET_BC) = ghosted_id_up
          bc_connection_set%dist(:,INLET_BC) = 0.d0
          bc_connection_set%dist(0,INLET_BC) = dist_up
          bc_connection_set%dist(-1,INLET_BC) = 0.d0 ! fraction upwind
          bc_connection_set%dist(1,INLET_BC) = 1.d0
          bc_connection_set%area(INLET_BC) = pipe_cross_sectional_area
          len_(1,ghosted_id_up) = len_(1,ghosted_id_up) + dist_up
        endif
        if (iconn == num_connections) then
          bc_connection_set%id_dn(OUTLET_BC) = ghosted_id_dn
          bc_connection_set%dist(:,OUTLET_BC) = 0.d0
          bc_connection_set%dist(0,OUTLET_BC) = dist_dn
          bc_connection_set%dist(-1,OUTLET_BC) = 1.d0
          bc_connection_set%dist(1,OUTLET_BC) = -1.d0
          bc_connection_set%area(OUTLET_BC) = pipe_cross_sectional_area
          len_(1,ghosted_id_dn) = len_(1,ghosted_id_dn) + dist_dn
        endif
      enddo
    case(2,3)
      select case(pm_well%iscenario)
        case(2)
          i = 22
          j = 7
          kstart = 102
          kend = 33
          dk = -1
        case(3)
          i = int(0.25001d0*(grid%structured_grid%nx))+1
          j = int(0.5d0*(grid%structured_grid%ny))+1
          kstart = grid%structured_grid%nz
          kend = int(0.25001d0*(grid%structured_grid%nz))+2
          dk = -1
      end select
      iconn = 0
      do k = kstart, kend, dk
        iconn = iconn + 1
        ghosted_id_up = nx*ny*(k-1) + nx*(j-1) + i
        ghosted_id_dn = ghosted_id_up - nx*ny
        connection_set%id_up(iconn) = ghosted_id_up
        connection_set%id_dn(iconn) = ghosted_id_dn
        connection_set%dist(:,iconn) = 0.d0
        dist_up = 0.5d0*grid%structured_grid%dz(ghosted_id_up)
        dist_dn = 0.5d0*grid%structured_grid%dz(ghosted_id_dn)
        connection_set%dist(0,iconn) = dist_up + dist_dn
        connection_set%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connection_set%dist(3,iconn) = -1.d0
        connection_set%area(iconn) = pipe_cross_sectional_area
        len_(3,ghosted_id_up) = len_(3,ghosted_id_up) + dist_up
        len_(3,ghosted_id_dn) = len_(3,ghosted_id_dn) + dist_dn
        if (k == kstart) then
          bc_connection_set%id_dn(INLET_BC) = ghosted_id_up
          bc_connection_set%dist(:,INLET_BC) = 0.d0
          bc_connection_set%dist(0,INLET_BC) = dist_up
          bc_connection_set%dist(-1,INLET_BC) = 0.d0 ! fraction upwind
          bc_connection_set%dist(3,INLET_BC) = -1.d0
          bc_connection_set%area(INLET_BC) = pipe_cross_sectional_area
          len_(3,ghosted_id_up) = len_(3,ghosted_id_up) + dist_up
        endif
      enddo
      select case(pm_well%iscenario)
        case(2)
          j = 7
          k = 32
          istart = 22
          iend = 260
        case(3)
          k = int(0.25001d0*(grid%structured_grid%nz))+1
          j = int(0.5d0*(grid%structured_grid%ny))+1
          istart = int(0.25001d0*(grid%structured_grid%nx))+1
          iend = int(0.75001d0*(grid%structured_grid%nx))-1
      end select
      do i = istart, iend
        iconn = iconn + 1
        ghosted_id_up = nx*ny*(k-1) + nx*(j-1) + i
        ghosted_id_dn = ghosted_id_up + 1
        connection_set%id_up(iconn) = ghosted_id_up
        connection_set%id_dn(iconn) = ghosted_id_dn
        connection_set%dist(:,iconn) = 0.d0
        dist_up = 0.5d0*grid%structured_grid%dx(i)
        dist_dn = 0.5d0*grid%structured_grid%dx(i+1)
        connection_set%dist(0,iconn) = dist_up + dist_dn
        connection_set%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connection_set%dist(1,iconn) = 1.d0
        connection_set%area(iconn) = pipe_cross_sectional_area
        len_(1,ghosted_id_up) = len_(1,ghosted_id_up) + dist_up
        len_(1,ghosted_id_dn) = len_(1,ghosted_id_dn) + dist_dn
      enddo
      select case(pm_well%iscenario)
        case(2)
          i = 261
          j = 7
          kstart = 32
          kend = 101
        case(3)
          i = int(0.75001d0*(grid%structured_grid%nx))
          j = int(0.5d0*(grid%structured_grid%ny))+1
          kstart = int(0.25001d0*(grid%structured_grid%nz))+1
          kend = grid%structured_grid%nz-1
      end select
      do k = kstart, kend
        iconn = iconn + 1
        ghosted_id_up = nx*ny*(k-1) + nx*(j-1) + i
        ghosted_id_dn = ghosted_id_up + nx*ny
        connection_set%id_up(iconn) = ghosted_id_up
        connection_set%id_dn(iconn) = ghosted_id_dn
        connection_set%dist(:,iconn) = 0.d0
        dist_up = 0.5d0*grid%structured_grid%dz(k)
        dist_dn = 0.5d0*grid%structured_grid%dz(k+1)
        connection_set%dist(0,iconn) = dist_up + dist_dn
        connection_set%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
        connection_set%dist(3,iconn) = 1.d0
        connection_set%area(iconn) = pipe_cross_sectional_area
        len_(3,ghosted_id_up) = len_(3,ghosted_id_up) + dist_up
        len_(3,ghosted_id_dn) = len_(3,ghosted_id_dn) + dist_dn
        if (k == kend) then
          bc_connection_set%id_dn(OUTLET_BC) = ghosted_id_dn
          bc_connection_set%dist(:,OUTLET_BC) = 0.d0
          bc_connection_set%dist(0,OUTLET_BC) = dist_dn
          bc_connection_set%dist(-1,OUTLET_BC) = 1.d0 ! fraction upwind
          bc_connection_set%dist(3,OUTLET_BC) = -1.d0
          bc_connection_set%area(OUTLET_BC) = pipe_cross_sectional_area
          len_(3,ghosted_id_dn) = len_(3,ghosted_id_dn) + dist_dn
        endif
      enddo
      if (iconn /= num_connections) then
        option%io_buffer = 'Mismatch on number of connections in THWellsetup: ' // &
          StringWrite(iconn) // ' vs ' // StringWrite(num_connections)
        call PrintErrMsg(option)
      endif
    case(4)
      iconn = 0
      ! num_connections set to pm_well%well_grid%nsegments above
      do i = 1, num_connections
        iconn = iconn + 1
        ghosted_id_up = pm_well%well_grid%h_ghosted_id(i)
        ghosted_id_dn = pm_well%well_grid%h_ghosted_id(i+1)
        connection_set%id_up(iconn) = ghosted_id_up
        connection_set%id_dn(iconn) = ghosted_id_dn
        dx = pm_well%well_grid%h(i+1)%x - pm_well%well_grid%h(i)%x
        dy = pm_well%well_grid%h(i+1)%y - pm_well%well_grid%h(i)%y
        dz = pm_well%well_grid%h(i+1)%z - pm_well%well_grid%h(i)%z
        tempreal = sqrt(dx*dx+dy*dy+dz*dz)
        connection_set%dist(-1,iconn) = UNINITIALIZED_DOUBLE
        connection_set%dist(0,iconn) = tempreal
        connection_set%dist(1,iconn) = dx/tempreal
        connection_set%dist(2,iconn) = dy/tempreal
        connection_set%dist(3,iconn) = dz/tempreal
        connection_set%area(iconn) = pipe_cross_sectional_area
        if (i == 1) then
          dx = pm_well%well_grid%h(i)%x - pm_well%well_grid%tophole(1)
          dy = pm_well%well_grid%h(i)%y - pm_well%well_grid%tophole(2)
          dz = pm_well%well_grid%h(i)%z - pm_well%well_grid%tophole(3)
          tempreal = sqrt(dx*dx+dy*dy+dz*dz)
          bc_connection_set%id_dn(INLET_BC) = ghosted_id_up
          bc_connection_set%dist(-1,INLET_BC) = UNINITIALIZED_DOUBLE
          bc_connection_set%dist(0,INLET_BC) = tempreal
          bc_connection_set%dist(1:3,INLET_BC) = [dx,dy,dz]/tempreal
          bc_connection_set%area(INLET_BC) = pipe_cross_sectional_area
        endif
        if (i == num_connections) then
          dx = pm_well%well_grid%bottomhole(1) - pm_well%well_grid%h(i+1)%x
          dy = pm_well%well_grid%bottomhole(2) - pm_well%well_grid%h(i+1)%y
          dz = pm_well%well_grid%bottomhole(3) - pm_well%well_grid%h(i+1)%z
          tempreal = sqrt(dx*dx+dy*dy+dz*dz)
          bc_connection_set%id_dn(OUTLET_BC) = ghosted_id_dn
          bc_connection_set%dist(-1,OUTLET_BC) = UNINITIALIZED_DOUBLE
          bc_connection_set%dist(0,OUTLET_BC) = tempreal
          bc_connection_set%dist(1:3,OUTLET_BC) = [dx,dy,dz]/tempreal
          bc_connection_set%area(OUTLET_BC) = pipe_cross_sectional_area
        endif
      enddo
      do i = 1, num_well_cells
        ghosted_id = pm_well%well_grid%h_ghosted_id(i)
        len_(:,ghosted_id) = [pm_well%well_grid%dx(i), &
                              pm_well%well_grid%dy(i), &
                              pm_well%well_grid%dz(i)]
      enddo
  end select
  pm_well%connection_set => connection_set
  pm_well%bc_connection_set => bc_connection_set

  allocate(pm_well%well_cells(num_well_cells))
  if (num_well_cells > 1) then
    do iconn = 1, num_connections
      local_id = grid%nG2L(connection_set%id_up(iconn))
      pm_well%well_cells(iconn) = local_id
      if (iconn == num_connections) then
        local_id = grid%nG2L(connection_set%id_dn(iconn))
        pm_well%well_cells(iconn+1) = local_id
      endif
    enddo
  else
    pm_well%well_cells(1) = 1
    len_(1,1) = grid%structured_grid%dx(1)
  endif

  ! allocate auxvar data structures for all grid cells
  allocate(auxvars_well(num_well_cells))
  do i = 1, num_well_cells
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    auxvars(ghosted_id)%iwellaux = i
    call THWellAuxVarInit(auxvars_well(i),option, &
                          th_numerical_derivatives)
    auxvars_well(i)%local_id = pm_well%well_cells(i)
    segment_length = sqrt(len_(1,ghosted_id)**2 + &
                          len_(2,ghosted_id)**2 + &
                          len_(3,ghosted_id)**2)
    auxvars_well(i)%segment_length = segment_length
    auxvars_well(i)%volume = segment_length * &
                             pipe_cross_sectional_area
    surface_area = segment_length * pm_well%pipe_inner_diameter * PI
    rock_thermal_conductivity = th_parameter%ckwet(patch%cct_id(ghosted_id))
    call THWellWI(rock_thermal_conductivity, &
                  rock_thermal_conductivity, &
                  rock_thermal_conductivity, &
                  grid%structured_grid%dx(ghosted_id), &
                  grid%structured_grid%dy(ghosted_id), &
                  grid%structured_grid%dz(ghosted_id), &
                  len_(1,ghosted_id), &
                  len_(2,ghosted_id), &
                  len_(3,ghosted_id), &
                  pm_well%pipe_inner_diameter, &
                  0.d0, & ! s
                  auxvars_well(i)%well_index)
    auxvars_well(i)%therm_cond_borehole_to_cell = auxvars_well(i)%well_index / &
                                                  surface_area
    if (th_numerical_derivatives) then
      ! must copy parameters down to perturbed auxvars
      call THWellAuxVarCopyParamsToPert(auxvars_well(i))
    endif
  enddo
  deallocate(len_)

  patch%aux%TH%auxvars_well => auxvars_well
  patch%aux%TH%num_well_aux = num_well_cells
  allocate(auxvars_well(2))
  do i = 1, 2
    call THWellAuxVarInit(auxvars_well(i),option,th_numerical_derivatives)
  enddo
  patch%aux%TH%auxvars_well_bc => auxvars_well
  patch%aux%TH%num_well_aux_bc = 2

  ! generate matrix zeroing
  pm_well%matrix_zeroing => MatrixZeroingCreate()
  allocate(inactive_cells(grid%nlmax))
  inactive_cells = 1
  do i = 1, size(pm_well%well_cells)
    inactive_cells(pm_well%well_cells(i)) = 0
  enddo
  i = 0
  do local_id = 1, grid%nlmax
    if (inactive_cells(local_id) == 1) then
      i = i + 1
      inactive_cells(i) = local_id
    endif
  enddo
  if (i < grid%nlmax) then
    inactive_cells(i+1:) = UNINITIALIZED_INTEGER
  endif
  call MatrixZeroingAllocateArray(pm_well%matrix_zeroing,i,option)
  do i = 1, pm_well%matrix_zeroing%n_zero_rows
    local_id = inactive_cells(i)
    ghosted_id = grid%nL2G(local_id)
    idof = (local_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local(i) = idof ! one-based
    idof = (ghosted_id-1)*option%nflowdof + th_well_dof
    pm_well%matrix_zeroing%zero_rows_local_ghosted(i) = idof-1 ! zero-based
  enddo
  deallocate(inactive_cells)

end subroutine THWellSetup

! ************************************************************************** !

subroutine THWellUpdateAuxVars(xx,pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  PetscReal :: xx(:) ! ghosted solution vec
  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: i, ibc
  PetscInt :: istart, iend
  PetscInt :: num_well_cells
  PetscReal :: xxbc(3)

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  num_well_cells = size(pm_well%well_cells)
  do i = 1, num_well_cells
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    iend = ghosted_id * option%nflowdof
    istart = iend - option%nflowdof + 1
    call THWellAuxVarCompute(liquid_pressure,xx(istart:iend), &
                             auxvars_well(i), &
                             natural_id,option)
  enddo

  xxbc = UNINITIALIZED_DOUBLE
  do ibc = 1, 2
    select case(ibc)
      case(INLET_BC)
        xxbc(th_well_dof) = pm_well%inlet_temperature
      case(OUTLET_BC)
        ! to guarantee zero conduction at boundary
        xxbc(th_well_dof) = auxvars_well(size(pm_well%well_cells))%temp
    end select
    call THWellAuxVarCompute(liquid_pressure,xxbc, &
                             auxvars_well_bc(ibc), &
                             UNINITIALIZED_INTEGER,option)
  enddo

end subroutine THWellUpdateAuxVars

! ************************************************************************** !

subroutine THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  use Option_module
  use Petsc_Utility_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  Vec :: r
  Mat :: J
  class(pm_well_type), pointer :: pm_well_base
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscReal :: Res(3)
  PetscReal :: Jblock(3,3)
  PetscReal :: dt
  PetscInt :: ndof
  PetscInt :: local_id
  PetscInt :: ghosted_id
  PetscInt :: i
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars_well => patch%aux%TH%auxvars_well
  dt = option%flow_dt

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

  do i = 1, size(pm_well%well_cells)
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    call THWellAccumulationDerivative(auxvars_well(i),Res,Jblock, &
                                      ndof,calculate_derivatives)
    Res = Res/dt
    Jblock = Jblock/dt
    if (calculate_derivatives) then
      call PetUtilMatSVBL(J,ghosted_id,ghosted_id,Jblock,ndof)
    else
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
    endif
  enddo

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWellAccumulationTerms

! ************************************************************************** !

subroutine THWell(r,J,pm_well_base,calculate_derivatives)
  !
  ! Calculates contribution of well model to the TH residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25
  !
  use Option_module
  use Petsc_Utility_module
  use Connection_module
  use Patch_module
  use Grid_module

  implicit none

  Vec :: r
  Mat :: J
  class(pm_well_type), pointer :: pm_well_base
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  PetscReal, pointer :: r_p(:)
  type(th_auxvar_type), pointer :: auxvars(:)
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  type(th_well_auxvar_type), pointer :: auxvars_well_bc(:)
  type(connection_set_type), pointer :: connection_set
  type(th_well_auxvar_type), pointer :: auxvar_well_up, auxvar_well_dn
  PetscReal :: Res(3)
  PetscReal :: Jup(3,3), Jdn(3,3)
  PetscInt :: ndof
  PetscInt :: local_id, local_id_up, local_id_dn
  PetscInt :: ghosted_id, ghosted_id_up, ghosted_id_dn
  PetscInt :: i, iconn, ibc
  PetscReal :: velocity_scale
  PetscErrorCode :: ierr

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  ndof = option%nflowdof
  auxvars => patch%aux%TH%auxvars
  auxvars_well => patch%aux%TH%auxvars_well
  auxvars_well_bc => patch%aux%TH%auxvars_well_bc

  call THWellAccumulationTerms(r,J,pm_well_base,calculate_derivatives)

  if (.not.calculate_derivatives) then
    call VecGetArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

#if 1
  ! calculate heat transfer between well cells
  connection_set => pm_well%connection_set
  i = 0
  do iconn = 1, connection_set%num_connections
    ! cells are ghosted in connection set
    ghosted_id_up = connection_set%id_up(iconn)
    ghosted_id_dn = connection_set%id_dn(iconn)
    local_id_up = grid%nG2L(ghosted_id_up)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    i = i + 1
    call THWellFluxDerivative(auxvars_well(i), &
                    auxvars_well(i+1), &
                    pm_well%flow_velocity, &
                    connection_set%dist(0,iconn), &
                    connection_set%area(iconn), &
                    Res,Jup,Jdn,ndof,calculate_derivatives, &
                    PETSC_FALSE)
    if (calculate_derivatives) then
      if (local_id_up > 0) then
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_up,Jup,ndof)
        call PetUtilMatSVBL(J,ghosted_id_up,ghosted_id_dn,Jdn,ndof)
      endif
      if (local_id_dn > 0) then
        Jup = -Jup
        Jdn = -Jdn
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
        call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_up,Jup,ndof)
      endif
    else
      if (local_id_up > 0) then
        call PetUtilVecSVBL(r_p,local_id_up,Res,ndof,PETSC_FALSE)
      endif
      if (local_id_dn > 0) then
        Res = -Res
        call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
      endif
    endif

  enddo
#endif

#if 1
  do ibc = 1, 2
    auxvar_well_up => auxvars_well_bc(ibc)
    ghosted_id_dn = pm_well%bc_connection_set%id_dn(ibc)
    local_id_dn = grid%nG2L(ghosted_id_dn)
    select case(ibc)
      case(INLET_BC)
        auxvar_well_dn => auxvars_well(1)
        velocity_scale = 1.d0
      case(OUTLET_BC) ! outflow
        ! note that the temperature is currently equal to the last well cell,
        ! so this is effectively a zero conduction bc
        auxvar_well_dn => auxvars_well(size(pm_well%well_cells))
        velocity_scale = -1.d0
    end select
    call THWellFluxDerivative(auxvar_well_up, &
                    auxvar_well_dn, &
                    pm_well%flow_velocity*velocity_scale, &
                    pm_well%bc_connection_set%dist(0,ibc), &
                    pm_well%bc_connection_set%area(ibc), &
                    Res,Jup,Jdn,ndof,calculate_derivatives, &
                    PETSC_TRUE)
    if (calculate_derivatives) then
      Jdn = -Jdn
      call PetUtilMatSVBL(J,ghosted_id_dn,ghosted_id_dn,Jdn,ndof)
    else
      Res = -Res
      pm_well%heat_flux(1,ibc) = -Res(3)
      call PetUtilVecSVBL(r_p,local_id_dn,Res,ndof,PETSC_FALSE)
    endif
  enddo
#endif

#if 1
  do i = 1, size(pm_well%well_cells)
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    call THWellHeatExchangeDerivative(auxvars(ghosted_id), &
                                      auxvars_well(i), &
                                      pm_well, &
                                      Res,Jdn, &
                                      ndof,calculate_derivatives)
    if (calculate_derivatives) then
      call PetUtilMatSVBL(J,ghosted_id,ghosted_id,Jdn,ndof)
    else
      call PetUtilVecSVBL(r_p,local_id,Res,ndof,PETSC_FALSE)
    endif
  enddo
#endif

  ! calculate heat transfer between well cells and subsurface cells
  if (.not.calculate_derivatives) then
    call VecRestoreArray(r,r_p,ierr);CHKERRQ(ierr)
  endif

end subroutine THWell

! ************************************************************************** !

subroutine THWellMatrixZeroing(pm_well_base,r,J,calculate_derivatives)
  !
  ! Zeros dof for non-well grid cells
  !
  ! Author: Glenn Hammond
  ! Date: 09/17/25
  !
  use Matrix_Zeroing_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base
  Vec :: r
  Mat :: J
  PetscBool :: calculate_derivatives

  class(pm_well_closed_loop_type), pointer :: pm_well

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  if (calculate_derivatives) then
    call MatrixZeroingZeroMatEntries(pm_well%matrix_zeroing,J)
  else
    call MatrixZeroingZeroVecEntries(pm_well%matrix_zeroing,r)
  endif

end subroutine THWellMatrixZeroing

! ************************************************************************** !

subroutine THWellPerturb(pm_well_base)
  !
  ! Configure the well
  !
  ! Author: Glenn Hammond
  ! Date: 09/12/25

  use Option_module
  use Patch_module
  use Grid_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(th_well_auxvar_type), pointer :: auxvars_well(:)
  PetscInt :: local_id, ghosted_id, natural_id
  PetscInt :: i
  PetscInt :: num_well_cells

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  patch => pm_well%realization%patch
  grid => patch%grid
  option => pm_well%option

  auxvars_well => patch%aux%TH%auxvars_well

  num_well_cells = size(pm_well%well_cells)
  do i = 1, num_well_cells
    local_id = pm_well%well_cells(i)
    ghosted_id = grid%nL2G(local_id)
    natural_id = grid%nG2A(ghosted_id)
    call THWellAuxVarPerturb(liquid_pressure,auxvars_well(i), &
                             natural_id,option)
  enddo

end subroutine THWellPerturb

! ************************************************************************** !

subroutine THWellAccumulation(auxvar_well,Res,Jac,ndof)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)

  PetscReal :: tempreal

  ! Res [MJ] = liquid_density_kmol [kmol liquid/m^3 liquid] *
  !            internal energy [MJ/kmol liquid] *
  !            volume [m^3 liquid]
  Res = 0.d0
  Jac = 0.d0
  tempreal = auxvar_well%volume * auxvar_well%den
  Res(ndof) = tempreal * auxvar_well%u
  Jac(ndof,ndof) = tempreal * auxvar_well%du_dT

end subroutine THWellAccumulation

! ************************************************************************** !

subroutine THWellAccumulationDerivative(auxvar_well,Res,Jac,ndof, &
                                        calculate_derivatives)
  !
  ! Calculates accumulation term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jac(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellAccumulation(auxvar_well,Res,Jac,ndof)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    do idof = 1, ndof
      call THWellAccumulation(auxvar_well%auxvar_pert(idof), &
                              Res_pert,Jdum,ndof)
      do ieq = 1, ndof
        Jac(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellAccumulationDerivative

! ************************************************************************** !

subroutine THWellFlux(auxvar_well_up,auxvar_well_dn, &
                      velocity,dist,area,Res,Jup,Jdn,ndof, &
                      calculate_derivatives)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: velocity
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: calculate_derivatives

  PetscReal :: tup, tdn
  PetscReal :: dh_dTup, dh_dTdn
  PetscReal :: enthalpy, q_kmol
  PetscReal :: diffusive_energy_flux, convective_energy_flux
  PetscReal :: tempreal

  Res = 0.d0

  tup = auxvar_well_up%temp
  tdn = auxvar_well_dn%temp

  ! convection
  dh_dTup = 0.d0
  dh_dTdn = 0.d0
  if (velocity > 0) then
    enthalpy = auxvar_well_up%h
    dh_dTup = auxvar_well_up%dh_dT
    q_kmol = velocity * area * auxvar_well_up%den
  else
    enthalpy = auxvar_well_dn%h
    dh_dTdn = auxvar_well_dn%dh_dT
    q_kmol = velocity * area * auxvar_well_dn%den
  endif

  convective_energy_flux = q_kmol * enthalpy

  ! conduction
  ! [MW] = [MW/m-K] * [m^2 bulk] [C] / [m bulk]
  tempreal = thermal_diffusivity * area / dist
  diffusive_energy_flux = tempreal * (tup-tdn)

  Res(th_well_eq) = convective_energy_flux + diffusive_energy_flux

  if (calculate_derivatives) then
    Jup = 0.d0
    Jup(th_well_eq,th_well_dof) = q_kmol * dh_dTup + tempreal
    Jdn = 0.d0
    Jdn(th_well_eq,th_well_dof) = q_kmol * dh_dTdn - tempreal
  endif

end subroutine THWellFlux

! ************************************************************************** !

subroutine THWellFluxDerivative(auxvar_well_up,auxvar_well_dn, &
                                velocity,dist,area,Res,Jup,Jdn,ndof, &
                                calculate_derivatives,is_bc)
  !
  ! Calculates internal flux term for TH well residual and Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 09/15/25
  !
  type(th_well_auxvar_type) :: auxvar_well_up, auxvar_well_dn
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: Jup(ndof,ndof)
  PetscReal :: Jdn(ndof,ndof)
  PetscReal :: velocity
  PetscReal :: dist
  PetscReal :: area
  PetscBool :: calculate_derivatives
  PetscBool :: is_bc

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof), Jdum1(ndof,ndof)
  PetscInt :: idof, ieq

  call THWellFlux(auxvar_well_up,auxvar_well_dn, &
                  velocity,dist,area,Res,Jup,Jdn,ndof, &
                  .not.th_numerical_derivatives)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    Jup = 0.d0
    if (.not.is_bc) then
      do idof = 1, ndof
        call THWellFlux(auxvar_well_up%auxvar_pert(idof),auxvar_well_dn, &
                        velocity,dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
        do ieq = 1, ndof
          Jup(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                          auxvar_well_up%auxvar_pert(idof)%pert
        enddo
      enddo
    endif
    Jdn = 0.d0
    do idof = 1, ndof
      call THWellFlux(auxvar_well_up,auxvar_well_dn%auxvar_pert(idof), &
                      velocity,dist,area,Res_pert,Jdum,Jdum1,ndof,PETSC_FALSE)
      do ieq = 1, ndof
        Jdn(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / &
                        auxvar_well_dn%auxvar_pert(idof)%pert
      enddo
    enddo
  endif

end subroutine THWellFluxDerivative

! ************************************************************************** !

function THWellResistanceWellIndex(auxvar_well,pm_well)
  !
  ! Calculates heat flux between pipe wall and insulator
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceWellIndex

  if (pm_well%use_well_index) then
    ! [K/MW] = 1 / [MW/K]
    THWellResistanceWellIndex = 1.d0 / auxvar_well%well_index
  else
    THWellResistanceWellIndex = 0.d0
  endif

end function THWellResistanceWellIndex

! ************************************************************************** !

function THWellResistanceInsulator(auxvar_well,pm_well)
  !
  ! Calculates heat flux between pipe wall and insulator
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceInsulator

  PetscReal :: r_inner
  PetscReal :: r_outer

  r_inner = 0.5d0 * pm_well%pipe_inner_diameter + pm_well%pipe_wall_thickness
  r_outer = r_inner + pm_well%insulator_thickness

  if (pm_well%insulator_thermal_conductivity > 0.d0 .and. &
       r_inner > 0.d0 .and. auxvar_well%segment_length > 0.d0) then
    ! Cylindrical conduction resistance: R = ln(r_o/r_i)/(2*pi*k*L) [K/MW]
    THWellResistanceInsulator = &
      log(r_outer/r_inner) / (2.d0 * PI * &
      pm_well%insulator_thermal_conductivity * auxvar_well%segment_length)
  elseif (pm_well%insulator_thermal_conductivity < 0.d0) then
    THWellResistanceInsulator = 0.d0
  else
    THWellResistanceInsulator = MAX_DOUBLE
  endif

end function THWellResistanceInsulator

! ************************************************************************** !

function THWellResistancePipeWall(auxvar_well,pm_well)
  !
  ! Calculates heat flux between fluid centerline and pipe wall
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistancePipeWall

  PetscReal :: r_inner
  PetscReal :: r_outer

  r_inner = 0.5d0 * pm_well%pipe_inner_diameter
  r_outer = r_inner + pm_well%pipe_wall_thickness

  if (pm_well%pipe_wall_thermal_conductivity > 0.d0 .and. &
      r_inner > 0.d0 .and. auxvar_well%segment_length > 0.d0) then
    ! Cylindrical conduction resistance: R = ln(r_o/r_i)/(2*pi*k*L) [K/MW]
    THWellResistancePipeWall = &
      log(r_outer/r_inner) / (2.d0 * PI * &
      pm_well%pipe_wall_thermal_conductivity * auxvar_well%segment_length)
  elseif (pm_well%pipe_wall_thermal_conductivity < 0.d0) then
    THWellResistancePipeWall = 0.d0
  else
    THWellResistancePipeWall = MAX_DOUBLE
  endif

end function THWellResistancePipeWall

! ************************************************************************** !

function THWellResistanceFluid(auxvar_well,pm_well)
  !
  ! Calculates heat flux between fluid centerline and pipe wall
  !
  ! Author: Glenn Hammond
  ! Date: 02/10/26
  !
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: THWellResistanceFluid
  PetscReal :: thermal_coeff_fluid
  PetscReal :: surface_area

  thermal_coeff_fluid = THWellComputeConvCoefficient(auxvar_well,pm_well)
  ! Inner wall surface area for convection
  surface_area = PI * pm_well%pipe_inner_diameter * auxvar_well%segment_length

  ! Convective resistance: R = 1/(h*A) [K/MW] as h [MW/m^2-K] and A [m^2]
  if (thermal_coeff_fluid > 0.d0 .and. surface_area > 0.d0) then
    THWellResistanceFluid = 1.d0 / (thermal_coeff_fluid * surface_area)
  elseif (thermal_coeff_fluid < 0.d0) then
    THWellResistanceFluid = 0.d0
  else
    THWellResistanceFluid = MAX_DOUBLE
  endif

end function THWellResistanceFluid

! ************************************************************************** !

subroutine THWellHeatExchange(auxvar,auxvar_well,pm_well, &
                              Res,J,ndof, &
                              calculate_derivatives)
  !
  ! Calculates heat flux between grid cell and pipe
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/25
  !
  type(th_auxvar_type) :: auxvar
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: J(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: tcell, tpipe
  PetscReal :: diffusive_energy_flux
  PetscReal :: total_thermal_resistance
  PetscReal :: tempreal

  Res = 0.d0

  tcell = auxvar%temp
  tpipe = auxvar_well%temp

  total_thermal_resistance = &  ! [K/MW]
    THWellResistanceFluid(auxvar_well,pm_well) + &
    THWellResistancePipeWall(auxvar_well,pm_well) + &
    THWellResistanceInsulator(auxvar_well,pm_well) + &
    THWellResistanceWellIndex(auxvar_well,pm_well)

  diffusive_energy_flux = (tcell-tpipe) / total_thermal_resistance

  Res(TH_ENERGY_EQUATION_INDEX) = diffusive_energy_flux
  Res(th_well_eq) = -diffusive_energy_flux

  if (calculate_derivatives) then
    J = 0.d0
    tempreal = 1.d0 / total_thermal_resistance
    J(TH_ENERGY_EQUATION_INDEX,th_well_dof) = -tempreal
    J(TH_ENERGY_EQUATION_INDEX,TH_TEMPERATURE_DOF) = tempreal
    J(th_well_eq,th_well_dof) = tempreal
    J(th_well_eq,TH_TEMPERATURE_DOF) = -tempreal
  endif

end subroutine THWellHeatExchange

! ************************************************************************** !
subroutine THWellHeatExchangeDerivative(auxvar,auxvar_well,pm_well, &
                                        Res,J,ndof, &
                                        calculate_derivatives)
  !
  ! Calculates heat flux between grid cell and pipe
  !
  ! Author: Glenn Hammond
  ! Date: 09/18/25
  !
  type(th_auxvar_type) :: auxvar
  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well
  PetscInt :: ndof
  PetscReal :: Res(ndof)
  PetscReal :: J(ndof,ndof)
  PetscBool :: calculate_derivatives

  PetscReal :: Res_pert(ndof)
  PetscReal :: Jdum(ndof,ndof)
  PetscReal :: pert
  PetscInt :: idof, ieq

  call THWellHeatExchange(auxvar,auxvar_well,pm_well, &
                          Res,J,ndof, &
                          .not.th_numerical_derivatives)
  if (calculate_derivatives .and. th_numerical_derivatives) then
    J = 0.d0
    do idof = 1, ndof
      call THWellHeatExchange(auxvar%auxvar_pert(idof), &
                              auxvar_well%auxvar_pert(idof),pm_well, &
                              Res_pert,Jdum,ndof,PETSC_FALSE)
      if (idof == th_well_dof) then
        pert = auxvar_well%auxvar_pert(idof)%pert
      else
        pert = auxvar%auxvar_pert(idof)%pert
      endif
      do ieq = 1, ndof
        J(ieq,idof) = (Res_pert(ieq) - Res(ieq)) / pert
      enddo
    enddo
  endif

end subroutine THWellHeatExchangeDerivative

! ************************************************************************** !

subroutine THWellWI(kx,ky,kz,dx,dy,dz,lenx,leny,lenz,pipe_diameter,s,wi)
  !
  ! Calculates well index for a well segment
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/25
  !
  implicit none

  PetscReal :: kx   ! conductivity in x
  PetscReal :: ky
  PetscReal :: kz
  PetscReal :: dx   ! grid spacing in x
  PetscReal :: dy
  PetscReal :: dz
  PetscReal :: lenx ! length of well segement in x
  PetscReal :: leny
  PetscReal :: lenz
  PetscReal :: pipe_diameter
  PetscReal :: s
  PetscReal :: wi

  PetscReal :: sqrt_kx_over_ky
  PetscReal :: sqrt_ky_over_kx
  PetscReal :: sqrt_ky_over_kz
  PetscReal :: sqrt_kz_over_ky
  PetscReal :: sqrt_kx_over_kz
  PetscReal :: sqrt_kz_over_kx
  PetscReal :: ro_x
  PetscReal :: ro_y
  PetscReal :: ro_z
  PetscReal :: wi_x
  PetscReal :: wi_y
  PetscReal :: wi_z
  PetscReal :: rw

  rw = 0.5d0 * pipe_diameter

  sqrt_kx_over_ky = sqrt(kx/ky)
  sqrt_ky_over_kx = sqrt(ky/kx)
  sqrt_ky_over_kz = sqrt(ky/kz)
  sqrt_kz_over_ky = sqrt(kz/ky)
  sqrt_kx_over_kz = sqrt(kx/kz)
  sqrt_kz_over_kx = sqrt(kz/kx)

  ro_x = 0.28d0 * sqrt(sqrt_ky_over_kz*dz**2 + sqrt_kz_over_ky*dy**2) / &
                  (sqrt(sqrt_ky_over_kz) + sqrt(sqrt_kz_over_ky))
  ro_y = 0.28d0 * sqrt(sqrt_kz_over_kx*dx**2 + sqrt_kx_over_kz*dz**2) / &
                  (sqrt(sqrt_kz_over_kx) + sqrt(sqrt_kx_over_kz))
  ro_z = 0.28d0 * sqrt(sqrt_ky_over_kx*dx**2 + sqrt_kx_over_ky*dy**2) / &
                  (sqrt(sqrt_ky_over_kx) + sqrt(sqrt_kx_over_ky))
!print *, 'kx: ', kx, ky, kz
!print *, 'dx: ', dx, dy, dz
!print *, 'ro_z: ', ro_z

  wi_x = 2.d0*PI*sqrt(ky*kz)*lenx/(log(ro_x/rw)+s)
  wi_y = 2.d0*PI*sqrt(kx*kz)*leny/(log(ro_y/rw)+s)
  wi_z = 2.d0*PI*sqrt(kx*ky)*lenz/(log(ro_z/rw)+s)

!print *, 'wi_z: ', wi_z, lenz, rw

  wi = sqrt(wi_x**2 + wi_y**2 + wi_z**2)

end subroutine THWellWI

! ************************************************************************** !

subroutine THWellComputePipeCoefficient(auxvar_well,pm_well)
  !
  ! Computes pipe wall conduction coefficient for cylindrical wall
  ! h_pipe = k_wall / (r_outer * ln(r_outer/r_inner))
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: thermal_cond_pipe_wall ! wall thermal conductivity [MW/m-K]
  PetscReal :: r_inner                ! inner radius [m]
  PetscReal :: r_outer                ! outer radius [m]
  PetscReal :: h_pipe_wall            ! [MW/m^2-K]

  thermal_cond_pipe_wall = pm_well%pipe_wall_thermal_conductivity
  r_inner = 0.5d0 * pm_well%pipe_inner_diameter
  r_outer = r_inner + pm_well%pipe_wall_thickness

  if (r_outer > r_inner .and. r_inner > 0.d0) then
    ! For cylindrical wall: R_wall = ln(r_o/r_i)/(2*pi*L*k)
    ! Per unit area (dividing by 2*pi*r_o*L): h = k/(r_o*ln(r_o/r_i))
    h_pipe_wall = thermal_cond_pipe_wall / (r_outer * log(r_outer/r_inner))
  else if (r_outer <= r_inner) then
    ! No wall thickness - infinite conductivity (no resistance)
    h_pipe_wall = 1.d20
  else
    h_pipe_wall = 0.d0
  endif

end subroutine THWellComputePipeCoefficient

! ************************************************************************** !

function THWellComputeConvCoefficient(auxvar_well,pm_well)
  !
  ! Computes convection heat transfer coefficient from fluid flow in well
  ! h = Nu * k / D
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  type(th_well_auxvar_type) :: auxvar_well
  class(pm_well_closed_loop_type) :: pm_well

  PetscReal :: density
  PetscReal :: velocity
  PetscReal :: pipe_diameter
  PetscReal :: viscosity
  PetscReal :: heat_capacity
  PetscReal :: thermal_cond
  PetscReal :: THWellComputeConvCoefficient

  PetscReal :: reynolds
  PetscReal :: prandtl
  PetscReal :: nusselt

  reynolds = 0.d0
  prandtl = 0.d0
  nusselt = 0.d0

  velocity = pm_well%flow_velocity
  pipe_diameter = pm_well%pipe_inner_diameter

  density = auxvar_well%den_kg
  viscosity = auxvar_well%vis
  heat_capacity = auxvar_well%spec_heat_fluid
  thermal_cond = auxvar_well%therm_cond_fluid

  !Compute Reynolds number
  reynolds = THWellComputeReynoldNumber(density,velocity,pipe_diameter, &
                                        viscosity)
  !Compute Prandtl number
  prandtl = THWellComputePrandtlNumber(heat_capacity,viscosity,thermal_cond)
  !Compute Nusselt number
  nusselt = THWellComputeNusseltNumber(reynolds,prandtl)

  if (pipe_diameter > 0.d0) then
    THWellComputeConvCoefficient = nusselt * thermal_cond / pipe_diameter
  else
    THWellComputeConvCoefficient = 0.d0
  endif

end function THWellComputeConvCoefficient

! ************************************************************************** !

function THWellComputeReynoldNumber(density,velocity,pipe_diameter, &
                                    viscosity)
  !
  ! Computes Reynolds number for borehole flow
  ! Re = rho * u * D / mu
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/16/26
  !
  implicit none

  PetscReal :: density
  PetscReal :: velocity
  PetscReal :: pipe_diameter
  PetscReal :: viscosity
  PetscReal :: THWellComputeReynoldNumber

  if (viscosity > 0.d0) then
    THWellComputeReynoldNumber = density * dabs(velocity) * pipe_diameter / &
                                 viscosity
  else
    THWellComputeReynoldNumber = 0.d0
  endif

end function THWellComputeReynoldNumber

! ************************************************************************** !

function THWellComputePrandtlNumber(heat_capacity,viscosity,thermal_cond)
  !
  ! Computes Prandtl number for fluid
  ! Pr = c_p * mu / k
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/17/26
  !
  implicit none

  PetscReal :: heat_capacity
  PetscReal :: viscosity
  PetscReal :: thermal_cond
  PetscReal :: THWellComputePrandtlNumber

  if (thermal_cond > 0.d0) then
    THWellComputePrandtlNumber = heat_capacity * viscosity / thermal_cond
  else
    THWellComputePrandtlNumber = 0.d0
  endif

end function THWellComputePrandtlNumber

! ************************************************************************** !

function THWellComputeNusseltNumber(reynolds,prandtl)
  !
  ! Computes Nusselt number using appropriate correlation
  !
  ! Author: Piyoosh Jaysaval
  ! Date: 02/18/26
  !
  implicit none

  PetscReal :: reynolds
  PetscReal :: prandtl
  PetscReal :: THWellComputeNusseltNumber

  !TODO: maybe add OGS approach OR Gnielinski correlation?
  if (reynolds > 2300.d0) then
    ! turbulent flow,
    ! Zhang et al. 2015. "The analytical solution of the water-rock heat
    ! transfer coefficient and sensitivity analyses of parameters."
    ! Proceedings World Geothermal Congress 2015, Melbourne, Australia,
    ! 19-25, April 2015, pg 6
    THWellComputeNusseltNumber = 8.7d-5 * (reynolds**0.92d0) * &
                                 (prandtl**1.89d0)
  elseif (reynolds > 0.d0) then
    ! Laminar flow - use Nu = 4.36 for constant heat flux boundary condition
    ! (Nu = 3.66 for constant wall temperature)
    THWellComputeNusseltNumber = 4.36d0
  else
    THWellComputeNusseltNumber = 0.d0
  endif

end function THWellComputeNusseltNumber

! ************************************************************************** !

subroutine THWellUpdateHeatFlux(pm_well_base)
  !
  ! Update the cumulative heat flux for the well model
  !
  ! Author: Glenn Hammond
  ! Date: 09/24/25
  !
  use Option_module

  implicit none

  class(pm_well_type), pointer :: pm_well_base

  class(pm_well_closed_loop_type), pointer :: pm_well
!  PetscInt :: i, ii

  if (.not.associated(pm_well_base)) return

  pm_well => PMWellCastToClosedLoop(pm_well_base)

  pm_well%heat_flux(2,:) = pm_well%heat_flux(2,:) + &
                           pm_well%heat_flux(1,:) * pm_well%option%flow_dt

#if 0
  print *, 'TIME: ', pm_well%option%time, pm_well%option%flow_dt
  print *, 'INSTANTANEOUS HEAT_FLUX: ', pm_well%heat_flux(1,:)
  print *, 'CUMULATIVE HEAT_FLUX: ', pm_well%heat_flux(2,:)
  i = size(pm_well%realization%patch%aux%TH%auxvars_well)
  print *, 'PIPE TEMPERATURES: ', &
    pm_well%realization%patch%aux%TH%auxvars_well(1)%temp, &
    pm_well%realization%patch%aux%TH%auxvars_well(i)%temp
  i = UNINITIALIZED_INTEGER
  II = UNINITIALIZED_INTEGER
  select case(pm_well%iscenario)
    case(1)
      i = 1
      ii = 2
    case(3)
      i = 3506
      ii = 3515
  end select
  if (i > 0 .and. ii > 0) then
  print *, 'CELL TEMPERATURES: ', &
    pm_well%realization%patch%aux%TH%auxvars(i)%temp, &
    pm_well%realization%patch%aux%TH%auxvars(ii)%temp
  endif
#endif

end subroutine THWellUpdateHeatFlux

end module TH_Well_module

