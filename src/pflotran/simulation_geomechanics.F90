module Simulation_Geomechanics_class

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Simulation_Subsurface_class
  use Geomechanics_Regression_module
  use PMC_Base_class
  use PMC_Subsurface_class
  use PMC_Geomechanics_class
  use Realization_Subsurface_class
  use Geomechanics_Realization_class
  use PFLOTRAN_Constants_module
  use Waypoint_module
  use Simulation_Aux_module
  use Output_Aux_module
  use Utility_module, only : Equal

  implicit none

  private

  type, public, extends(simulation_subsurface_type) :: &
                simulation_geomechanics_type
  contains
    procedure, public :: InitializeRun => GeomechanicsSimulationInitializeRun
    procedure, public :: InputRecord => GeomechanicsSimInputRecord
    procedure, public :: ExecuteRun => GeomechanicsSimulationExecuteRun
    procedure, public :: FinalizeRun => GeomechanicsSimulationFinalizeRun
    procedure, public :: Strip => GeomechanicsSimulationStrip
  end type simulation_geomechanics_type

  public :: GeomechanicsSimulationCreate, &
            GeomechanicsSimulationDestroy

contains

! ************************************************************************** !

function GeomechanicsSimulationCreate(driver,option)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !
  use Driver_class
  use Option_module

  implicit none

  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  class(simulation_geomechanics_type), pointer :: GeomechanicsSimulationCreate

#ifdef GEOMECH_DEBUG
  print *,'GeomechanicsSimulationCreate'
#endif

  allocate(GeomechanicsSimulationCreate)
  call GeomechanicsSimulationInit(GeomechanicsSimulationCreate,driver,option)

end function GeomechanicsSimulationCreate

! ************************************************************************** !

subroutine GeomechanicsSimulationInit(this,driver,option)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified: Satish Karra, 06/01/2016
  !
  use Waypoint_module
  use Driver_class
  use Option_module
  use Geomechanics_Attr_module

  implicit none

  class(simulation_geomechanics_type) :: this
  class(driver_type), pointer :: driver
  type(option_type), pointer :: option

  call SimSubsurfInit(this,driver,option)
  this%geomech => GeomechAttrCreate()
  this%geomech%waypoint_list => WaypointListCreate()

end subroutine GeomechanicsSimulationInit

! ************************************************************************** !

subroutine GeomechanicsSimulationInitializeRun(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  use Output_module
  use PMC_Geomechanics_class

  implicit none

  class(simulation_geomechanics_type) :: this

#ifdef GEOMECH_DEBUG
  call PrintMsg(this%option,'Simulation%InitializeRun()')
#endif

  call SimSubsurfInitializeRun(this)

end subroutine GeomechanicsSimulationInitializeRun

! ************************************************************************** !

subroutine GeomechanicsSimInputRecord(this)
  !
  ! Writes ingested information to the input record file.
  !
  ! Author: Jenn Frederick, SNL
  ! Date: 03/17/2016
  !
  use Output_module

  implicit none

  class(simulation_geomechanics_type) :: this

  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') 'simulation type: '
  write(id,'(a)') 'geomechanics'

  ! print output file information
  call OutputInputRecord(this%output_option,this%geomech%waypoint_list)

end subroutine GeomechanicsSimInputRecord

! ************************************************************************** !

subroutine GeomechanicsSimulationExecuteRun(this)
  !
  ! Executes the geomechanics simulation with checkpoint/restart support.
  !
  ! For the decoupled path (no geomech realization), delegates to
  ! SimSubsurfExecuteRun which handles checkpoints via waypoint_list_outer.
  !
  ! For the coupled path, integrates checkpoint handling into the
  ! dt_coupling loop: after each RunToTime call, any checkpoint waypoints
  ! that have been reached trigger a checkpoint. A final checkpoint is
  ! written on successful completion.
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified: J.A. Angeles, 2025 - added checkpoint/restart support
  !

  use Waypoint_module
  use Timestepper_Base_class
  use Checkpoint_module

  implicit none

  class(simulation_geomechanics_type) :: this

  PetscReal :: time
  PetscReal :: final_time
  PetscReal :: dt
  PetscReal :: target_time
  PetscReal, parameter :: tolerance = 1.d-3
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXSTRINGLENGTH) :: append_name

  time = this%option%time
  final_time = SimSubsurfGetFinalWaypointTime(this)

#ifdef GEOMECH_DEBUG
  call PrintMsg(this%option,'GeomechanicsSimulationExecuteRun()')
#endif

  if (.not.associated(this%geomech%realization)) then
    ! Decoupled path: delegate to parent which handles checkpoints
    call SimSubsurfExecuteRun(this)

  else

    dt = this%geomech%realization%dt_coupling

    ! If simulation is decoupled subsurface-geomech simulation, set
    ! dt_coupling to be dt_max
    if (Equal(this%geomech%realization%dt_coupling,0.d0)) then
      this%option%io_buffer = 'Set non-zero COUPLING_TIME_SIZE in GEOMECHANICS_TIME.'
      call PrintErrMsg(this%option)
    else

      ! Clear checkpoint flags from subsurface waypoints. In the coupled
      ! geomechanics path, checkpoints are handled explicitly in this
      ! loop via waypoint_list_outer. If we leave the flags on the
      ! subsurface waypoints, the flow PMC (which has is_master=TRUE)
      ! will also trigger its own checkpoint inside PMCBaseRunToTime,
      ! resulting in duplicate checkpoint writes.
      cur_waypoint => this%waypoint_list_subsurface%first
      do while (associated(cur_waypoint))
        cur_waypoint%print_checkpoint = PETSC_FALSE
        cur_waypoint => cur_waypoint%next
      enddo

      ! Initialize waypoint pointer for checkpoint tracking
      cur_waypoint => this%waypoint_list_outer%first
      ! Checkpoint at initial time if the first waypoint says so (restart case)
      if (associated(cur_waypoint)) then
        if (cur_waypoint%print_checkpoint .and. &
            Equal(cur_waypoint%time,this%option%time)) then
          append_name = &
            CheckpointAppendNameAtTime( &
              this%process_model_coupler_list%option%time, &
              this%process_model_coupler_list%option)
          call this%process_model_coupler_list%Checkpoint(append_name)
        endif
      endif
      call WaypointSkipToTime(cur_waypoint,this%option%time)

      do
        dt = this%geomech%realization%dt_coupling

        ! Compute next target time from dt_coupling
        if (time + dt*(1.d0+tolerance) >= final_time) then
          target_time = final_time
        else
          target_time = time + dt
        endif

        ! If a checkpoint waypoint falls before our coupling target,
        ! run to the waypoint time first so checkpoint is written at
        ! the correct time.
        do while (associated(cur_waypoint))
          if (cur_waypoint%time > target_time) exit
          ! Waypoint is at or before our target — run to it
          this%geomech%process_model_coupler%timestepper%dt = &
            cur_waypoint%time - time
          if (cur_waypoint%time - time > 0.d0) then
            call this%RunToTime(cur_waypoint%time)
            time = cur_waypoint%time
          endif
          if (this%stop_flag /= TS_CONTINUE) exit
          ! Write checkpoint if this waypoint requests it
          if (cur_waypoint%print_checkpoint) then
            append_name = &
              CheckpointAppendNameAtTime( &
                this%process_model_coupler_list%option%time, &
                this%process_model_coupler_list%option)
            call this%process_model_coupler_list%Checkpoint(append_name)
          endif
          cur_waypoint => cur_waypoint%next
        enddo

        if (this%stop_flag /= TS_CONTINUE) exit

        ! Run remaining interval to the coupling target
        if (target_time > time) then
          this%geomech%process_model_coupler%timestepper%dt = &
            target_time - time
          call this%RunToTime(target_time)
          time = target_time
        endif

        if (this%stop_flag /= TS_CONTINUE) exit ! end simulation

        if (time >= final_time) exit
      enddo

      ! Final checkpoint: only checkpoint successful simulations
      if (this%stop_flag /= TS_STOP_FAILURE) then
        select case(this%stop_flag)
          case(TS_STOP_MAX_TIME_STEP)
            append_name = '-restart-max-ts'
          case(TS_STOP_WALLCLOCK_EXCEEDED)
            append_name = '-restart-max-wc'
          case default ! TS_STOP_END_SIMULATION
            append_name = '-restart'
        end select
        if (associated(this%option%checkpoint)) then
          call this%process_model_coupler_list%Checkpoint(append_name)
        endif
      endif

    endif
  endif

end subroutine GeomechanicsSimulationExecuteRun


! ************************************************************************** !

subroutine GeomechanicsSimulationFinalizeRun(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified by Satish Karra, 06/22/16

  use Timestepper_Base_class

  implicit none

  class(simulation_geomechanics_type) :: this

#ifdef GEOMECH_DEBUG
  call PrintMsg(this%option,'GeomechanicsSimulationFinalizeRun')
#endif

  select case(this%stop_flag)
    case(TS_STOP_END_SIMULATION,TS_STOP_MAX_TIME_STEP)
      call SimSubsurfWriteRegression(this)
  end select

  call SimSubsurfFinalizeRun(this)

end subroutine GeomechanicsSimulationFinalizeRun

! ************************************************************************** !

subroutine GeomechanicsSimulationStrip(this)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  ! Modified by Satish Karra, 06/01/16
  !

  implicit none

  class(simulation_geomechanics_type) :: this

#ifdef GEOMECH_DEBUG
  call PrintMsg(this%option,'GeomechanicsSimulationStrip()')
#endif

  call GeomechanicsRegressionDestroy(this%geomech%regression)
  call WaypointListDestroy(this%geomech%waypoint_list)
  call SimSubsurfStrip(this)
  call WaypointListDestroy(this%waypoint_list_subsurface)

end subroutine GeomechanicsSimulationStrip

! ************************************************************************** !

subroutine GeomechanicsSimulationDestroy(simulation)
  !
  ! This routine
  !
  ! Author: Gautam Bisht, LBNL
  ! Date: 01/01/14
  !

  implicit none

  class(simulation_geomechanics_type), pointer :: simulation

  if (.not.associated(simulation)) return

  call simulation%Strip()
  deallocate(simulation)
  nullify(simulation)

end subroutine GeomechanicsSimulationDestroy

end module Simulation_Geomechanics_class
