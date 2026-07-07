module Richards_Aux_module

#include "petsc/finclude/petscsys.h"
  use petscsys

#ifdef BUFFER_MATRIX
  use Matrix_Buffer_module
#endif

  use PFLOTRAN_Constants_module
  use Matrix_Zeroing_module
  use Richards_Aniso_module
  use Geomechanics_Linear_Aux_module

  implicit none

  private

  PetscReal, public :: richards_itol_scaled_res = 1.d-5
  PetscReal, public :: richards_itol_rel_update = UNINITIALIZED_DOUBLE
  PetscReal, public :: richards_density_kmol_to_kg = FMWH2O
  PetscInt, public :: richards_ni_count
  PetscInt, public :: richards_ts_cut_count
  PetscInt, public :: richards_ts_count

  PetscInt, parameter, public :: RICHARDS_PRESSURE_DOF = 1
  PetscInt, parameter, public :: RICHARDS_CONDUCTANCE_DOF = 2
  PetscInt, public :: RICHARDS_WELL_DOF = UNINITIALIZED_INTEGER

  PetscInt, parameter, public :: RICHARDS_UPDATE_FOR_FIXED_ACCUM = 0
  PetscInt, parameter, public :: RICHARDS_UPDATE_FOR_ACCUM = 1
  PetscInt, parameter, public :: RICHARDS_UPDATE_FOR_BOUNDARY = 2
  PetscInt, parameter, public :: RICHARDS_UPDATE_FOR_PERTURBATION = 3

  ! Well
  PetscInt, public :: richards_well_coupling = UNINITIALIZED_INTEGER
  PetscInt, parameter, public :: RICHARDS_FULLY_IMPLICIT_WELL = ONE_INTEGER

  type, public :: richards_well_aux_type
    PetscReal :: pl   ! liquid pressure
    PetscReal :: sl   ! liquid saturation
    PetscReal :: dpl  ! reservoir-well liquid pressure differential
    PetscReal :: Ql   ! liquid exchange flux [kg/s]
    PetscReal :: bh_p ! bottom hole pressure
    PetscReal :: pressure_bump ! pressure change for initialization  ! TODO: JOE: Can I remove this as I am not creating a numerical jacobian?
  end type richards_well_aux_type

  type, public :: richards_auxvar_type

    PetscReal :: pc
    PetscReal :: kvr
    PetscReal :: kr
    PetscReal :: visc
    PetscReal :: dkvr_dp
    PetscReal :: dsat_dp
    PetscReal :: dden_dp
#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    PetscReal :: bc_alpha  ! Brooks Corey parameterization: alpha
    PetscReal :: bc_lambda ! Brooks Corey parameterization: lambda
#endif

    PetscReal :: d2sat_dp2
    PetscReal :: d2den_dp2
    PetscReal :: mass
    PetscReal :: dpres_dtime
    PetscReal :: dmass_dtime

    PetscReal :: effective_porosity
    PetscReal :: dpor_dp

    type(richards_well_aux_type), pointer :: well
  end type richards_auxvar_type

  type, public :: richards_parameter_type
    type(geomech_linear_parameter_type), pointer :: geomech_parameter
  end type richards_parameter_type

  type, public :: richards_type
    PetscBool :: auxvars_up_to_date
    PetscBool :: auxvars_cell_pressures_up_to_date
    PetscInt :: num_aux, num_aux_bc, num_aux_ss
    type(richards_auxvar_type), pointer :: auxvars(:)
    type(richards_auxvar_type), pointer :: auxvars_bc(:)
    type(richards_auxvar_type), pointer :: auxvars_ss(:)
    type(aniso_richards_data_type), pointer :: aniso_richards_data(:) ! per connection, when needed
#ifdef BUFFER_MATRIX
    type(matrix_buffer_type), pointer :: matrix_buffer
#endif
    type(matrix_zeroing_type), pointer :: matrix_zeroing
    type(richards_parameter_type), pointer :: richards_parameter
  end type richards_type

  PetscReal, parameter :: perturbation_tolerance = 1.d-6

  public :: RichardsAuxCreate, RichardsAuxDestroy, &
            RichardsAuxVarCompute, RichardsAuxVarInit, &
            RichardsAuxVarCopy, &
            RichardsAuxVarCompute2ndOrderDeriv

contains

! ************************************************************************** !

function RichardsAuxCreate()
  !
  ! Allocate and initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  use Option_module

  implicit none

  type(richards_type), pointer :: RichardsAuxCreate

  type(richards_type), pointer :: aux

  allocate(aux)
  aux%auxvars_up_to_date = PETSC_FALSE
  aux%auxvars_cell_pressures_up_to_date = PETSC_FALSE
  aux%num_aux = 0
  aux%num_aux_bc = 0
  aux%num_aux_ss = 0
  nullify(aux%auxvars)
  nullify(aux%auxvars_bc)
  nullify(aux%auxvars_ss)
  nullify(aux%aniso_richards_data)
#ifdef BUFFER_MATRIX
  nullify(aux%matrix_buffer)
#endif
  nullify(aux%matrix_zeroing)
  allocate(aux%richards_parameter)
  nullify(aux%richards_parameter%geomech_parameter)

  RichardsAuxCreate => aux

end function RichardsAuxCreate

! ************************************************************************** !

subroutine RichardsAuxVarInit(auxvar,option)
  !
  ! Initialize auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  use Option_module

  implicit none

  type(richards_auxvar_type) :: auxvar
  type(option_type) :: option

  auxvar%pc = 0.d0

  auxvar%kvr = 0.d0
  auxvar%visc = 0.d0
  auxvar%kr = 0.d0
  auxvar%dkvr_dp = 0.d0

  auxvar%dsat_dp = 0.d0
  auxvar%dden_dp = 0.d0

  auxvar%d2sat_dp2 = 0.d0
  auxvar%d2den_dp2 = 0.d0
  auxvar%mass = 0.d0
  auxvar%dpres_dtime = 0.d0
  auxvar%dmass_dtime = 0.0d0

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar%bc_alpha  = 0.0d0
  auxvar%bc_lambda  = 0.0d0
#endif

  auxvar%effective_porosity = 0.0d0
  auxvar%dpor_dp = 0.0d0

  if (richards_well_coupling > ZERO_INTEGER) then
    allocate(auxvar%well)
    auxvar%well%pl = UNINITIALIZED_DOUBLE
    auxvar%well%sl = UNINITIALIZED_DOUBLE
    auxvar%well%dpl = UNINITIALIZED_DOUBLE
    auxvar%well%Ql = UNINITIALIZED_DOUBLE
    auxvar%well%bh_p = UNINITIALIZED_DOUBLE
    auxvar%well%pressure_bump = 0.d0
  else
    nullify(auxvar%well)
  endif

end subroutine RichardsAuxVarInit

! ************************************************************************** !

subroutine RichardsAuxVarCopy(auxvar,auxvar2,option)
  !
  ! Copies an auxiliary variable
  !
  ! Author: Glenn Hammond
  ! Date: 12/13/07
  !

  use Option_module

  implicit none

  type(richards_auxvar_type) :: auxvar, auxvar2
  type(option_type) :: option

  auxvar2%pc = auxvar%pc

  auxvar2%kvr = auxvar%kvr
  auxvar2%visc = auxvar%visc
  auxvar2%kr = auxvar%kr
  auxvar2%dkvr_dp = auxvar%dkvr_dp

  auxvar2%dsat_dp = auxvar%dsat_dp
  auxvar2%dden_dp = auxvar%dden_dp

  auxvar2%d2sat_dp2 = auxvar%d2sat_dp2
  auxvar2%d2den_dp2 = auxvar%d2den_dp2
  auxvar2%mass = auxvar%mass
  auxvar2%dpres_dtime = auxvar%dpres_dtime
  auxvar2%dmass_dtime = auxvar%dmass_dtime

#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
  auxvar2%bc_alpha  = auxvar%bc_alpha
  auxvar2%bc_lambda = auxvar%bc_lambda
#endif

  if (richards_well_coupling > ZERO_INTEGER) auxvar2%well = auxvar%well

end subroutine RichardsAuxVarCopy

! ************************************************************************** !

subroutine RichardsAuxVarCompute(x,auxvar,global_auxvar,material_auxvar, &
                                 richards_parameter, &
                                 characteristic_curves,natural_id, &
                                 update_porosity,option)
  !
  ! Computes auxiliary variables for each grid cell
  !
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  !

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_module
  use Utility_module, only : Equal

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  PetscReal :: x(option%nflowdof)
  type(richards_auxvar_type) :: auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscInt :: natural_id
  PetscBool :: update_porosity

  PetscBool :: saturated
  PetscErrorCode :: ierr
  PetscReal :: pw,dw_kg,dw_mol,sat_pressure,visl
  PetscReal :: kr, ds_dp, dkr_dp
  PetscReal :: dvis_dt, dvis_dp
  PetscReal :: dw_dp, dw_dt, hw_dp
  PetscReal :: dkr_sat
  PetscReal :: aux(1)

  type(richards_parameter_type), pointer :: richards_parameter

  ierr = 0
  global_auxvar%sat = 0.d0
  global_auxvar%den = 0.d0
  global_auxvar%den_kg = 0.d0

  auxvar%kvr = 0.d0
  auxvar%visc = 0.d0
  auxvar%kr = 0.d0

  kr = 0.d0

  global_auxvar%pres = x(1)
  global_auxvar%temp = option%flow%reference_temperature

  call RichardsAuxPorosity(auxvar,material_auxvar,global_auxvar, &
                       update_porosity,richards_parameter,option)

  if (richards_well_coupling == RICHARDS_FULLY_IMPLICIT_WELL) then
    ! This is an initialization hack:
    if(x(RICHARDS_PRESSURE_DOF) /= x(RICHARDS_WELL_DOF))then
      auxvar%well%bh_p = x(RICHARDS_WELL_DOF)
    endif
  endif

  ! For a very large negative liquid pressure (e.g. -1.d18), the capillary
  ! pressure can go near infinite, resulting in ds_dp being < 1.d-40 below
  ! and flipping the cell to saturated, when it is really far from saturated.
  ! The large negative liquid pressure is then passed to the EOS causing it
  ! to blow up.  Therefore, we truncate to the max capillary pressure here.
  auxvar%pc = min(option%flow%reference_pressure - global_auxvar%pres(1), &
                  characteristic_curves%saturation_function%pcmax)

  if (option%flow%disable_capillarity) then
    auxvar%pc = 0.d0
  endif

!***************  Liquid phase properties **************************
  pw = option%flow%reference_pressure
  ds_dp = 0.d0
  dkr_dp = 0.d0

  if (auxvar%pc > 0.d0) then
#if defined(CLM_PFLOTRAN) || defined(CLM_OFFLINE)
    if (auxvar%bc_alpha > 0.d0) then
      select type(sf => characteristic_curves%saturation_function)
        class is(sat_func_vg_type)
          sf%m     = auxvar%bc_lambda
          sf%alpha = auxvar%bc_alpha
        class is(sat_func_bc_type)
            sf%lambda = auxvar%bc_lambda
            sf%alpha  = auxvar%bc_alpha
        class default
          option%io_buffer = 'CLM-PFLOTRAN only supports ' // &
            'sat_func_vg_type and sat_func_bc_type'
          call PrintErrMsg(option)
      end select

      select type(rpf => characteristic_curves%liq_rel_perm_function)
        class is(rpf_mualem_vg_liq_type)
          rpf%m = auxvar%bc_lambda
        class is(rpf_burdine_bc_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_mualem_bc_liq_type)
          rpf%lambda = auxvar%bc_lambda
        class is(rpf_burdine_vg_liq_type)
          rpf%m = auxvar%bc_lambda
        class default
          option%io_buffer = 'Unsupported LIQUID-REL-PERM-FUNCTION'
          call PrintErrMsg(option)
      end select
    endif
#endif
    saturated = PETSC_FALSE
    call characteristic_curves%saturation_function% &
                               Saturation(auxvar%pc,global_auxvar%sat(1), &
                                          ds_dp,option)
    ! if ds_dp is 0, we consider the cell saturated.
    if (ds_dp < 1.d-40) then
      saturated = PETSC_TRUE
    else
      call characteristic_curves%liq_rel_perm_function% &
                       RelativePermeability(global_auxvar%sat(1), &
                                            kr,dkr_sat,option)
      dkr_dp = ds_dp * dkr_sat
    endif
  else
    saturated = PETSC_TRUE
  endif

  ! the purpose for splitting this condition from the 'else' statement
  ! above is due to SaturationFunctionCompute switching a cell to
  ! saturated to prevent unstable (potentially infinite) derivatives when
  ! capillary pressure is very small
  if (saturated) then
    auxvar%pc = 0.d0
    global_auxvar%sat = 1.d0
    kr = 1.d0
    pw = max(global_auxvar%pres(1),pw)
  endif

  if (.not.option%flow%density_depends_on_salinity) then
    call EOSWaterDensity(global_auxvar%temp,pw,dw_kg,dw_mol, &
                         dw_dp,dw_dt,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                          'Error in RichardsAuxVarCompute->EOSWaterDensity')
    endif
    ! may need to compute dpsat_dt to pass to VISW
    call EOSWaterSaturationPressure(global_auxvar%temp,sat_pressure,ierr)
    !geh: 0.d0 passed in for derivative of pressure w/respect to temp
    call EOSWaterViscosity(global_auxvar%temp,pw,sat_pressure,0.d0, &
                           visl,dvis_dt,dvis_dp,ierr)
  else
    if (option%iflag == RICHARDS_UPDATE_FOR_FIXED_ACCUM) then
      ! For the computation of fixed accumulation term use NaCl
      ! value, m_nacl(2), from the previous time step.
      aux(1) = global_auxvar%m_nacl(2)
    else
      ! Use NaCl value for the current time step, m_nacl(1), for computing
      ! the accumulation term
      aux(1) = global_auxvar%m_nacl(1)
    endif
    call EOSWaterDensityExt(global_auxvar%temp,pw,aux, &
                            dw_kg,dw_mol,dw_dp,dw_dt,ierr)
    if (ierr /= 0) then
      call PrintMsgByCell(option,natural_id, &
                          'Error in RichardsAuxVarCompute->EOSWaterDensityExt')
    endif
    call EOSWaterViscosityExt(global_auxvar%temp,pw,sat_pressure,0.d0,aux, &
                              visl,dvis_dt,dvis_dp,ierr)
  endif
  if (.not.saturated) then !kludge since pw is constant in the unsat zone
    dvis_dp = 0.d0
    dw_dp = 0.d0
    hw_dp = 0.d0
  endif

  ! richards_density_kmol_to_kg = 1 or FMWH2O (no other value)
  if (richards_density_kmol_to_kg > 1.1d0) then
    ! cannot convert her as some EOS use kg as the origial density
    global_auxvar%den = dw_mol
    auxvar%dden_dp = dw_dp
  else
    global_auxvar%den = dw_kg
    auxvar%dden_dp = dw_dp*FMWH2O
  endif
  global_auxvar%den_kg = dw_kg
  auxvar%dsat_dp = ds_dp
  auxvar%kr = kr  ! stored solely for output purposes and use in the wellbore
  auxvar%kvr = kr/visl
  auxvar%visc = visl
  auxvar%dkvr_dp = dkr_dp/visl - kr/(visl*visl)*dvis_dp

  if (size(global_auxvar%sat) > 1) then
    global_auxvar%sat(2) = 1.d0 - global_auxvar%sat(1)
  endif

end subroutine RichardsAuxVarCompute

! ************************************************************************** !
subroutine RichardsAuxVarCompute2ndOrderDeriv(rich_auxvar,global_auxvar, &
                                              material_auxvar, &
                                              richards_parameter, &
                                              characteristic_curves, &
                                              option)

  ! Computes 2nd order derivatives auxiliary variables for each grid cell
  !
  ! Author: Gautam Bisht, Satish Karra
  ! Date: 07/02/18, 06/26/19
  !

  use Option_module
  use Global_Aux_module

  use EOS_Water_module
  use Characteristic_Curves_module
  use Characteristic_Curves_Common_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  class(characteristic_curves_type) :: characteristic_curves
  type(richards_auxvar_type) :: rich_auxvar, rich_auxvar_pert
  type(global_auxvar_type) :: global_auxvar, global_auxvar_pert
  type(material_auxvar_type) :: material_auxvar
  type(material_auxvar_type) :: material_auxvar_pert
  type(richards_parameter_type), pointer :: richards_parameter

  PetscReal :: x(option%nflowdof), x_pert(option%nflowdof), pert
  PetscInt :: ideriv

  rich_auxvar%d2sat_dp2 = 0.d0
  rich_auxvar%d2den_dp2 = 0.d0

  call GlobalAuxVarInit(global_auxvar_pert,option)
  call MaterialAuxVarInit(material_auxvar_pert,option)
  call RichardsAuxVarCopy(rich_auxvar,rich_auxvar_pert,option)
  call GlobalAuxVarCopy(global_auxvar,global_auxvar_pert,option)
  call MaterialAuxVarCopy(material_auxvar,material_auxvar_pert,option)
  x(1) = global_auxvar%pres(1)

  ideriv = 1
  pert = max(dabs(x(ideriv)*perturbation_tolerance),0.1d0)
  x_pert = x
  if (x_pert(ideriv) < option%flow%reference_pressure) pert = -1.d0*pert
  x_pert(ideriv) = x_pert(ideriv) + pert

  option%iflag = RICHARDS_UPDATE_FOR_PERTURBATION
  call RichardsAuxVarCompute(x_pert(1),rich_auxvar_pert,global_auxvar_pert, &
                       material_auxvar_pert, &
                       richards_parameter, &
                       characteristic_curves, &
                       -999, PETSC_TRUE, &
                       option)

  rich_auxvar%d2den_dp2 = (rich_auxvar_pert%dden_dp - rich_auxvar%dden_dp)/pert
  if (rich_auxvar%pc > 0.d0) then
    call characteristic_curves%saturation_function% &
                               D2SatDP2(rich_auxvar%pc, &
                                          rich_auxvar%d2sat_dp2,option)
  endif
  option%iflag = UNINITIALIZED_INTEGER

  call GlobalAuxVarStrip(global_auxvar_pert)
  call MaterialAuxVarStrip(material_auxvar_pert)

end subroutine RichardsAuxVarCompute2ndOrderDeriv

! ************************************************************************** !

subroutine RichardsAuxPorosity(auxvar,material_auxvar,global_auxvar, &
                         update_porosity,richards_parameter,option)
  !
  ! Calculates the update to porosity for RICHARDS
  !
  ! Author: Jumanah Al Kubaisy
  ! Date: 03/16/26
  !

  use Material_Aux_module
  use Global_Aux_module
  use Option_module
  use Geomechanics_Linear_Aux_module

  implicit none

  type(richards_auxvar_type) :: auxvar
  type(material_auxvar_type) :: material_auxvar
  type(global_auxvar_type) :: global_auxvar
  PetscBool :: update_porosity
  type(richards_parameter_type) :: richards_parameter
  type(option_type) :: option

  type(geomech_linear_parameter_type), pointer :: geomech_param

  PetscReal :: cell_pressure
  PetscReal :: biot_coeff
  PetscReal :: youngs_mod, poissons_ratio
  PetscReal :: dr_bulk_modulus
  PetscReal :: por

  PetscInt :: id_press_0
  PetscInt :: id_vstrain_0, id_vstrain, id_press_mech
  PetscInt :: id_porosity_mech, id_porosity_flow
  PetscInt :: mat_id
  PetscReal :: porosity_0, C1, C2, press_0
  PetscReal :: vstrain_0, vstrain, del_vstrain
  PetscReal :: press_mech

  geomech_param => richards_parameter%geomech_parameter

  if (option%iflag /= RICHARDS_UPDATE_FOR_BOUNDARY) then
    if (update_porosity) then
      material_auxvar%porosity = material_auxvar%porosity_base
      material_auxvar%dporosity_dp = 0.d0
      if (soil_compressibility_index > 0) then
        cell_pressure = global_auxvar%pres(1)
        call MaterialCompressSoil(material_auxvar,cell_pressure, &
                                  auxvar%effective_porosity, &
                                  auxvar%dpor_dp)
      else
        auxvar%effective_porosity = material_auxvar%porosity_base
        auxvar%dpor_dp =  0.d0
      endif
    else
      auxvar%effective_porosity = material_auxvar%porosity_base
      auxvar%dpor_dp = 0.d0
    endif

    if (associated(geomech_param)) then
      select case(option%geomechanics%flow_coupling)
        case(GEOMECH_TWO_WAY_COUPLED)
          select case(option%geomechanics%split_scheme)
            case(GEOMECH_FIXED_STRESS_SPLIT)
              mat_id = material_auxvar%id
              ! get geomech and flow material properties
              youngs_mod = geomech_param%youngs_modulus(mat_id)
              poissons_ratio = geomech_param%poissons_ratio(mat_id)
              biot_coeff = geomech_param%biot_coeff(mat_id)
              porosity_0 = material_auxvar%porosity_0
              ! 3D bulk modulus
              dr_bulk_modulus = youngs_mod / &
                                (3.d0 * (1.d0 - (2.d0 * poissons_ratio)))
              ! C1 constant
              C1 = (biot_coeff-porosity_0)*(1.d0 - biot_coeff)/dr_bulk_modulus
              ! C2 constant
              C2 = biot_coeff**2/dr_bulk_modulus
              ! get stored values
              id_press_0 = geomech_param%press_0_id
              id_vstrain_0 = geomech_param%vol_strain_0_id
              id_vstrain = geomech_param%vol_strain_id
              id_press_mech = geomech_param%stored_pressure_id
              id_porosity_mech = geomech_param%stored_porosity_id
              id_porosity_flow = geomech_param%flow_porosity_id

              press_0 = global_auxvar%parameters(id_press_0)
              press_mech = global_auxvar%parameters(id_press_mech)
              vstrain_0 = global_auxvar%parameters(id_vstrain_0)
              vstrain = global_auxvar%parameters(id_vstrain)
              if (option%iflag == RICHARDS_UPDATE_FOR_FIXED_ACCUM) then
                por = global_auxvar%parameters(id_porosity_mech)
                auxvar%effective_porosity = por
              elseif (option%iflag == RICHARDS_UPDATE_FOR_ACCUM) then
                ! delta vstrain
                del_vstrain = vstrain - vstrain_0
                ! Burghardt 2017 paper
                por = porosity_0 + (biot_coeff * del_vstrain) + &
                      ( C1 * ( press_mech - press_0 )) + &
                      ( (C2 + C1) * ( global_auxvar%pres(1) - press_mech ))
                ! store new (mass conserved) flow porosity for threshold check
                global_auxvar%parameters(id_porosity_flow) = por

                auxvar%effective_porosity = por
                ! Burghardt 2017 paper
                auxvar%dpor_dp = (biot_coeff**2)/dr_bulk_modulus + &
                    ((biot_coeff-porosity_0)*(1.d0-biot_coeff))/dr_bulk_modulus
              endif
          end select
      end select
    else
      material_auxvar%dporosity_dp = UNINITIALIZED_DOUBLE
    endif

    if (option%iflag /= RICHARDS_UPDATE_FOR_PERTURBATION) then
      material_auxvar%porosity = auxvar%effective_porosity
    endif
  endif

end subroutine RichardsAuxPorosity

! ************************************************************************** !

subroutine AuxVarDestroy(auxvar)
  !
  ! Deallocates a richards auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !

  implicit none

  type(richards_auxvar_type) :: auxvar

  if (associated(auxvar%well)) then
    deallocate(auxvar%well)
    nullify(auxvar%well)
  endif

end subroutine AuxVarDestroy

! ************************************************************************** !

subroutine RichardsAuxDestroy(aux)
  !
  ! Deallocates a richards auxiliary object
  !
  ! Author: Glenn Hammond
  ! Date: 02/14/08
  !
  use Utility_module, only : DeallocateArray

  implicit none

  type(richards_type), pointer :: aux
  PetscInt :: iaux

  if (.not.associated(aux)) return

  if (associated(aux%auxvars)) then
    do iaux = 1, aux%num_aux
      call AuxVarDestroy(aux%auxvars(iaux))
    enddo
    deallocate(aux%auxvars)
  endif
  nullify(aux%auxvars)
  if (associated(aux%auxvars_bc)) then
    do iaux = 1, aux%num_aux_bc
      call AuxVarDestroy(aux%auxvars_bc(iaux))
    enddo
    deallocate(aux%auxvars_bc)
  endif
  nullify(aux%auxvars_bc)
  if (associated(aux%auxvars_ss)) then
    do iaux = 1, aux%num_aux_ss
      call AuxVarDestroy(aux%auxvars_ss(iaux))
    enddo
    deallocate(aux%auxvars_ss)
  endif
  nullify(aux%auxvars_ss)

  if (associated(aux%richards_parameter)) then
    ! solely nullify geomech_parameter as it is destroyed elsewhere
    nullify(aux%richards_parameter%geomech_parameter)
  endif
  nullify(aux%richards_parameter)

  call MatrixZeroingDestroy(aux%matrix_zeroing)

#ifdef BUFFER_MATRIX
  if (associated(aux%matrix_buffer)) then
    call MatrixBufferDestroy(aux%matrix_buffer)
  endif
  nullify(aux%matrix_buffer)
#endif

  deallocate(aux)
  nullify(aux)

end subroutine RichardsAuxDestroy

end module Richards_Aux_module

