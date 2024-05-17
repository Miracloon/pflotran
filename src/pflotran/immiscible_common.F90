module Immiscible_Common_module

#include "petsc/finclude/petscsys.h"

  use Immiscible_Aux_module
  use Global_Aux_module

  use PFLOTRAN_Constants_module
  use petscsys

  implicit none

  private

#define SOLUTE_TRANSPORT

! Cutoff parameters
  PetscReal, parameter :: eps       = 1.d-8
  PetscReal, parameter :: floweps   = 0.d0

  procedure(FluxDummy), pointer :: XXFlux => null()
  procedure(BCFluxDummy), pointer :: XXBCFlux => null()

  interface
    subroutine FluxDummy(immiscible_auxvar_up,global_auxvar_up, &
                         material_auxvar_up, &
                         immiscible_auxvar_dn,global_auxvar_dn, &
                         material_auxvar_dn, &
                         area, dist, &
                         immiscible_parameter, &
                         option,v_darcy,Res,Jup,Jdn, &
                         calculate_derivatives, &
                         debug_connection)
      use petscsys
      use Immiscible_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(option_type) :: option
      PetscReal :: v_darcy(IMMIS_MAX_DOF)
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(immiscible_parameter_type) :: immiscible_parameter
      PetscReal :: Res(IMMIS_MAX_DOF)
      PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine FluxDummy
    subroutine BCFluxDummy(ibndtype,auxvar_mapping,auxvars, &
                           immiscible_auxvar_up,global_auxvar_up, &
                           immiscible_auxvar_dn,global_auxvar_dn, &
                           material_auxvar_dn, &
                           area,dist, &
                           immiscible_parameter, &
                           option,v_darcy,Res,Jdn, &
                           calculate_derivatives, &
                           debug_connection)
      use petscsys
      use Immiscible_Aux_module
      use Global_Aux_module
      use Option_module
      use Material_Aux_module
      implicit none
      type(option_type) :: option
      PetscInt :: ibndtype(2)
      PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
      PetscReal :: auxvars(:) ! from aux_real_var array
      type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(material_auxvar_type) :: material_auxvar_dn
      PetscReal :: area
      PetscReal :: dist(-1:3)
      type(immiscible_parameter_type) :: immiscible_parameter
      PetscReal :: v_darcy(IMMIS_MAX_DOF)
      PetscReal :: Res(IMMIS_MAX_DOF)
      PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
      PetscBool :: calculate_derivatives
      PetscBool :: debug_connection
    end subroutine BCFluxDummy
  end interface

  public :: ImmiscibleAccumulation, &
            ImmiscibleFlux, &
            ImmiscibleBCFlux, &
            ImmiscibleAccumDerivative, &
            XXFluxDerivative, &
            XXBCFluxDerivative, &
            ImmiscibleSrcSinkDerivative

  public :: XXFlux, &
            XXBCFlux

contains

! ************************************************************************** !

subroutine ImmiscibleAccumulation(immiscible_auxvar,global_auxvar, &
                                  material_auxvar,option,Res,Jac, &
                                  calculate_derivatives)
  !
  ! Computes the non-fixed portion of the accumulation
  ! term for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscInt :: iphase

  Res = 0.d0
  Jac = 0.d0

  ! accumulation term units = kg comp/s
  ! Res[kg comp/sec] = den[kg comp/m^3 phase] * sat[m^3 phase/m^3 void] *
  !                    por[m^3 void/m^3 bulk] * vol[m^3 bulk] / dt[sec]
  do iphase = 1, 2
    Res(iphase) = immiscible_auxvar%den_kg(iphase) * &
                  immiscible_auxvar%sat(iphase)

  enddo
  Res(:) = Res(:) * material_auxvar%porosity * &
                    material_auxvar%volume / option%flow_dt

  if (calculate_derivatives) then
    stop 'Implement derivative in ImmiscibleAccumulation'
  endif

end subroutine ImmiscibleAccumulation

! ************************************************************************** !

subroutine ImmiscibleFlux(immiscible_auxvar_up,global_auxvar_up, &
                          material_auxvar_up, &
                          immiscible_auxvar_dn,global_auxvar_dn, &
                          material_auxvar_dn, &
                          area, dist,immiscible_parameter,option,v_darcy, &
                          Res,Jup,Jdn,calculate_derivatives, &
                          debug_connection)
  !
  ! Computes the internal flux terms for the residual based on harmonic
  ! intrinsic permeability
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use Connection_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: v_darcy(IMMIS_MAX_DOF)
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscReal :: dist_gravity  ! distance along gravity vector
  PetscReal :: dist_up, dist_dn
  PetscReal :: upweight

  PetscReal :: perm_ave_over_dist
  PetscReal :: perm_up, perm_dn
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: density_ave
  PetscReal :: mobility, q
  PetscReal :: dq_dpup, dq_dpdn
  PetscReal :: dq_dKup, dq_dKdn
  PetscInt :: iphase

  Res = 0.d0
  Jup = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpup = 0.d0
  dq_dpdn = 0.d0
  dq_dKup = 0.d0
  dq_dKdn = 0.d0
  delta_pressure = 0.d0

  call ConnectionCalculateDistances(dist,option%gravity,dist_up,dist_dn, &
                                    dist_gravity,upweight)

  call PermeabilityTensorToScalar(material_auxvar_up,dist,perm_up)
  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  perm_ave_over_dist = perm_up * perm_dn / (dist_up*perm_dn + dist_dn*perm_up)

  do iphase = 1, 2

    if (immiscible_auxvar_up%kr(iphase) + &
        immiscible_auxvar_dn%kr(iphase) > floweps) then

      density_ave = 0.5d0*(immiscible_auxvar_up%den_kg(iphase) + &
                           immiscible_auxvar_dn%den_kg(iphase))
      gravity_term = density_ave * dist_gravity
      delta_pressure = immiscible_auxvar_up%pres(iphase) - &
                       immiscible_auxvar_dn%pres(iphase) + &
                       gravity_term
      if (delta_pressure >= 0.d0) then
        mobility = immiscible_auxvar_up%kr(iphase) / &
                   immiscible_auxvar_up%mu(iphase)
      else
        mobility = immiscible_auxvar_dn%kr(iphase) / &
                   immiscible_auxvar_dn%mu(iphase)
      endif

      if (mobility > floweps) then
        ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
        !                    dP[Pa]]
        v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
        ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area
        ! Res[m^3 liquid/sec]
        Res(iphase) = Res(iphase) + q * density_ave

        if (calculate_derivatives) then
          stop 'Implement derivative in ImmiscibleFlux'
        endif
      endif
    endif
    if (debug_connection) then
      write(*,'(9x,"IFF(res,kr,dp): ",8es12.4)') Res(iphase), &
                                                 mobility, delta_pressure
    endif
  enddo

end subroutine ImmiscibleFlux

! ************************************************************************** !

subroutine ImmiscibleBCFlux(ibndtype,auxvar_mapping,auxvars, &
                            immiscible_auxvar_up,global_auxvar_up, &
                            immiscible_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            area,dist,immiscible_parameter,option,v_darcy, &
                            Res,Jdn,calculate_derivatives,debug_connection)
  !
  ! Computes the boundary flux terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use String_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(2)
  PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy(IMMIS_MAX_DOF)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives
  PetscBool :: debug_connection

  PetscInt :: bc_type
  PetscInt :: idof
  PetscReal :: perm_ave_over_dist
  PetscReal :: dist_gravity
  PetscReal :: delta_pressure
  PetscReal :: gravity_term
  PetscReal :: kr, q
  PetscReal :: perm_dn
  PetscReal :: boundary_pressure
  PetscReal :: density_ave
  PetscReal :: mobility
  PetscBool :: derivative_toggle
  PetscReal :: dq_dpdn
  PetscInt :: iphase

  Res = 0.d0
  Jdn = 0.d0
  v_darcy = 0.d0
  q = 0.d0
  dq_dpdn = 0.d0
  delta_pressure = 0.d0

  call PermeabilityTensorToScalar(material_auxvar_dn,dist,perm_dn)

  do iphase = 1, 2

    if (immiscible_auxvar_up%kr(iphase) + &
        immiscible_auxvar_dn%kr(iphase) > floweps) then

      kr = 0.d0
      bc_type = ibndtype(iphase)
      derivative_toggle = PETSC_TRUE
      select case(bc_type)
        ! figure out the direction of flow
        case(DIRICHLET_BC,DIRICHLET_SEEPAGE_BC,DIRICHLET_CONDUCTANCE_BC, &
            HYDROSTATIC_BC,HYDROSTATIC_SEEPAGE_BC,HYDROSTATIC_CONDUCTANCE_BC)
          if (immiscible_auxvar_up%kr(iphase) + &
              immiscible_auxvar_dn%kr(iphase) > floweps) then
            ! dist(0) = scalar - magnitude of distance
            ! gravity = vector(3)
            ! dist(1:3) = vector(3) - unit vector
            dist_gravity = dist(0) * dot_product(option%gravity,dist(1:3))

            perm_ave_over_dist = perm_dn / dist(0)

            density_ave = 0.5d0*(immiscible_auxvar_up%den_kg(iphase) + &
                                immiscible_auxvar_dn%den_kg(iphase))
            gravity_term = density_ave * dist_gravity
            boundary_pressure = immiscible_auxvar_up%pres(iphase)
            delta_pressure = boundary_pressure - &
                            immiscible_auxvar_dn%pres(iphase) + &
                            gravity_term
            if (delta_pressure >= 0.d0) then
              mobility = immiscible_auxvar_up%kr(iphase) / &
                        immiscible_auxvar_up%mu(iphase)
            else
              mobility = immiscible_auxvar_dn%kr(iphase) / &
                        immiscible_auxvar_dn%mu(iphase)
            endif

            ! v_darcy[m/sec] = perm[m^2] / dist[m] * kr[-] / mu[Pa-sec]
            !                    dP[Pa]]
            v_darcy(iphase) = perm_ave_over_dist * mobility * delta_pressure
          endif

        case(NEUMANN_BC)
          idof = auxvar_mapping(IMMIS_GAS_SATURATION_INDEX+iphase)
          if (dabs(auxvars(idof)) > floweps) then
            v_darcy(iphase) = auxvars(idof)
          endif
          if (v_darcy(iphase) > 0.d0) then
            density_ave = immiscible_auxvar_up%den_kg(iphase)
          else
            density_ave = immiscible_auxvar_up%den_kg(iphase)
          endif
          derivative_toggle = PETSC_FALSE
        case default
          option%io_buffer = &
            'Boundary condition type (' // trim(StringWrite(bc_type)) // &
            ') not recognized in ImmiscibleBCFlux phase loop.'
          call PrintErrMsg(option)
      end select
      if (dabs(v_darcy(iphase)) > 0.d0 .or. mobility > 0.d0) then
        ! q[m^3 liquid/sec] = v_darcy[m/sec] * area[m^2]
        q = v_darcy(iphase) * area
        ! Res[m^3 liquid/sec]
        Res(iphase) = Res(iphase) + q * density_ave
        if (calculate_derivatives .and. derivative_toggle) then
          ! derivative toggle takes care of the NEUMANN side
          stop 'Implement derivative in ImmiscibleBCFlux'
        endif
      endif
      if (debug_connection) then
        write(*,'(9x,"BFF(res,kr,dp): ",8es12.4)') Res(iphase), &
                                                  kr, delta_pressure
      endif
    endif
  enddo

end subroutine ImmiscibleBCFlux

! ************************************************************************** !

subroutine ImmiscibleSrcSink(option,source_sink,scaling_factor, &
                             immiscible_auxvar,global_auxvar,material_auxvar, &
                             ss_flow_vol_flux,Res,Jdn, &
                             calculate_derivatives)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use EOS_Water_module
  use EOS_Gas_module
  use Coupler_module

  implicit none

  type(option_type) :: option
  type(coupler_type) :: source_sink
  PetscReal :: scaling_factor
  type(immiscible_auxvar_type) :: immiscible_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux(:)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: calculate_derivatives

  PetscInt :: iphase
  PetscReal :: qsrc_kg(IMMIS_MAX_DOF)

  qsrc_kg = source_sink%flow_condition%general%rate%dataset%rarray(:)

  Res = 0.d0
  Jdn = 0.d0

  select case(source_sink%flow_condition%general%rate%itype)
    case(MASS_RATE_SS)
    case(SCALED_MASS_RATE_SS)
      qsrc_kg = qsrc_kg*scaling_factor
    case default
      option%io_buffer = 'src_sink_type not supported in ImmiscibleSrcSink'
      call PrintErrMsg(option)
  end select

  do iphase = 1, 2
    ! Res[kg/sec]
    Res(iphase) = qsrc_kg(iphase)
    ss_flow_vol_flux(iphase) = qsrc_kg(iphase) / &
                               immiscible_auxvar%den_kg(iphase)
    if (calculate_derivatives) then
    endif
  enddo

end subroutine ImmiscibleSrcSink

! ************************************************************************** !

subroutine ImmiscibleAccumDerivative(immiscible_auxvar,global_auxvar, &
                                material_auxvar, &
                                option,Res,Jac)
  !
  ! Computes derivatives of the accumulation
  ! term for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  type(option_type) :: option
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  call ImmiscibleAccumulation(immiscible_auxvar(ZERO_INTEGER), &
                         global_auxvar, &
                         material_auxvar, &
                         option,Res,Jac, &
                         .not.immis_numerical_derivatives)

  if (immis_numerical_derivatives) then
    do idof = 1, option%nflowdof
      call ImmiscibleAccumulation(immiscible_auxvar(idof), &
                             global_auxvar, &
                             material_auxvar, &
                             option,res_pert,Jdum, &
                             PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jac(ieq,idof) = (res_pert(ieq)-Res(ieq))/immiscible_auxvar(idof)%pert
      enddo
    enddo

  endif

end subroutine ImmiscibleAccumDerivative

! ************************************************************************** !

subroutine XXFluxDerivative(immiscible_auxvar_up,global_auxvar_up, &
                            material_auxvar_up, &
                            immiscible_auxvar_dn,global_auxvar_dn, &
                            material_auxvar_dn, &
                            area, dist,immiscible_parameter,option,v_darcy, &
                            Res,Jup,Jdn, &
                            debug_connection)
  !
  ! Computes the derivatives of the internal flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module

  implicit none

  type(immiscible_auxvar_type) :: immiscible_auxvar_up(0:), immiscible_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
  type(option_type) :: option
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy(IMMIS_MAX_DOF)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jup(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: v_darcy_dum(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jup = 0.d0
  Jdn = 0.d0

  call XXFlux(immiscible_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
              material_auxvar_up, &
              immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
              material_auxvar_dn, &
              area,dist,immiscible_parameter,option,v_darcy, &
              Res,Jup,Jdn, &
              .not.immis_numerical_derivatives, &
              debug_connection)

  if (immis_numerical_derivatives) then
    ! upgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(immiscible_auxvar_up(idof),global_auxvar_up, &
                  material_auxvar_up, &
                  immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jup(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_up(idof)%pert
      enddo
    enddo
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXFlux(immiscible_auxvar_up(ZERO_INTEGER),global_auxvar_up, &
                  material_auxvar_up, &
                  immiscible_auxvar_dn(idof),global_auxvar_dn, &
                  material_auxvar_dn, &
                  area,dist,immiscible_parameter,option,v_darcy_dum, &
                  res_pert,Jdum,Jdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXFluxDerivative

! ************************************************************************** !

subroutine XXBCFluxDerivative(ibndtype,auxvar_mapping,auxvars, &
                              immiscible_auxvar_up, &
                              global_auxvar_up, &
                              immiscible_auxvar_dn,global_auxvar_dn, &
                              material_auxvar_dn, &
                              area,dist,immiscible_parameter,option,v_darcy, &
                              Res,Jdn,debug_connection)
  !
  ! Computes the derivatives of the boundary flux terms
  ! for the Jacobian
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !

  use Option_module
  use Material_Aux_module

  implicit none

  type(option_type) :: option
  PetscInt :: ibndtype(2)
  PetscInt :: auxvar_mapping(IMMIS_MAX_INDEX)
  PetscReal :: auxvars(:) ! from aux_real_var array
  type(immiscible_auxvar_type) :: immiscible_auxvar_up, immiscible_auxvar_dn(0:)
  type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
  type(material_auxvar_type) :: material_auxvar_dn
  PetscReal :: area
  PetscReal :: dist(-1:3)
  type(immiscible_parameter_type) :: immiscible_parameter
  PetscReal :: v_darcy(IMMIS_MAX_DOF)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jdn(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscBool :: debug_connection

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: v_darcy_dum(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jdn = 0.d0

  call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                immiscible_auxvar_up,global_auxvar_up, &
                immiscible_auxvar_dn(ZERO_INTEGER),global_auxvar_dn, &
                material_auxvar_dn, &
                area,dist,immiscible_parameter,option,v_darcy, &
                Res,Jdn,.not.immis_numerical_derivatives, &
                debug_connection)

  if (immis_numerical_derivatives) then
    ! downgradient derivatives
    do idof = 1, option%nflowdof
      call XXBCFlux(ibndtype,auxvar_mapping,auxvars, &
                    immiscible_auxvar_up,global_auxvar_up, &
                    immiscible_auxvar_dn(idof),global_auxvar_dn, &
                    material_auxvar_dn, &
                    area,dist,immiscible_parameter,option,v_darcy_dum, &
                    res_pert,Jdum,PETSC_FALSE,PETSC_FALSE)
      do ieq = 1, option%nflowdof
        Jdn(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                        immiscible_auxvar_dn(idof)%pert
      enddo
    enddo
  endif

end subroutine XXBCFluxDerivative

! ************************************************************************** !

subroutine ImmiscibleSrcSinkDerivative(option,source_sink,scaling_factor, &
                                       immiscible_auxvar,global_auxvar, &
                                       material_auxvar,ss_flow_vol_flux, &
                                       Res,Jac)
  !
  ! Computes the source/sink terms for the residual
  !
  ! Author: Glenn Hammond
  ! Date: 05/16/24
  !
  use Option_module
  use Material_Aux_module
  use Coupler_module

  implicit none

  type(option_type) :: option
  type(coupler_type) :: source_sink
  PetscReal :: scaling_factor
  type(immiscible_auxvar_type) :: immiscible_auxvar(0:)
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  PetscReal :: ss_flow_vol_flux(:)
  PetscReal :: Res(IMMIS_MAX_DOF)
  PetscReal :: Jac(IMMIS_MAX_DOF,IMMIS_MAX_DOF)

  PetscReal :: res_pert(IMMIS_MAX_DOF)
  PetscReal :: dummy_real(IMMIS_MAX_DOF)
  PetscReal :: Jdum(IMMIS_MAX_DOF,IMMIS_MAX_DOF)
  PetscInt :: idof, ieq

  Jac = 0.d0
  ! unperturbed immiscible_auxvars value
  call ImmiscibleSrcSink(option,source_sink,scaling_factor, &
                         immiscible_auxvar(ZERO_INTEGER),global_auxvar, &
                         material_auxvar,ss_flow_vol_flux, &
                         Res,Jac,.not.immis_numerical_derivatives)

  if (immis_numerical_derivatives) then
    ! perturbed immiscible_auxvars values
    do idof = 1, option%nflowdof
      call ImmiscibleSrcSink(option,source_sink,scaling_factor, &
                             immiscible_auxvar(idof),global_auxvar, &
                             material_auxvar, &
                             dummy_real, &
                             res_pert,Jdum,PETSC_FALSE)
      do ieq = 1, option%nflowdof
      Jac(ieq,idof) = (res_pert(ieq)-Res(ieq)) / &
                      immiscible_auxvar(idof)%pert
      enddo
    enddo
  endif

end subroutine ImmiscibleSrcSinkDerivative

end module Immiscible_Common_module
