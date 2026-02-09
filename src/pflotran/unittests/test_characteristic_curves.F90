module Test_Characteristic_Curves_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use pFUnit_mod
  use Characteristic_Curves_module
  use Characteristic_Curves_Base_module
  use Characteristic_Curves_Common_module
  use Characteristic_Curves_WIPP_module
  use Option_module
  use PFLOTRAN_Constants_module, only : MAXSTRINGLENGTH

  implicit none


  public :: Test_Characteristic_Curves

  interface Test_Characteristic_Curves
     module procedure newTest_Characteristic_Curves
  end interface Test_Characteristic_Curves

! ************************************************************************** !
!  @TestCase
  type, extends(TestCase) :: Test_Characteristic_Curves
      type(option_type), pointer :: option
      class(characteristic_curves_type), pointer :: cc_const
      class(characteristic_curves_type), pointer :: cc_bcb
      class(characteristic_curves_type), pointer :: cc_bcm
      class(characteristic_curves_type), pointer :: cc_vgm
      class(characteristic_curves_type), pointer :: cc_vgb
      class(characteristic_curves_type), pointer :: cc_vgt
      class(characteristic_curves_type), pointer :: cc_lb
      class(characteristic_curves_type), pointer :: cc_lm
      class(characteristic_curves_type), pointer :: cc_mk3
      class(characteristic_curves_type), pointer :: cc_mk4
      class(characteristic_curves_type), pointer :: cc_krp1
      class(characteristic_curves_type), pointer :: cc_krp2
      class(characteristic_curves_type), pointer :: cc_krp3
      class(characteristic_curves_type), pointer :: cc_krp4
      class(characteristic_curves_type), pointer :: cc_krp5
      class(characteristic_curves_type), pointer :: cc_krp8
      class(characteristic_curves_type), pointer :: cc_krp9
      class(characteristic_curves_type), pointer :: cc_krp11
      class(characteristic_curves_type), pointer :: cc_krp12
      class(characteristic_curves_type), pointer :: cc_modbc
      procedure(runMethod), pointer :: userMethod => null()
    contains
      procedure :: setUp
      procedure :: tearDown
      procedure :: runMethod
  end type Test_Characteristic_Curves

contains

! ************************************************************************** !

  function newTest_Characteristic_Curves(name, userMethod) result(test)

    implicit none

    character(len=*), intent(in) :: name
    procedure(runMethod) :: userMethod

    type(Test_Characteristic_Curves) :: test

    call test%setName(name)
    test%userMethod => userMethod

  end function newTest_Characteristic_Curves

! ************************************************************************** !

  subroutine setUp(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this
    character(len=MAXSTRINGLENGTH) :: error_string

    this%option => OptionCreate()

  ! Setting up the characteristic curve for Brooks Corey Burdine
    this%cc_bcb => CharacteristicCurvesCreate()
    this%cc_bcb%saturation_function => SFBCCreate()
    this%cc_bcb%saturation_function%Sr = 0.2d0
    this%cc_bcb%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        sf%lambda = 0.7d0
        sf%alpha = 9.869d-6
      class default
        print *, 'not BC type in setUp'
    end select
    this%cc_bcb%liq_rel_perm_function => RPFBurdineBCLiqCreate()
    this%cc_bcb%liq_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcb%liq_rel_perm_function)
      class is(rpf_Burdine_BC_liq_type)
        sf%lambda = 0.7d0
      class default
        print *, 'not Burdine BC Liq type'
    end select
    this%cc_bcb%gas_rel_perm_function => RPFBurdineBCGasCreate()
    this%cc_bcb%gas_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcb%gas_rel_perm_function)
      class is(rpf_Burdine_BC_gas_type)
        sf%lambda = 0.7d0
        sf%Srg = 1.d-5
      class default
        print *, 'not Burdine BC Gas type'
    end select
    error_string = 'pFUnit Brooks Corey saturation function'
    call this%cc_bcb%saturation_function%SetupPolynomials(this%option, &
                                                          error_string)

  ! Setting up the characteristic curve for Brooks Corey Mualem
    this%cc_bcm => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_bcb
    this%cc_bcm%liq_rel_perm_function => RPFMualemBCLiqCreate()
    this%cc_bcm%liq_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcm%liq_rel_perm_function)
      class is(rpf_Mualem_BC_liq_type)
        sf%lambda = 0.7d0
      class default
        print *, 'not Mualem BC Liq type'
    end select
    this%cc_bcm%gas_rel_perm_function => RPFMualemBCGasCreate()
    this%cc_bcm%gas_rel_perm_function%Sr = 0.2d0
    select type(sf=>this%cc_bcm%gas_rel_perm_function)
      class is(rpf_Mualem_BC_gas_type)
        sf%lambda = 0.7d0
        sf%Srg = 1.d-5
      class default
        print *, 'not Mualem BC Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten Mualem
    this%cc_vgm => CharacteristicCurvesCreate()
    this%cc_vgm%saturation_function => SFVGCreate()
    this%cc_vgm%saturation_function%Sr = 0.143d0
    this%cc_vgm%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_vgm%saturation_function)
      class is(sat_func_vg_type)
        sf%m = 0.527d0
        sf%alpha = 5.1054d-5
      class default
        print *, 'not VG type in setUp'
    end select
    this%cc_vgm%liq_rel_perm_function => RPFMualemVGLiqCreate()
    this%cc_vgm%liq_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%liq_rel_perm_function)
      class is(rpf_mualem_vg_liq_type)
        sf%m = 0.527d0
      class default
        print *, 'not Mualem VG Liq type'
    end select
    this%cc_vgm%gas_rel_perm_function => RPFMualemVGGasCreate()
    this%cc_vgm%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgm%gas_rel_perm_function)
      class is(rpf_Mualem_VG_gas_type)
        sf%m = 0.527d0
        sf%Srg = 0.d0
      class default
        print *, 'not Mualem VG Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten Burdine
    this%cc_vgb => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_vgm
    this%cc_vgb%liq_rel_perm_function => RPFBurdineVGLiqCreate()
    this%cc_vgb%liq_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgb%liq_rel_perm_function)
      class is(rpf_Burdine_VG_liq_type)
        sf%m = 0.527d0
      class default
        print *, 'not Burdine VG Liq type'
    end select
    this%cc_vgb%gas_rel_perm_function => RPFBurdineVGGasCreate()
    this%cc_vgb%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgb%gas_rel_perm_function)
      class is(rpf_Burdine_VG_gas_type)
        sf%m = 0.527d0
        sf%Srg = 0.01d0
      class default
        print *, 'not Burdine VG Gas type'
    end select

  ! Setting up the characteristic curve for van Genuchten TOUGH2 IRP7
    this%cc_vgt => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_vgm
    this%cc_vgt%gas_rel_perm_function => RPFTOUGH2IRP7GasCreate()
    this%cc_vgt%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_vgt%gas_rel_perm_function)
      class is(rpf_TOUGH2_IRP7_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not TOUGH2 IRP7 Gas type'
    end select

  ! Setting up the characteristic curve for Linear Mualem
    this%cc_lm => CharacteristicCurvesCreate()
    this%cc_lm%saturation_function => SFLinearCreate()
    this%cc_lm%saturation_function%Sr = 0.143d0
    this%cc_lm%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_lm%saturation_function)
      class is(sat_func_linear_type)
        sf%alpha = 5.1054d-5
      class default
        print *, 'not Linear type in setUp'
    end select
    this%cc_lm%liq_rel_perm_function => RPFMualemLinearLiqCreate()
    this%cc_lm%liq_rel_perm_function%Sr = 0.143d0
    this%cc_lm%gas_rel_perm_function => RPFMualemLinearGasCreate()
    this%cc_lm%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_lm%gas_rel_perm_function)
      class is(rpf_Mualem_Linear_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not Mualem Linear Gas type'
    end select

  ! Setting up the characteristic curve for Linear Burdine
    this%cc_lb => CharacteristicCurvesCreate()
    ! do not set up saturation function as it is tested with cc_lm
    this%cc_lb%liq_rel_perm_function => RPFBurdineLinearLiqCreate()
    this%cc_lb%liq_rel_perm_function%Sr = 0.143d0
    this%cc_lb%gas_rel_perm_function => RPFBurdineLinearGasCreate()
    this%cc_lb%gas_rel_perm_function%Sr = 0.143d0
    select type(sf=>this%cc_lb%gas_rel_perm_function)
      class is(rpf_Burdine_Linear_gas_type)
        sf%Srg = 0.01d0
      class default
        print *, 'not Burdine Linear Gas type'
    end select

  ! Setting up the characteristic curve for 3-param modified Kosugi
    this%cc_mk3 => CharacteristicCurvesCreate()
    this%cc_mk3%saturation_function => SFmKCreate()
    this%cc_mk3%saturation_function%Sr = 1.53D-1
    this%cc_mk3%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_mk3%saturation_function)
      class is(sat_func_mk_type)
        sf%nparam = 3
        sf%sigmaz = 3.37D-1
        sf%muz = -6.3D0
        sf%rmax = 2.52D-3
      class default
        print *, 'no mK3 saturation function type in setUp'
    end select
    this%cc_mk3%liq_rel_perm_function => RPFmKLiqCreate()
    this%cc_mk3%liq_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk3%liq_rel_perm_function)
      class is(rpf_mK_liq_type)
        sf%sigmaz = 3.37D-1
      class default
        print *, 'no mK3 Liq rel perm type in setUp'
    end select
    this%cc_mk3%gas_rel_perm_function => RPFmKGasCreate()
    this%cc_mk3%gas_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk3%gas_rel_perm_function)
      class is(rpf_mk_gas_type)
        sf%sigmaz = 3.37D-1
        sf%Srg = 1.0D-3
      class default
        print *, 'no mK3 Gas rel perm type in setUp'
    end select

  ! Setting up the characteristic curve for 4-param modified Kosugi
    this%cc_mk4 => CharacteristicCurvesCreate()
    this%cc_mk4%saturation_function => SFmKCreate()
    this%cc_mk4%saturation_function%Sr = 1.53D-1
    this%cc_mk4%saturation_function%pcmax = 0.999d8
    select type(sf=>this%cc_mk4%saturation_function)
      class is(sat_func_mk_type)
        sf%nparam = 4
        sf%sigmaz = 3.37D-1
        sf%muz = -6.3D0
        sf%rmax = 2.52D-3
        sf%r0 = 1.07D-4
      class default
        print *, 'no mK4 saturation function type in setUp'
    end select
    this%cc_mk4%liq_rel_perm_function => RPFmKLiqCreate()
    this%cc_mk4%liq_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk4%liq_rel_perm_function)
      class is(rpf_mK_liq_type)
        sf%sigmaz = 3.37D-1
      class default
        print *, 'no mK4 Liq rel perm type in setUp'
    end select
    this%cc_mk4%gas_rel_perm_function => RPFmKGasCreate()
    this%cc_mk4%gas_rel_perm_function%Sr = 1.53D-1
    select type(sf=>this%cc_mk4%gas_rel_perm_function)
      class is(rpf_mk_gas_type)
        sf%sigmaz = 3.37D-1
        sf%Srg = 1.0D-3
      class default
        print *, 'no mK4 Gas rel perm type in setUp'
    end select

  ! Setting up the characteristic curve constant capillary pressure
  ! saturation function
    this%cc_const => CharacteristicCurvesCreate()
    this%cc_const%saturation_function => SFConstantCreate()
    select type(sf=>this%cc_const%saturation_function)
      class is(sat_func_constant_type)
        sf%constant_capillary_pressure = 1.d5
        sf%constant_saturation = 0.5d0
      class default
        print *, 'no constant saturation function type in setUp'
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP1
    this%cc_krp1 => CharacteristicCurvesCreate()
    this%cc_krp1%saturation_function => SFKRP1Create()
    this%cc_krp1%saturation_function%Sr = 0.20d0
    this%cc_krp1%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp1%saturation_function)
      class is(sat_func_KRP1_type)
        sf%Srg = 0.25d0
        sf%m = 0.60d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 1
    end select
    this%cc_krp1%liq_rel_perm_function => RPFKRP1LiqCreate()
    this%cc_krp1%liq_rel_perm_function%Sr = 0.20d0
    select type(sf=>this%cc_krp1%liq_rel_perm_function)
      class is(rpf_KRP1_liq_type)
        sf%m = 0.60d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp1%gas_rel_perm_function => RPFKRP1GasCreate()
    this%cc_krp1%gas_rel_perm_function%Sr = 0.25d0
    select type(sf=>this%cc_krp1%gas_rel_perm_function)
      class is(rpf_KRP1_gas_type)
        sf%m = 0.60d0
        sf%Srg = 0.25d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP2
    this%cc_krp2 => CharacteristicCurvesCreate()
    this%cc_krp2%saturation_function => SFKRP2Create()
    this%cc_krp2%saturation_function%Sr = 0.10d0
    this%cc_krp2%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp2%saturation_function)
      class is(sat_func_krp2_type)
        sf%lambda = 0.25d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp2%liq_rel_perm_function => RPFKRP2LiqCreate()
    this%cc_krp2%liq_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp2%liq_rel_perm_function)
      class is(rpf_KRP2_liq_type)
        sf%lambda = 0.25d0
    end select
    this%cc_krp2%gas_rel_perm_function => RPFKRP2GasCreate()
    this%cc_krp2%gas_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp2%gas_rel_perm_function)
      class is(rpf_KRP2_gas_type)
        sf%lambda = 0.25d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP3
    this%cc_krp3 => CharacteristicCurvesCreate()
    this%cc_krp3%saturation_function => SFKRP3Create()
    this%cc_krp3%saturation_function%Sr = 0.15d0
    this%cc_krp3%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp3%saturation_function)
      class is(sat_func_krp3_type)
        sf%Srg = 0.25d0
        sf%lambda = 0.35d0
        sf%alpha = 5.0d-3
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp3%liq_rel_perm_function => RPFKRP3LiqCreate()
    this%cc_krp3%liq_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp3%liq_rel_perm_function)
      class is(rpf_KRP3_liq_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp3%gas_rel_perm_function => RPFKRP3GasCreate()
    this%cc_krp3%gas_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp3%gas_rel_perm_function)
      class is(rpf_KRP3_gas_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP4
    this%cc_krp4 => CharacteristicCurvesCreate()
    this%cc_krp4%saturation_function => SFKRP4Create()
    this%cc_krp4%saturation_function%Sr = 0.15d0
    this%cc_krp4%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp4%saturation_function)
      class is(sat_func_krp4_type)
        sf%Srg = 0.25d0
        sf%lambda = 0.35d0
        sf%alpha = 5.0d-3
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp4%liq_rel_perm_function => RPFKRP4LiqCreate()
    this%cc_krp4%liq_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp4%liq_rel_perm_function)
      class is(rpf_KRP4_liq_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select
    this%cc_krp4%gas_rel_perm_function => RPFKRP4GasCreate()
    this%cc_krp4%gas_rel_perm_function%Sr = 0.15d0
    select type(sf=>this%cc_krp4%gas_rel_perm_function)
      class is(rpf_KRP4_gas_type)
        sf%lambda = 0.35d0
        sf%Srg = 0.25d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP5
    this%cc_krp5 => CharacteristicCurvesCreate()
    this%cc_krp5%saturation_function => SFKRP5Create()
    this%cc_krp5%saturation_function%Sr = 0.05d0
    this%cc_krp5%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp5%saturation_function)
      class is(sat_func_krp5_type)
        sf%alpha = 5.0d-4
        sf%Srg = 0.20d0
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp5%liq_rel_perm_function => RPFKRP5LiqCreate()
    this%cc_krp5%liq_rel_perm_function%Sr = 0.05d0
    select type(sf=>this%cc_krp5%liq_rel_perm_function)
      class is(rpf_KRP5_liq_type)
        sf%Srg = 0.20d0
    end select
    this%cc_krp5%gas_rel_perm_function => RPFKRP5GasCreate()
    this%cc_krp5%gas_rel_perm_function%Sr = 0.05d0
    select type(sf=>this%cc_krp5%gas_rel_perm_function)
      class is(rpf_KRP5_gas_type)
        sf%Srg = 0.20d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP8
    this%cc_krp8 => CharacteristicCurvesCreate()
    this%cc_krp8%saturation_function => SFKRP8Create()
    this%cc_krp8%saturation_function%Sr = 0.10d0
    this%cc_krp8%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp8%saturation_function)
      class is(sat_func_KRP8_type)
        sf%Srg = 0.25d0
        sf%m = 0.45d0
        sf%alpha = 5.0d-5
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp8%liq_rel_perm_function => RPFKRP8LiqCreate()
    this%cc_krp8%liq_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp8%liq_rel_perm_function)
      class is(rpf_KRP8_liq_type)
        sf%m = 0.45d0
    end select
    this%cc_krp8%gas_rel_perm_function => RPFKRP8GasCreate()
    this%cc_krp8%gas_rel_perm_function%Sr = 0.10d0
    select type(sf=>this%cc_krp8%gas_rel_perm_function)
      class is(rpf_KRP8_gas_type)
        sf%m = 0.45d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP9
    this%cc_krp9 => CharacteristicCurvesCreate()
    this%cc_krp9%saturation_function => SFKRP9Create()
    this%cc_krp9%saturation_function%Sr = 0.05d0
    this%cc_krp9%saturation_function%pcmax = 0.999d7
    this%cc_krp9%liq_rel_perm_function => RPFKRP9LiqCreate()
    this%cc_krp9%liq_rel_perm_function%Sr = 0.05d0
    this%cc_krp9%gas_rel_perm_function => RPFKRP9GasCreate()
    this%cc_krp9%gas_rel_perm_function%Sr = 0.05d0

  ! Setting up the characteristic curve for BRAGFLO_KRP11
    this%cc_krp11 => CharacteristicCurvesCreate()
    this%cc_krp11%saturation_function => SFKRP11Create()
    this%cc_krp11%saturation_function%Sr = 0.08d0
    this%cc_krp11%saturation_function%pcmax = 0.999d7
    this%cc_krp11%liq_rel_perm_function => RPFKRP11LiqCreate()
    this%cc_krp11%liq_rel_perm_function%Sr = 0.08d0
    select type(sf=>this%cc_krp11%liq_rel_perm_function)
      class is(rpf_KRP11_liq_type)
        sf%tolc = 0.10
        sf%Srg = 0.18d0
    end select
    this%cc_krp11%gas_rel_perm_function => RPFKRP11GasCreate()
    this%cc_krp11%gas_rel_perm_function%Sr = 0.08d0
    select type(sf=>this%cc_krp11%gas_rel_perm_function)
      class is(rpf_KRP11_gas_type)
        sf%tolc = 0.10
        sf%Srg = 0.18d0
    end select

  ! Setting up the characteristic curve for BRAGFLO_KRP12
    this%cc_krp12 => CharacteristicCurvesCreate()
    this%cc_krp12%saturation_function => SFKRP12Create()
    this%cc_krp12%saturation_function%Sr = 0.12d0
    this%cc_krp12%saturation_function%pcmax = 0.999d7
    select type(sf=>this%cc_krp12%saturation_function)
      class is(sat_func_krp12_type)
        sf%lambda = 0.2d0
        sf%alpha = 5.0d-4
        sf%s_min = 0.015d0
        sf%s_effmin = 0.005d0
        sf%ignore_permeability = PETSC_TRUE
        sf%kpc = 2
    end select
    this%cc_krp12%liq_rel_perm_function => RPFKRP12LiqCreate()
    this%cc_krp12%liq_rel_perm_function%Sr = 0.12d0
    select type(sf=>this%cc_krp12%liq_rel_perm_function)
      class is(rpf_KRP12_liq_type)
        sf%lambda = 0.20d0
        sf%Srg = 0.20d0
    end select
    this%cc_krp12%gas_rel_perm_function => RPFKRP12GasCreate()
    this%cc_krp12%gas_rel_perm_function%Sr = 0.12d0
    select type(sf=>this%cc_krp12%gas_rel_perm_function)
      class is(rpf_KRP12_gas_type)
        sf%lambda = 0.20d0
        sf%Srg = 0.20d0
    end select

  ! Setting up the characteristic curve for modified Brooks Corey
  ! relative permeability
    this%cc_modbc => CharacteristicCurvesCreate()
    this%cc_modbc%liq_rel_perm_function => RPFModBrooksCoreyLiqCreate()
    this%cc_modbc%liq_rel_perm_function%Sr = 0.153d0
    this%cc_modbc%liq_rel_perm_function%Srg = 0.d0
    select type(rpf=>this%cc_modbc%liq_rel_perm_function)
      class is(rpf_mod_Brooks_Corey_liq_type)
        rpf%kr_max = 0.888d0
        rpf%n = 3.2d0
      class default
        print *, 'not modified BC Liq type'
    end select

    this%cc_modbc%gas_rel_perm_function => RPFModBrooksCoreyGasCreate()
    this%cc_modbc%gas_rel_perm_function%Sr = 0.143d0
    this%cc_modbc%gas_rel_perm_function%Srg = 0.1d0
    select type(rpf=>this%cc_modbc%gas_rel_perm_function)
      class is(rpf_mod_Brooks_Corey_gas_type)
        rpf%kr_max = 0.777d0
        rpf%n = 2.2d0
      class default
        print *, 'not modified BC Gas type'
    end select

  end subroutine setUp

! ************************************************************************** !

  subroutine tearDown(this)

    implicit none
    class (Test_Characteristic_Curves), intent(inout) :: this

    call OptionDestroy(this%option)
    call CharacteristicCurvesDestroy(this%cc_bcb)

  end subroutine tearDown

! ************************************************************************** !

  subroutine runMethod(this)
    implicit none
    class (Test_Characteristic_Curves), intent(inout) :: this
    call this%userMethod()
  end subroutine runMethod

! ************************************************************************** !

!  @Test
  subroutine testSF_BC_SetupPolynomials(this)

    implicit none

    class(Test_characteristic_Curves), intent(inout) :: this

    PetscReal :: values(4)
    PetscReal, parameter :: tolerance = 1.d-8
    PetscInt :: i
    character(len=128) :: string

    ! pressure polynomial
    values = [-4.6122570934041036d0, 1.4087313882738163d-4, &
              -1.1098865178492886d-9, 2.6190024815166203d-15]
    do i = 1, 4
      write(string,*) i
      string = 'Brooks-Corey-Burdine pressure polynomial coefficient #' // &
               trim(adjustl(string))
#line 599 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%pres_poly%coefficients(i), dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 599) )
  if (anyExceptions()) return
#line 600 "test_characteristic_curves.pf"
    enddo

    ! saturation spline
    values = [-83508464.603000879d0, 173197055.36650354d0, &
              -89688590.763502657d0, 0.d0]
    do i = 1, 3
      write(string,*) i
      string = 'Brooks-Corey-Burdine saturation spline coefficient #' // &
               trim(adjustl(string))
#line 609 "test_characteristic_curves.pf"
  call assertEqual(values(i), this%cc_bcb%saturation_function%sat_poly%coefficients(i),  dabs(values(i))*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 609) )
  if (anyExceptions()) return
#line 610 "test_characteristic_curves.pf"
    enddo

  end subroutine testSF_BC_SetupPolynomials

! ************************************************************************** !

!  @Test
  subroutine testsf_Brooks_Corey(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) below the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        capillary_pressure = 0.94d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure below polynomial fit'
    value = 1.d0
#line 650 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 650) )
  if (anyExceptions()) return
#line 651 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure below polynomial fit'
    value = 1.d0
#line 654 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 654) )
  if (anyExceptions()) return
#line 655 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure below polynomial fit'
    value = 0.d0
#line 658 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 658) )
  if (anyExceptions()) return
#line 659 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure below polynomial fit'
    value = 0.d0
#line 662 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 662) )
  if (anyExceptions()) return
#line 663 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) within the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        capillary_pressure = 0.96d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure within polynomial fit'
    value = 0.99971176979312304d0
#line 681 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 681) )
  if (anyExceptions()) return
#line 682 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure within polynomial fit'
    value = 0.99789158871529349d0
#line 685 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 685) )
  if (anyExceptions()) return
#line 686 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure within polynomial fit'
    value = 5.6675690490728353d-7
#line 689 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 689) )
  if (anyExceptions()) return
#line 690 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure within polynomial fit'
    value = 4.1422137957785640d-006
#line 693 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 693) )
  if (anyExceptions()) return
#line 694 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) above the polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        capillary_pressure = 1.06d0/sf%alpha
      class default
        print *, 'not bc type in testsf_Brooks_Corey'
    end select
    call this%cc_bcb%saturation_function% &
         Saturation(capillary_pressure, &
                    liquid_saturation, dsat_pres,this%option)
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Brooks-Corey-Burdine saturation as a function of capillary &
             &pressure above polynomial fit'
    value = 0.96802592722174041d0
#line 712 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 712) )
  if (anyExceptions()) return
#line 713 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine relative permeability as a function of &
             &capillary pressure above polynomial fit'
    value = 0.78749164071142996d0
#line 716 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 716) )
  if (anyExceptions()) return
#line 717 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of saturation as a function &
             &of capillary pressure above polynomial fit'
    value = 5.0054278424773111d-6
#line 720 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 720) )
  if (anyExceptions()) return
#line 721 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of relative permeability &
             &as a function of pressure above polynomial fit'
    value = 3.0060561800889172d-005
#line 724 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 724) )
  if (anyExceptions()) return
#line 725 "test_characteristic_curves.pf"

  end subroutine testsf_Brooks_Corey

! ************************************************************************** !

!  @Test
  subroutine testcp_Brooks_Corey(this)

    implicit none

    class (Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation) well within polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        liquid_saturation = 1.00001d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation barely within polynomial fit'
    value = 54.068777590990067d0
#line 758 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 758) )
  if (anyExceptions()) return
#line 759 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation barely within polynomial fit'
    value = -7.7231957793806121d6
#line 762 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 762) )
  if (anyExceptions()) return
#line 763 "test_characteristic_curves.pf"

    ! capillary pressure = f(saturation) slightly within polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        liquid_saturation = 1.04d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation well within polynomial fit'
    value = 106436.99642977261d0
#line 777 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 777) )
  if (anyExceptions()) return
#line 778 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation well within polynomial fit'
    value = -1.9672548074671443d5
#line 781 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 781) )
  if (anyExceptions()) return
#line 782 "test_characteristic_curves.pf"

    ! capillary pressure = f(saturation) above polynomial fit
    select type(sf=>this%cc_bcb%saturation_function)
      class is(sat_func_bc_type)
        liquid_saturation = 1.06d0**(-sf%lambda)
      class default
        print *, 'not bc type in capillary pressure Brooks Corey'
    end select
    call this%cc_bcb%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Brooks-Corey capillary pressure as a function of &
             &liquid saturation above within polynomial fit'
    value = 109024.42772683989d0
#line 796 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 796) )
  if (anyExceptions()) return
#line 797 "test_characteristic_curves.pf"
    string = 'Brooks-Corey derivative of capillary pressure as a function of &
             &liquid saturation above within polynomial fit'
    value = -2.0492439614025093d5
#line 800 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 800) )
  if (anyExceptions()) return
#line 801 "test_characteristic_curves.pf"

  end subroutine testcp_Brooks_Corey

! ************************************************************************** !

!  @Test
  subroutine testrpf_BC_Burdine(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_bcb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Burdine liquid relative permeability as a &
             &function of liquid saturation'
    value = 3.1991918327000197d-3
#line 828 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 828) )
  if (anyExceptions()) return
#line 829 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 6.2460411971762310d-2
#line 832 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 832) )
  if (anyExceptions()) return
#line 833 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.38173220142506209d0
#line 841 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 841) )
  if (anyExceptions()) return
#line 842 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.6412199910843073d0
#line 845 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 845) )
  if (anyExceptions()) return
#line 846 "test_characteristic_curves.pf"

  end subroutine testrpf_BC_Burdine

! ************************************************************************** !

!  @Test
  subroutine testrpf_BC_Mualem(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: relative_permeability
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_bcm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 5.2242583862629442d-3
#line 873 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 873) )
  if (anyExceptions()) return
#line 874 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 9.3290328326124022d-2
#line 877 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 877) )
  if (anyExceptions()) return
#line 878 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_bcm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Brooks-Corey-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.65126653365343257d0
#line 886 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 886) )
  if (anyExceptions()) return
#line 887 "test_characteristic_curves.pf"
    string = 'Brooks-Corey-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.7243442005604193d0
#line 890 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, 1.d-8, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 890) )
  if (anyExceptions()) return
#line 891 "test_characteristic_curves.pf"

  end subroutine testrpf_BC_Mualem

! ************************************************************************** !

!  @Test
  subroutine testsf_van_Genuchten(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) at low capillary pressure
    capillary_pressure = 10.d0
    call this%cc_vgm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99999995045230206d0
#line 925 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 925) )
  if (anyExceptions()) return
#line 926 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99957025105913566d0
#line 929 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 929) )
  if (anyExceptions()) return
#line 930 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 1.0475199529417896d-8
#line 933 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 933) )
  if (anyExceptions()) return
#line 934 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at low capillary pressure'
    value = 4.7878857031474202d-005
#line 937 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 937) )
  if (anyExceptions()) return
#line 938 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) at high capillary pressure
    select type(sf=>this%cc_vgm%saturation_function)
      class is(sat_func_vg_type)
        capillary_pressure = 10.d0/sf%alpha
      class default
        print *, 'not vg type in testsf_van_Genuchten'
    end select
    call this%cc_vgm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.20862404282784081d0
#line 956 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 956) )
  if (anyExceptions()) return
#line 957 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 4.4900562293186444d-6
#line 960 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 960) )
  if (anyExceptions()) return
#line 961 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 3.7043838142841442d-7
#line 964 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 964) )
  if (anyExceptions()) return
#line 965 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of relative permeability &
             &as a function of capillary pressure at high capillary pressure'
    value = 1.0903614584398361d-010
#line 968 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 968) )
  if (anyExceptions()) return
#line 969 "test_characteristic_curves.pf"

  end subroutine testsf_van_Genuchten

! ************************************************************************** !

!  @Test
  subroutine testcp_van_Genuchten(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgm%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'van Genuchten-Mualem capillary pressure as a function of &
             &saturation'
    value = 38910.985405751228d0
#line 997 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 997) )
  if (anyExceptions()) return
#line 998 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of capillary pressure as a &
              &function of saturation'
    value = -1.2074621994359563d5
#line 1001 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1001) )
  if (anyExceptions()) return
#line 1002 "test_characteristic_curves.pf"

  end subroutine testcp_van_Genuchten

! ************************************************************************** !

!  @Test
  subroutine testrpf_van_Genuchten_Mualem(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'van Genuchten-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 7.1160141309814171d-3
#line 1029 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1029) )
  if (anyExceptions()) return
#line 1030 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.9580035202641822d-2
#line 1033 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1033) )
  if (anyExceptions()) return
#line 1034 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'van Genuchten-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.61184154078016839d0
#line 1042 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1042) )
  if (anyExceptions()) return
#line 1043 "test_characteristic_curves.pf"
    string = 'van Genuchten-Mualem derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4149310375033495d0
#line 1046 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1046) )
  if (anyExceptions()) return
#line 1047 "test_characteristic_curves.pf"

  end subroutine testrpf_van_Genuchten_Mualem

! ************************************************************************** !

!  @Test
  subroutine testrpf_van_Genuchten_Burdine(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_vgb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'van Genuchten-Burdine liquid relative permeability as a &
             &function of liquid saturation'
    value = 1.8220963608953099d-2
#line 1074 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1074) )
  if (anyExceptions()) return
#line 1075 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 0.20400586616752553d0
#line 1078 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1078) )
  if (anyExceptions()) return
#line 1079 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_vgb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'van Genuchten-Burdine gas relative permeability as a &
             &function of liquid saturation'
    value = 0.29870096712333277d0
#line 1087 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1087) )
  if (anyExceptions()) return
#line 1088 "test_characteristic_curves.pf"
    string = 'van Genuchten-Burdine derivative of gas relative &
             &permeability as a function of liquid saturation'
    value = -1.4207000510364920d0
#line 1091 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1091) )
  if (anyExceptions()) return
#line 1092 "test_characteristic_curves.pf"

  end subroutine testrpf_van_Genuchten_Burdine

! ************************************************************************** !

!  @Test
  subroutine testrpf_TOUGH2_IRP7_gas(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0

    ! gas relative permeability = f(saturation)
    call this%cc_vgt%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'TOUGH2 IRP7 gas relative permeability as a &
             &function of liquid saturation'
    value = 0.27522069402853439d0
#line 1121 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1121) )
  if (anyExceptions()) return
#line 1122 "test_characteristic_curves.pf"
    string = 'TOUGH2 IRP7 derivative of gas relative permeability as a &
             &function of liquid saturation'
    value = -1.4564360410147879d0
#line 1125 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1125) )
  if (anyExceptions()) return
#line 1126 "test_characteristic_curves.pf"
  end subroutine testrpf_TOUGH2_IRP7_gas

! ************************************************************************** !

!  @Test
  subroutine testsf_Linear(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure) at low capillary pressure
    capillary_pressure = 10.d0
    call this%cc_lm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'Linear-Mualem saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0001679766584224d0
#line 1159 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1159) )
  if (anyExceptions()) return
#line 1160 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0000000000000000d0
#line 1163 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1163) )
  if (anyExceptions()) return
#line 1164 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at low capillary pressure'
    value = 8.5802608854958100d-9
#line 1167 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1167) )
  if (anyExceptions()) return
#line 1168 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at low capillary pressure'
    value = 0.d0
#line 1171 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1171) )
  if (anyExceptions()) return
#line 1172 "test_characteristic_curves.pf"

    ! saturation = f(capillary_pressure) at high capillary pressure
    select type(sf=>this%cc_lm%saturation_function)
      class is(sat_func_linear_type)
        capillary_pressure = 10.d0/sf%alpha
      class default
        print *, 'not linear type in testsf_Linear'
    end select
    call this%cc_lm%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = 'van Genuchten-Mualem saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.99848743785071759d0
#line 1190 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1190) )
  if (anyExceptions()) return
#line 1191 "test_characteristic_curves.pf"
    string = 'Linear-Mualem relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.21040639641042236d0
#line 1194 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1194) )
  if (anyExceptions()) return
#line 1195 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of saturation as a function &
             &of capillary pressure at high capillary pressure'
    value = 8.5802608854958100d-9
#line 1198 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1198) )
  if (anyExceptions()) return
#line 1199 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of relative permeability as a &
             &function of capillary pressure at high capillary pressure'
    value = 3.7741597839123353d-007
#line 1202 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1202) )
  if (anyExceptions()) return
#line 1203 "test_characteristic_curves.pf"

  end subroutine testsf_Linear

! ************************************************************************** !

!  @Test
  subroutine testcp_Linear(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lm%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = 'Linear capillary pressure as a function of &
             &saturation'
    value = 58292873.507671818d0
#line 1231 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1231) )
  if (anyExceptions()) return
#line 1232 "test_characteristic_curves.pf"
    string = 'Linear derivative of capillary pressure as a function of &
             &saturation'
    value = -1.1654657280764198d8
#line 1235 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1235) )
  if (anyExceptions()) return
#line 1236 "test_characteristic_curves.pf"

  end subroutine testcp_Linear

! ************************************************************************** !

!  @Test
  subroutine testrpf_Linear_Mualem(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lm%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Linear-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 9.8205932575543323d-4
#line 1263 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1263) )
  if (anyExceptions()) return
#line 1264 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 8.6657420154635095d-3
#line 1267 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1267) )
  if (anyExceptions()) return
#line 1268 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lm%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.70258636467899827d0
#line 1276 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1276) )
  if (anyExceptions()) return
#line 1277 "test_characteristic_curves.pf"

  end subroutine testrpf_Linear_Mualem

! ************************************************************************** !

!  @Test
  subroutine testrpf_Linear_Burdine(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_lb%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Linear-Mualem liquid relative permeability as a &
             &function of liquid saturation'
    value = 0.41656942823803966d0
#line 1304 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1304) )
  if (anyExceptions()) return
#line 1305 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = 1.1668611435239207d0
#line 1308 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1308) )
  if (anyExceptions()) return
#line 1309 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_lb%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = 'Linear-Mualem gas relative permeability as a &
             &function of liquid saturation'
    value = 0.57851239669421495d0
#line 1317 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1317) )
  if (anyExceptions()) return
#line 1318 "test_characteristic_curves.pf"
    string = 'Linear-Mualem derivative of liquid relative &
             &permeability as a function of liquid saturation'
    value = -1.1806375442739079d0
#line 1321 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1321) )
  if (anyExceptions()) return
#line 1322 "test_characteristic_curves.pf"

  end subroutine testrpf_Linear_Burdine

! ************************************************************************** !

!  @Test
  subroutine testsf_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal :: dkr_p
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure)
    capillary_pressure = 1.0D+4
    call this%cc_mk3%saturation_function%Saturation(capillary_pressure, &
                                         liquid_saturation, &
                                         dsat_pres,this%option)
    ! relative permeability = f(saturation)
    call this%cc_mk3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, &
                              relative_permeability,dkr_sat,this%option)
    dkr_p = dsat_pres * dkr_sat
    string = '3-parameter modified Kosugi saturation as a function of '//&
         &'capillary pressure'
    value = 0.97481347586415668d0
#line 1357 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1357) )
  if (anyExceptions()) return
#line 1358 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi relative permeability as a '//&
         &'function of capillary pressure'
    value = 0.92508223462292838d0
#line 1361 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1361) )
  if (anyExceptions()) return
#line 1362 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of saturation as '//&
         &'a function of capillary pressure'
    value = 4.7607337955577978D-5
#line 1365 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1365) )
  if (anyExceptions()) return
#line 1366 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of relative perm as'// &
         &' a function of capillary pressure'
    value = 1.2551629198938072D-4
#line 1369 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_p, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1369) )
  if (anyExceptions()) return
#line 1370 "test_characteristic_curves.pf"

  end subroutine testsf_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testcp_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 2.5d-1
    call this%cc_mk3%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = '3-parameter modified Kosugi capillary pressure as a '//&
         &'function of saturation'
    value = 17707.010438883957d0
#line 1398 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1398) )
  if (anyExceptions()) return
#line 1399 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -24500.406539923177d0
#line 1402 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1402) )
  if (anyExceptions()) return
#line 1403 "test_characteristic_curves.pf"

  end subroutine testcp_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testrpf_modified_kosugi_3param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_mk3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '3-parameter modified Kosugi liquid relative permeability'//&
         &' as a function of liquid saturation'
    value = 0.18300229972367688d0
#line 1430 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1430) )
  if (anyExceptions()) return
#line 1431 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0
#line 1434 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1434) )
  if (anyExceptions()) return
#line 1435 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '3-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0
#line 1443 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1443) )
  if (anyExceptions()) return
#line 1444 "test_characteristic_curves.pf"
    string = '3-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0
#line 1447 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1447) )
  if (anyExceptions()) return
#line 1448 "test_characteristic_curves.pf"

  end subroutine testrpf_modified_kosugi_3param

! ************************************************************************** !

!  @Test
  subroutine testcp_modified_kosugi_4param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string
    PetscReal, parameter :: temperature = 25.d0

    ! capillary pressure = f(saturation)
    liquid_saturation = 2.5d-1
    call this%cc_mk4%saturation_function% &
         CapillaryPressure(liquid_saturation, &
                           capillary_pressure, dpc_dsatl, this%option)
    string = '4-parameter modified Kosugi capillary pressure as a '//&
         &'function of saturation'
    value = 15671.942098006928d0
#line 1476 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1476) )
  if (anyExceptions()) return
#line 1477 "test_characteristic_curves.pf"
    string = 'r-parameter modified Kosugi derivative of capillary ' // &
         'pressure as a function of saturation'
    value = -19192.362637570313d0
#line 1480 "test_characteristic_curves.pf"
  call assertEqual(value, dpc_dsatl, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1480) )
  if (anyExceptions()) return
#line 1481 "test_characteristic_curves.pf"

  end subroutine testcp_modified_kosugi_4param

! ************************************************************************** !

!  @Test
  subroutine testrpf_modified_kosugi_4param(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dkr_sat
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! liquid relative permeability = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_mk4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '4-parameter modified Kosugi liquid relative permeability'//&
         &' as a function of liquid saturation'
    value = 0.18300229972367688d0  ! same as 3-param
#line 1508 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1508) )
  if (anyExceptions()) return
#line 1509 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of liquid relative'//&
         &' permeability as a function of liquid saturation'
    value = 0.92477440290277491d0  ! same as 3-param
#line 1512 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1512) )
  if (anyExceptions()) return
#line 1513 "test_characteristic_curves.pf"

    ! gas relative permeability = f(saturation)
    call this%cc_mk4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation, relative_permeability, &
                              dkr_sat, this%option)
    string = '4-parameter modified Kosugi gas relative permeability as'//&
         &' a function of liquid saturation'
    value = 0.35040465818058608d0  ! same as 3-param
#line 1521 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1521) )
  if (anyExceptions()) return
#line 1522 "test_characteristic_curves.pf"
    string = '4-parameter modified Kosugi derivative of gas relative'//&
         &' permeability as a function of liquid saturation'
    value = -1.2770282015052450d0   ! same as 3-param
#line 1525 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1525) )
  if (anyExceptions()) return
#line 1526 "test_characteristic_curves.pf"

  end subroutine testrpf_modified_kosugi_4param

! ************************************************************************** !

!  @Test
  subroutine testsf_constant(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: dpc_dsatl
    PetscReal :: liquid_saturation
    PetscReal :: dsat_pres
    PetscReal :: value
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    ! saturation = f(capillary_pressure)
    capillary_pressure = 1.d5
    call this%cc_const%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'Constant saturation function liquid saturation'
    value = 0.5d0
#line 1553 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1553) )
  if (anyExceptions()) return
#line 1554 "test_characteristic_curves.pf"
    string = 'Constant saturation function liquid saturation derivative'
    value = 0.d0
#line 1556 "test_characteristic_curves.pf"
  call assertEqual(value, dsat_pres, tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1556) )
  if (anyExceptions()) return
#line 1557 "test_characteristic_curves.pf"

    ! capillary_pressure = f(saturation)
    liquid_saturation = 0.5d0
    call this%cc_const%saturation_function% &
                                    CapillaryPressure(liquid_saturation, &
                                                      capillary_pressure, &
                                                      dpc_dsatl, &
                                                      this%option)
    string = 'Constant saturation function capillary pressure'
    value = 1.d5
#line 1567 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1567) )
  if (anyExceptions()) return
#line 1568 "test_characteristic_curves.pf"

  end subroutine testsf_constant

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP1(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.22
    call this%cc_krp1%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP1 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 211645.56118006655d0
#line 1597 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1597) )
  if (anyExceptions()) return
#line 1598 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp1%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP1 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.749999998736306000
#line 1605 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1605) )
  if (anyExceptions()) return
#line 1606 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.10878895018640269d0
#line 1612 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1612) )
  if (anyExceptions()) return
#line 1613 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 5.9765155008000600d-4
#line 1619 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1619) )
  if (anyExceptions()) return
#line 1620 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp1%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP1 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 1628 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1628) )
  if (anyExceptions()) return
#line 1629 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp1%saturation_function)
      class is(sat_func_KRP1_type)
        capillary_pressure = 10.d0/sf%alpha
    end select
    call this%cc_krp1%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP1 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.20000256369851988d0
#line 1639 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1639) )
  if (anyExceptions()) return
#line 1640 "test_characteristic_curves.pf"
    call this%cc_krp1%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 3.1269006686492169d-22
#line 1646 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1646) )
  if (anyExceptions()) return
#line 1647 "test_characteristic_curves.pf"
    call this%cc_krp1%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP1 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 1653 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1653) )
  if (anyExceptions()) return
#line 1654 "test_characteristic_curves.pf"

  end subroutine testsf_KRP1

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP2(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.05
    call this%cc_krp2%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP2 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 9990000.0d0
#line 1683 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1683) )
  if (anyExceptions()) return
#line 1684 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1691 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1691) )
  if (anyExceptions()) return
#line 1692 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1698 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1698) )
  if (anyExceptions()) return
#line 1699 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1705 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1705) )
  if (anyExceptions()) return
#line 1706 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp2%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP2 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 51637.669399930157d0
#line 1714 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1714) )
  if (anyExceptions()) return
#line 1715 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp2%saturation_function)
      class is(sat_func_krp2_type)
        capillary_pressure = 15.d0/sf%alpha
    end select
    call this%cc_krp2%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP2 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.55731947333915333d0
#line 1724 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1724) )
  if (anyExceptions()) return
#line 1725 "test_characteristic_curves.pf"
    call this%cc_krp2%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8310805074531403d-4
#line 1731 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1731) )
  if (anyExceptions()) return
#line 1732 "test_characteristic_curves.pf"
    call this%cc_krp2%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP2 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.24138701885980715d0
#line 1738 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1738) )
  if (anyExceptions()) return
#line 1739 "test_characteristic_curves.pf"

  end subroutine testsf_KRP2

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP3(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.21
    call this%cc_krp3%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP3 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 143937.17953954436d0
#line 1768 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1768) )
  if (anyExceptions()) return
#line 1769 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1775 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1775) )
  if (anyExceptions()) return
#line 1776 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1782 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1782) )
  if (anyExceptions()) return
#line 1783 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1789 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1789) )
  if (anyExceptions()) return
#line 1790 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp3%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP3 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 200.0d0
#line 1798 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1798) )
  if (anyExceptions()) return
#line 1799 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp3%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP3 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1805 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1805) )
  if (anyExceptions()) return
#line 1806 "test_characteristic_curves.pf"
    call this%cc_krp3%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 5.8632915084341688d-9
#line 1812 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1812) )
  if (anyExceptions()) return
#line 1813 "test_characteristic_curves.pf"
    call this%cc_krp3%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP3 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1819 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1819) )
  if (anyExceptions()) return
#line 1820 "test_characteristic_curves.pf"

  end subroutine testsf_KRP3

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP4(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp4%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP4 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 9990000.0d0
#line 1849 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1849) )
  if (anyExceptions()) return
#line 1850 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 1856 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1856) )
  if (anyExceptions()) return
#line 1857 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1863 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1863) )
  if (anyExceptions()) return
#line 1864 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1870 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1870) )
  if (anyExceptions()) return
#line 1871 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp4%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP4 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 152.32289193581161d0
#line 1879 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1879) )
  if (anyExceptions()) return
#line 1880 "test_characteristic_curves.pf"
    capillary_pressure = 1.d5
    call this%cc_krp4%saturation_function%Saturation(capillary_pressure, &
                                   liquid_saturation,dsat_pres,this%option)
    string = 'KRP4 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.21815720128839766d0
#line 1886 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1886) )
  if (anyExceptions()) return
#line 1887 "test_characteristic_curves.pf"
    call this%cc_krp4%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.8180697920854255d-10
#line 1893 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1893) )
  if (anyExceptions()) return
#line 1894 "test_characteristic_curves.pf"
    call this%cc_krp4%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP4 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.78571287226820719d0
#line 1900 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1900) )
  if (anyExceptions()) return
#line 1901 "test_characteristic_curves.pf"

  end subroutine testsf_KRP4

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP5(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp5%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP5 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 8525093.3809598293d0
#line 1930 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1930) )
  if (anyExceptions()) return
#line 1931 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.80014642571085304d0
#line 1938 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1938) )
  if (anyExceptions()) return
#line 1939 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 1945 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1945) )
  if (anyExceptions()) return
#line 1946 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.0d0
#line 1952 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1952) )
  if (anyExceptions()) return
#line 1953 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp5%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP5 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2000.0d0
#line 1961 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1961) )
  if (anyExceptions()) return
#line 1962 "test_characteristic_curves.pf"
    capillary_pressure = 8.d5
    call this%cc_krp5%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP5 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.74007809371245503d0
#line 1969 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1969) )
  if (anyExceptions()) return
#line 1970 "test_characteristic_curves.pf"
    call this%cc_krp5%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.92010412494993998d0
#line 1976 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1976) )
  if (anyExceptions()) return
#line 1977 "test_characteristic_curves.pf"
    call this%cc_krp5%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP5 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 7.9895875050060017d-2
#line 1983 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 1983) )
  if (anyExceptions()) return
#line 1984 "test_characteristic_curves.pf"

  end subroutine testsf_KRP5

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP8(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp8%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP8 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 390287.05422819848d0
#line 2013 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2013) )
  if (anyExceptions()) return
#line 2014 "test_characteristic_curves.pf"
    capillary_pressure = 10.d0
    call this%cc_krp8%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP8 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99813934029313034d0
#line 2021 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2021) )
  if (anyExceptions()) return
#line 2022 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.82967396002016347d0
#line 2028 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2028) )
  if (anyExceptions()) return
#line 2029 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 3.5744473041460765d-4
#line 2035 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2035) )
  if (anyExceptions()) return
#line 2036 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp8%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP8 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 11673.176741409796d0
#line 2044 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2044) )
  if (anyExceptions()) return
#line 2045 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp8%saturation_function)
      class is(sat_func_KRP8_type)
        capillary_pressure = 15.d0/sf%alpha
    end select
    call this%cc_krp8%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP8 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.10220037076344601d0
#line 2055 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2055) )
  if (anyExceptions()) return
#line 2056 "test_characteristic_curves.pf"
    call this%cc_krp8%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 2.4705248145444649d-14
#line 2062 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2062) )
  if (anyExceptions()) return
#line 2063 "test_characteristic_curves.pf"
    call this%cc_krp8%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP8 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.99877541173464723d0
#line 2069 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2069) )
  if (anyExceptions()) return
#line 2070 "test_characteristic_curves.pf"

  end subroutine testsf_KRP8

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP9(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.16
    call this%cc_krp9%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP9 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 6701.4503478003835d0
#line 2099 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2099) )
  if (anyExceptions()) return
#line 2100 "test_characteristic_curves.pf"
    capillary_pressure = 50.d0
    call this%cc_krp9%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP9 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 0.99999644138299326d0
#line 2107 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2107) )
  if (anyExceptions()) return
#line 2108 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.99999998839687010d0
#line 2114 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2114) )
  if (anyExceptions()) return
#line 2115 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.1603129901338605d-8
#line 2121 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2121) )
  if (anyExceptions()) return
#line 2122 "test_characteristic_curves.pf"

    liquid_saturation = 0.81
    call this%cc_krp9%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP9 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 2294.5062171637132d0
#line 2130 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2130) )
  if (anyExceptions()) return
#line 2131 "test_characteristic_curves.pf"
    select type(sf=>this%cc_krp9%saturation_function)
      class is(sat_func_krp9_type)
        capillary_pressure = 8.d5
    end select
    call this%cc_krp9%saturation_function%Saturation(capillary_pressure, &
                                                     liquid_saturation, &
                                                     dsat_pres,this%option)
    string = 'KRP9 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 1.8062138845808708d-7
#line 2141 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2141) )
  if (anyExceptions()) return
#line 2142 "test_characteristic_curves.pf"
    call this%cc_krp9%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.d0
#line 2148 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2148) )
  if (anyExceptions()) return
#line 2149 "test_characteristic_curves.pf"
    call this%cc_krp9%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP9 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.d0
#line 2155 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2155) )
  if (anyExceptions()) return
#line 2156 "test_characteristic_curves.pf"

  end subroutine testsf_KRP9

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP11(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.15
    call this%cc_krp11%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP11 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 0.0d0
#line 2185 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2185) )
  if (anyExceptions()) return
#line 2186 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2193 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2193) )
  if (anyExceptions()) return
#line 2194 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2200 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2200) )
  if (anyExceptions()) return
#line 2201 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2207 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2207) )
  if (anyExceptions()) return
#line 2208 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp11%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP11 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 0.0d0
#line 2216 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2216) )
  if (anyExceptions()) return
#line 2217 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp11%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP11 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 1.0d0
#line 2224 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2224) )
  if (anyExceptions()) return
#line 2225 "test_characteristic_curves.pf"
    call this%cc_krp11%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.0d0
#line 2231 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2231) )
  if (anyExceptions()) return
#line 2232 "test_characteristic_curves.pf"
    call this%cc_krp11%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP11 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.0d0
#line 2238 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2238) )
  if (anyExceptions()) return
#line 2239 "test_characteristic_curves.pf"

  end subroutine testsf_KRP11

! ************************************************************************** !

!  @Test
  subroutine testsf_KRP12(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: capillary_pressure
    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: dsat_pres
    PetscReal :: dpc_dsatl
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=128) :: string

    liquid_saturation = 0.15
    call this%cc_krp12%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP12 capillary pressure as a function of saturation &
             &at low liquid saturation'
    value = 35364351.440591194d0
#line 2268 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2268) )
  if (anyExceptions()) return
#line 2269 "test_characteristic_curves.pf"
    capillary_pressure = 100.d0
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at low capillary pressure'
    value = 1.0d0
#line 2276 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2276) )
  if (anyExceptions()) return
#line 2277 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 1.0d0
#line 2283 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2283) )
  if (anyExceptions()) return
#line 2284 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at low capillary pressure'
    value = 0.d0
#line 2290 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2290) )
  if (anyExceptions()) return
#line 2291 "test_characteristic_curves.pf"

    liquid_saturation = 0.85
    call this%cc_krp12%saturation_function% &
         CapillaryPressure(liquid_saturation,capillary_pressure, &
                           dpc_dsatl,this%option)
    string = 'KRP12 capillary pressure as a function of saturation &
             &at high liquid saturation'
    value = 4547.8850244305122d0
#line 2299 "test_characteristic_curves.pf"
  call assertEqual(value, capillary_pressure, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2299) )
  if (anyExceptions()) return
#line 2300 "test_characteristic_curves.pf"
    capillary_pressure = 1.d6
    call this%cc_krp12%saturation_function%Saturation(capillary_pressure, &
                                                      liquid_saturation, &
                                                      dsat_pres,this%option)
    string = 'KRP12 saturation as a function of capillary &
             &pressure at high capillary pressure'
    value = 0.30276918155781385d0
#line 2307 "test_characteristic_curves.pf"
  call assertEqual(value, liquid_saturation, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2307) )
  if (anyExceptions()) return
#line 2308 "test_characteristic_curves.pf"
    call this%cc_krp12%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 liquid relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 1.3380237859265290d-9
#line 2314 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2314) )
  if (anyExceptions()) return
#line 2315 "test_characteristic_curves.pf"
    call this%cc_krp12%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'KRP12 gas relative permeability as a function of &
             &capillary pressure at high capillary pressure'
    value = 0.53468502607262869d0
#line 2321 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2321) )
  if (anyExceptions()) return
#line 2322 "test_characteristic_curves.pf"

  end subroutine testsf_KRP12

! ************************************************************************** !

!  @Test
  subroutine testcc_KRPModBrooksCorey(this)

    implicit none

    class(Test_Characteristic_Curves), intent(inout) :: this

    PetscReal :: liquid_saturation
    PetscReal :: relative_permeability
    PetscReal :: value
    PetscReal :: dkr_sat
    PetscReal, parameter :: tolerance = 1.d-8
    character(len=MAXSTRINGLENGTH) :: string

    string = ''
    call this%cc_modbc%liq_rel_perm_function%Verify(string,this%option)
    call this%cc_modbc%gas_rel_perm_function%Verify(string,this%option)

    liquid_saturation = 0.13d0
    call this%cc_modbc%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey liquid relative permeability at &
             &very low liquid saturation'
    value = 0.d0
#line 2352 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2352) )
  if (anyExceptions()) return
#line 2353 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey liquid relative permeability derivative at &
             &very low liquid saturation'
    value = 0.d0
#line 2356 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2356) )
  if (anyExceptions()) return
#line 2357 "test_characteristic_curves.pf"

    call this%cc_modbc%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey gas relative permeability at &
             &very low liquid saturation'
    value = 0.777d0
#line 2364 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2364) )
  if (anyExceptions()) return
#line 2365 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey gas relative permeability derivative at &
             &very low liquid saturation'
    value = 0.d0
#line 2368 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2368) )
  if (anyExceptions()) return
#line 2369 "test_characteristic_curves.pf"

    liquid_saturation = 1.d0
    call this%cc_modbc%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey liquid relative permeability at &
             &liquid saturation'
    value = 0.888d0
#line 2377 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2377) )
  if (anyExceptions()) return
#line 2378 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey liquid relative permeability derivative at &
             &very low liquid saturation'
    value = 0.d0
#line 2381 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2381) )
  if (anyExceptions()) return
#line 2382 "test_characteristic_curves.pf"

    call this%cc_modbc%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey gas relative permeability at &
             &liquid saturation'
    value = 0.d0
#line 2389 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2389) )
  if (anyExceptions()) return
#line 2390 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey gas relative permeability derivative at &
             &liquid saturation'
    value = 0.d0
#line 2393 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2393) )
  if (anyExceptions()) return
#line 2394 "test_characteristic_curves.pf"

    liquid_saturation = 0.85d0
    call this%cc_modbc%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey liquid relative permeability at &
             &high liquid saturation'
    value = 0.47591669791650282d0
#line 2402 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2402) )
  if (anyExceptions()) return
#line 2403 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey liquid relative permeability derivative at &
             &high liquid saturation'
    value = 2.1849834050685928d0
#line 2406 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2406) )
  if (anyExceptions()) return
#line 2407 "test_characteristic_curves.pf"

    call this%cc_modbc%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey gas relative permeability at &
             &high liquid saturation'
    value = 1.9685393179233227d-3
#line 2414 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2414) )
  if (anyExceptions()) return
#line 2415 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey gas relative permeability derivative at &
             &high liquid saturation'
    value = -8.6615729988626183d-2
#line 2418 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2418) )
  if (anyExceptions()) return
#line 2419 "test_characteristic_curves.pf"

    liquid_saturation = 0.2d0
    call this%cc_modbc%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey liquid relative permeability at &
             &low liquid saturation'
    value = 8.5094050780208850d-5
#line 2427 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2427) )
  if (anyExceptions()) return
#line 2428 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey liquid relative permeability derivative at &
             &low liquid saturation'
    value = 5.7936374999291110d-3
#line 2431 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2431) )
  if (anyExceptions()) return
#line 2432 "test_characteristic_curves.pf"

    call this%cc_modbc%gas_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey gas relative permeability at &
             &low liquid saturation'
    value = 0.65407232230065315d0
#line 2439 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2439) )
  if (anyExceptions()) return
#line 2440 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey gas relative permeability derivative at &
             &low liquid saturation'
    value = -2.0556558700877670d0
#line 2443 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2443) )
  if (anyExceptions()) return
#line 2444 "test_characteristic_curves.pf"

    ! set Srg > 0 for liquid rel perm function
    liquid_saturation = 0.6d0
    this%cc_modbc%liq_rel_perm_function%Srg = 0.05d0
    call this%cc_modbc%liq_rel_perm_function% &
         RelativePermeability(liquid_saturation,relative_permeability, &
                              dkr_sat,this%option)
    string = 'Modified Brooks Corey liquid relative permeability at &
             &mid liquid saturation and Srg > 0 for liquid'
    value = 0.13955018331922028d0
#line 2454 "test_characteristic_curves.pf"
  call assertEqual(value, relative_permeability, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2454) )
  if (anyExceptions()) return
#line 2455 "test_characteristic_curves.pf"
    string = 'Modified Brooks Corey liquid relative permeability derivative at &
             &mid liquid saturation and Srg > 0'
    value = 0.99901697230761743d0
#line 2458 "test_characteristic_curves.pf"
  call assertEqual(value, dkr_sat, dabs(value)*tolerance, string, &
 & location=SourceLocation( &
 & 'test_characteristic_curves.pf', &
 & 2458) )
  if (anyExceptions()) return
#line 2459 "test_characteristic_curves.pf"

  end subroutine testcc_KRPModBrooksCorey

! ************************************************************************** !

end module Test_Characteristic_Curves_module




function Test_Characteristic_Curves_module_suite() result(suite)
   use pFUnit_mod
   use Test_Characteristic_Curves_module
   implicit none

   type (TestSuite) :: suite
   suite = newTestSuite('Test_Characteristic_Curves_module_suite')

   call suite%addTest(Test_Characteristic_Curves('testSF_BC_SetupPolynomials',testSF_BC_SetupPolynomials))
   call suite%addTest(Test_Characteristic_Curves('testsf_Brooks_Corey',testsf_Brooks_Corey))
   call suite%addTest(Test_Characteristic_Curves('testcp_Brooks_Corey',testcp_Brooks_Corey))
   call suite%addTest(Test_Characteristic_Curves('testrpf_BC_Burdine',testrpf_BC_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testrpf_BC_Mualem',testrpf_BC_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testsf_van_Genuchten',testsf_van_Genuchten))
   call suite%addTest(Test_Characteristic_Curves('testcp_van_Genuchten',testcp_van_Genuchten))
   call suite%addTest(Test_Characteristic_Curves('testrpf_van_Genuchten_Mualem',testrpf_van_Genuchten_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testrpf_van_Genuchten_Burdine',testrpf_van_Genuchten_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testrpf_TOUGH2_IRP7_gas',testrpf_TOUGH2_IRP7_gas))
   call suite%addTest(Test_Characteristic_Curves('testsf_Linear',testsf_Linear))
   call suite%addTest(Test_Characteristic_Curves('testcp_Linear',testcp_Linear))
   call suite%addTest(Test_Characteristic_Curves('testrpf_Linear_Mualem',testrpf_Linear_Mualem))
   call suite%addTest(Test_Characteristic_Curves('testrpf_Linear_Burdine',testrpf_Linear_Burdine))
   call suite%addTest(Test_Characteristic_Curves('testsf_modified_kosugi_3param',testsf_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testcp_modified_kosugi_3param',testcp_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testrpf_modified_kosugi_3param',testrpf_modified_kosugi_3param))
   call suite%addTest(Test_Characteristic_Curves('testcp_modified_kosugi_4param',testcp_modified_kosugi_4param))
   call suite%addTest(Test_Characteristic_Curves('testrpf_modified_kosugi_4param',testrpf_modified_kosugi_4param))
   call suite%addTest(Test_Characteristic_Curves('testsf_constant',testsf_constant))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP1',testsf_KRP1))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP2',testsf_KRP2))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP3',testsf_KRP3))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP4',testsf_KRP4))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP5',testsf_KRP5))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP8',testsf_KRP8))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP9',testsf_KRP9))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP11',testsf_KRP11))
   call suite%addTest(Test_Characteristic_Curves('testsf_KRP12',testsf_KRP12))
   call suite%addTest(Test_Characteristic_Curves('testcc_KRPModBrooksCorey',testcc_KRPModBrooksCorey))

end function Test_Characteristic_Curves_module_suite

