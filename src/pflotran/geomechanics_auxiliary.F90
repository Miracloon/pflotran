module Geomechanics_Auxiliary_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Geomechanics_Global_Aux_module
  use Geomechanics_Linear_Aux_module
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: geomech_auxiliary_type
    type(geomech_global_type), pointer :: Global
    type(geomech_linear_type), pointer :: Linear
  end type geomech_auxiliary_type

  public :: GeomechAuxInit, &
            GeomechAuxDestroy

contains

! ************************************************************************** !

subroutine GeomechAuxInit(aux)
  !
  ! Nullifies pointers in geomech auxiliary type
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13

  implicit none

  type(geomech_auxiliary_type) :: aux

  nullify(aux%Global)
  nullify(aux%Linear)

end subroutine GeomechAuxInit

! ************************************************************************** !

subroutine GeomechAuxDestroy(aux)
  !
  ! Strips a geomech auxiliary type
  !
  ! Author: Satish Karra, LANL
  ! Date: 06/14/13

  implicit none

  type(geomech_auxiliary_type) :: aux

  call GeomechGlobalAuxDestroy(aux%Global)
  call GeomechLinearAuxDestroy(aux%Linear)

end subroutine GeomechAuxDestroy

end module Geomechanics_Auxiliary_module
