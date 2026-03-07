module Shape_Function_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Gauss_module
  use Grid_Unstructured_Cell_module
  use PFLOTRAN_Constants_module
  use Option_module

  implicit none

  private

  type, public :: shapefunction_type
    PetscInt :: element_type          ! element type
    PetscReal, pointer :: zeta(:)     ! coordinates of point in reference element
    PetscReal, pointer :: N(:)        ! shape function for all nodes evaluated at zeta (N is a row vector)
    PetscReal, pointer :: DN(:,:)     ! derivatives of shape function with respect to zeta
    PetscReal, pointer :: coord(:,:)  ! local coordinates of the vertices in the reference element
  end type shapefunction_type

  public :: ShapeFunctionInitialize, ShapeFunctionCalculate, &
            ShapeFunctionDestroy

  contains

! ************************************************************************** !

subroutine ShapeFunctionInitialize(shapefunction,option)
  !
  ! Allocate memory for shapefunction type
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !
  ! Updated by Piyoosh Jaysaval, PNNL
  ! Date: 2/16/2024
  !

  type(shapefunction_type) :: shapefunction
  type(option_type), optional, intent(inout) :: option
  PetscReal, pointer :: coord(:,:)

  select case(shapefunction%element_type)
    case(LINE_TYPE)
      allocate(shapefunction%N(TWO_INTEGER))
      allocate(shapefunction%DN(TWO_INTEGER,ONE_INTEGER))
      allocate(shapefunction%zeta(ONE_INTEGER))
      allocate(shapefunction%coord(TWO_INTEGER,ONE_INTEGER))
    case(TRI_TYPE)
      allocate(shapefunction%N(THREE_INTEGER))
      allocate(shapefunction%DN(THREE_INTEGER,TWO_INTEGER))
      allocate(shapefunction%zeta(TWO_INTEGER))
      allocate(shapefunction%coord(THREE_INTEGER,TWO_INTEGER))
    case(QUAD_TYPE)
      allocate(shapefunction%N(FOUR_INTEGER))
      allocate(shapefunction%DN(FOUR_INTEGER,TWO_INTEGER))
      allocate(shapefunction%zeta(TWO_INTEGER))
      allocate(shapefunction%coord(FOUR_INTEGER,TWO_INTEGER))
    case(WEDGE_TYPE)
      allocate(shapefunction%N(SIX_INTEGER))
      allocate(shapefunction%DN(SIX_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))
      allocate(shapefunction%coord(SIX_INTEGER,THREE_INTEGER))
    case(TET_TYPE)
      allocate(shapefunction%N(FOUR_INTEGER))
      allocate(shapefunction%DN(FOUR_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))
      allocate(shapefunction%coord(FOUR_INTEGER,THREE_INTEGER))
    case(PYR_TYPE)
      allocate(shapefunction%N(FIVE_INTEGER))
      allocate(shapefunction%DN(FIVE_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))
      allocate(shapefunction%coord(FIVE_INTEGER,THREE_INTEGER))
    case(HEX_TYPE)
      allocate(shapefunction%N(EIGHT_INTEGER))
      allocate(shapefunction%DN(EIGHT_INTEGER,THREE_INTEGER))
      allocate(shapefunction%zeta(THREE_INTEGER))
      allocate(shapefunction%coord(EIGHT_INTEGER,THREE_INTEGER))
    case default
      call ShapeFunctionError('Invalid element_type in ShapeFunctionInitialize.',option)
  end select

  shapefunction%zeta = 0.d0
  shapefunction%N  = 0.d0
  shapefunction%DN = 0.d0
  shapefunction%coord = 0.d0

  coord => shapefunction%coord

  select case (shapefunction%element_type)
    case (LINE_TYPE)
      coord(1,1) = -1.d0
      coord(2,1) =  1.d0

    case (TRI_TYPE)
      coord(1,:) = (/ 0.d0, 0.d0 /)
      coord(2,:) = (/ 1.d0, 0.d0 /)
      coord(3,:) = (/ 0.d0, 1.d0 /)

    case (QUAD_TYPE)
      coord(1,:) = (/ -1.d0, -1.d0 /)
      coord(2,:) = (/  1.d0, -1.d0 /)
      coord(3,:) = (/  1.d0,  1.d0 /)
      coord(4,:) = (/ -1.d0,  1.d0 /)

    case (WEDGE_TYPE)
      coord(1,:) = (/ 0.d0, 0.d0, -1.d0 /)
      coord(2,:) = (/ 1.d0, 0.d0, -1.d0 /)
      coord(3,:) = (/ 0.d0, 1.d0, -1.d0 /)
      coord(4,:) = (/ 0.d0, 0.d0,  1.d0 /)
      coord(5,:) = (/ 1.d0, 0.d0,  1.d0 /)
      coord(6,:) = (/ 0.d0, 1.d0,  1.d0 /)

    case (TET_TYPE)
      coord(1,:) = (/ 0.d0, 0.d0, 0.d0 /)
      coord(2,:) = (/ 1.d0, 0.d0, 0.d0 /)
      coord(3,:) = (/ 0.d0, 1.d0, 0.d0 /)
      coord(4,:) = (/ 0.d0, 0.d0, 1.d0 /)

    case (PYR_TYPE)
      coord(1,:) = (/ -1.d0, -1.d0, 0.d0 /)
      coord(2,:) = (/  1.d0, -1.d0, 0.d0 /)
      coord(3,:) = (/  1.d0,  1.d0, 0.d0 /)
      coord(4,:) = (/ -1.d0,  1.d0, 0.d0 /)
      coord(5,:) = (/  0.d0,  0.d0, 1.d0 /)

    case (HEX_TYPE)
      coord(1,:) = (/ -1.d0, -1.d0, -1.d0 /)
      coord(2,:) = (/  1.d0, -1.d0, -1.d0 /)
      coord(3,:) = (/  1.d0,  1.d0, -1.d0 /)
      coord(4,:) = (/ -1.d0,  1.d0, -1.d0 /)
      coord(5,:) = (/ -1.d0, -1.d0,  1.d0 /)
      coord(6,:) = (/  1.d0, -1.d0,  1.d0 /)
      coord(7,:) = (/  1.d0,  1.d0,  1.d0 /)
      coord(8,:) = (/ -1.d0,  1.d0,  1.d0 /)

    case default
      call ShapeFunctionError('Invalid element_type in ShapeFunctionInitialize.',option)
  end select

end subroutine ShapeFunctionInitialize

! ************************************************************************** !

subroutine ShapeFunctionCalculate(shapefunction,option)
  !
  ! Subroutine provides shape functions and its derivatives
  ! at a given spatial point (in the reference element) 'zeta' for various finite
  ! elements.
  ! Input variables
  ! element_type: element type
  ! L2: one-dimensional element -1 <= x <= +1
  ! Q4: four node quadrilateral element -1 <= x, y <= +1
  ! Q9: nine node quadrilateral element
  ! T3: three node right-angled triangle with h = b = 1
  ! T6: six node right-angled triangle
  ! B8: eight node right-angled tetrahedron element
  ! zeta: coordinates of a point in the reference element
  ! Output variables
  ! N: shape functions for all nodes evaluated at zeta (N is a row
  ! vector!)
  ! DN: derivatives of shape functions with respect to zeta
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !
  ! Updated by Piyoosh Jaysaval, PNNL
  ! Date: 2/16/2024
  !
!

  type(shapefunction_type) :: shapefunction
  type(option_type), optional, intent(inout) :: option
  PetscReal, pointer :: zeta(:)
  PetscReal, pointer :: N(:)
  PetscReal, pointer :: DN(:,:)
  PetscInt :: i
  PetscReal :: c4, c8
  PetscReal :: x1m, x1p, x2m, x2p, x3m, x3p
  PetscBool :: zeta_valid

  N => shapefunction%N
  DN => shapefunction%DN
  zeta => shapefunction%zeta

  ! Element-type-aware coordinate range validation
  zeta_valid = PETSC_TRUE
  select case(shapefunction%element_type)
    case(TRI_TYPE)
      ! Barycentric: each zeta(i) in [0,1], sum(zeta) <= 1
      do i = 1, size(zeta)
        if (zeta(i) < 0.d0 .or. zeta(i) > 1.d0) zeta_valid = PETSC_FALSE
      enddo
      if (zeta(1) + zeta(2) > 1.d0 + 1.d-12) zeta_valid = PETSC_FALSE
    case(TET_TYPE)
      ! Barycentric: each zeta(i) in [0,1], sum(zeta) <= 1
      do i = 1, size(zeta)
        if (zeta(i) < 0.d0 .or. zeta(i) > 1.d0) zeta_valid = PETSC_FALSE
      enddo
      if (zeta(1) + zeta(2) + zeta(3) > 1.d0 + 1.d-12) zeta_valid = PETSC_FALSE
    case(WEDGE_TYPE)
      ! Triangular face: zeta(1),zeta(2) in [0,1] with sum <= 1
      ! Axial: zeta(3) in [-1,1]
      if (zeta(1) < 0.d0 .or. zeta(1) > 1.d0) zeta_valid = PETSC_FALSE
      if (zeta(2) < 0.d0 .or. zeta(2) > 1.d0) zeta_valid = PETSC_FALSE
      if (zeta(1) + zeta(2) > 1.d0 + 1.d-12) zeta_valid = PETSC_FALSE
      if (zeta(3) < -1.d0 .or. zeta(3) > 1.d0) zeta_valid = PETSC_FALSE
    case(PYR_TYPE)
      ! Base: zeta(1),zeta(2) in [-1,1]; apex: zeta(3) in [0,1]
      if (zeta(1) < -1.d0 .or. zeta(1) > 1.d0) zeta_valid = PETSC_FALSE
      if (zeta(2) < -1.d0 .or. zeta(2) > 1.d0) zeta_valid = PETSC_FALSE
      if (zeta(3) < 0.d0  .or. zeta(3) > 1.d0) zeta_valid = PETSC_FALSE
    case default
      ! LINE_TYPE, QUAD_TYPE, HEX_TYPE: all coords in [-1,1]
      do i = 1, size(zeta)
        if (zeta(i) < -1.d0 .or. zeta(i) > 1.d0) zeta_valid = PETSC_FALSE
      enddo
  end select
  if (zeta_valid .eqv. PETSC_FALSE) then
    call ShapeFunctionError('Zeta coordinates out of valid range ' // &
      'for the element type.',option)
    return
  endif

  select case(shapefunction%element_type)
    case(LINE_TYPE)
      N(1) = 0.5d0 * (1.d0 - zeta(1))
      N(2) = 0.5d0 * (1.d0 + zeta(1))
      DN(1,1) = -0.5d0
      DN(2,1) =  0.5d0

    case(TRI_TYPE)
      N(1) = 1.d0 - zeta(1) - zeta(2)
      N(2) = zeta(1)
      N(3) = zeta(2)
      DN(1,1) = -1.d0
      DN(1,2) = -1.d0
      DN(2,1) =  1.d0
      DN(2,2) =  0.d0
      DN(3,1) =  0.d0
      DN(3,2) =  1.d0

    case(QUAD_TYPE)
      c4 = 0.25d0
      x1m = 1.d0 - zeta(1)
      x1p = 1.d0 + zeta(1)
      x2m = 1.d0 - zeta(2)
      x2p = 1.d0 + zeta(2)

      N(1) = c4 * x1m * x2m
      N(2) = c4 * x1p * x2m
      N(3) = c4 * x1p * x2p
      N(4) = c4 * x1m * x2p

      DN(1,1) = -c4 * x2m
      DN(1,2) = -c4 * x1m
      DN(2,1) =  c4 * x2m
      DN(2,2) = -c4 * x1p
      DN(3,1) =  c4 * x2p
      DN(3,2) =  c4 * x1p
      DN(4,1) = -c4 * x2p
      DN(4,2) =  c4 * x1m

    case(WEDGE_TYPE)
      x1m = 1.d0 - zeta(1) - zeta(2)
      x3m = 1.d0 - zeta(3)
      x3p = 1.d0 + zeta(3)

      ! Bottom triangle nodes at z = -1
      N(1) = 0.5d0 * x1m     * x3m
      N(2) = 0.5d0 * zeta(1) * x3m
      N(3) = 0.5d0 * zeta(2) * x3m

      ! Top triangle nodes at z = +1
      N(4) = 0.5d0 * x1m     * x3p
      N(5) = 0.5d0 * zeta(1) * x3p
      N(6) = 0.5d0 * zeta(2) * x3p

      DN(1,1) = -0.5d0 * x3m
      DN(1,2) = -0.5d0 * x3m
      DN(1,3) = -0.5d0 * x1m

      DN(2,1) =  0.5d0 * x3m
      DN(2,2) =  0.d0
      DN(2,3) = -0.5d0 * zeta(1)

      DN(3,1) =  0.d0
      DN(3,2) =  0.5d0 * x3m
      DN(3,3) = -0.5d0 * zeta(2)

      DN(4,1) = -0.5d0 * x3p
      DN(4,2) = -0.5d0 * x3p
      DN(4,3) =  0.5d0 * x1m

      DN(5,1) =  0.5d0 * x3p
      DN(5,2) =  0.d0
      DN(5,3) =  0.5d0 * zeta(1)

      DN(6,1) =  0.d0
      DN(6,2) =  0.5d0 * x3p
      DN(6,3) =  0.5d0 * zeta(2)

    case(TET_TYPE)
      N(1) = 1.d0 - zeta(1) - zeta(2) - zeta(3)
      N(2) = zeta(1)
      N(3) = zeta(2)
      N(4) = zeta(3)

      DN(1,1) = -1.d0
      DN(1,2) = -1.d0
      DN(1,3) = -1.d0
      DN(2,1) =  1.d0
      DN(2,2) =  0.d0
      DN(2,3) =  0.d0
      DN(3,1) =  0.d0
      DN(3,2) =  1.d0
      DN(3,3) =  0.d0
      DN(4,1) =  0.d0
      DN(4,2) =  0.d0
      DN(4,3) =  1.d0

    case(PYR_TYPE)
      c4 = 0.25d0
      x1m = 1.d0 - zeta(1)
      x1p = 1.d0 + zeta(1)
      x2m = 1.d0 - zeta(2)
      x2p = 1.d0 + zeta(2)
      x3m = 1.d0 - zeta(3)

      N(1) = c4 * x1m * x2m * x3m
      N(2) = c4 * x1p * x2m * x3m
      N(3) = c4 * x1p * x2p * x3m
      N(4) = c4 * x1m * x2p * x3m
      N(5) = zeta(3)

      DN(1,1) = -c4 * x2m * x3m
      DN(1,2) = -c4 * x1m * x3m
      DN(1,3) = -c4 * x1m * x2m

      DN(2,1) =  c4 * x2m * x3m
      DN(2,2) = -c4 * x1p * x3m
      DN(2,3) = -c4 * x1p * x2m

      DN(3,1) =  c4 * x2p * x3m
      DN(3,2) =  c4 * x1p * x3m
      DN(3,3) = -c4 * x1p * x2p

      DN(4,1) = -c4 * x2p * x3m
      DN(4,2) =  c4 * x1m * x3m
      DN(4,3) = -c4 * x1m * x2p

      DN(5,1) = 0.d0
      DN(5,2) = 0.d0
      DN(5,3) = 1.d0

    case(HEX_TYPE)
      c8 = 0.125d0
      x1m = 1.d0 - zeta(1)
      x1p = 1.d0 + zeta(1)
      x2m = 1.d0 - zeta(2)
      x2p = 1.d0 + zeta(2)
      x3m = 1.d0 - zeta(3)
      x3p = 1.d0 + zeta(3)

      N(1) = c8 * x1m * x2m * x3m
      N(2) = c8 * x1p * x2m * x3m
      N(3) = c8 * x1p * x2p * x3m
      N(4) = c8 * x1m * x2p * x3m
      N(5) = c8 * x1m * x2m * x3p
      N(6) = c8 * x1p * x2m * x3p
      N(7) = c8 * x1p * x2p * x3p
      N(8) = c8 * x1m * x2p * x3p

      DN(1,1) = -c8 * x2m * x3m
      DN(1,2) = -c8 * x1m * x3m
      DN(1,3) = -c8 * x1m * x2m

      DN(2,1) =  c8 * x2m * x3m
      DN(2,2) = -c8 * x1p * x3m
      DN(2,3) = -c8 * x1p * x2m

      DN(3,1) =  c8 * x2p * x3m
      DN(3,2) =  c8 * x1p * x3m
      DN(3,3) = -c8 * x1p * x2p

      DN(4,1) = -c8 * x2p * x3m
      DN(4,2) =  c8 * x1m * x3m
      DN(4,3) = -c8 * x1m * x2p

      DN(5,1) = -c8 * x2m * x3p
      DN(5,2) = -c8 * x1m * x3p
      DN(5,3) =  c8 * x1m * x2m

      DN(6,1) =  c8 * x2m * x3p
      DN(6,2) = -c8 * x1p * x3p
      DN(6,3) =  c8 * x1p * x2m

      DN(7,1) =  c8 * x2p * x3p
      DN(7,2) =  c8 * x1p * x3p
      DN(7,3) =  c8 * x1p * x2p

      DN(8,1) = -c8 * x2p * x3p
      DN(8,2) =  c8 * x1m * x3p
      DN(8,3) =  c8 * x1m * x2p

    case default
      call ShapeFunctionError('Invalid element_type in ShapeFunctionCalculate.',option)
  end select

end subroutine ShapeFunctionCalculate

! ************************************************************************** !

subroutine ShapeFunctionDestroy(shapefunction)
  !
  ! Deallocate memory for shapefunction type
  !
  ! Author: Satish Karra, LANL
  ! Date: 5/17/2013
  !

  type(shapefunction_type) :: shapefunction

  deallocate(shapefunction%N)
  nullify(shapefunction%N)
  deallocate(shapefunction%DN)
  nullify(shapefunction%DN)
  deallocate(shapefunction%zeta)
  nullify(shapefunction%zeta)
  deallocate(shapefunction%coord)
  nullify(shapefunction%coord)

end subroutine ShapeFunctionDestroy

! ************************************************************************** !

subroutine ShapeFunctionError(message,option)

  implicit none

  character(len=*), intent(in) :: message
  type(option_type), optional, intent(inout) :: option

  if (present(option)) then
    option%io_buffer = trim(message)
    call PrintErrMsg(option)
  else
    print *, trim(message)
    stop
  endif

end subroutine ShapeFunctionError

end module Shape_Function_module
