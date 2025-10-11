module Debug_module

#include "petsc/finclude/petscmat.h"
  use petscmat
  use PFLOTRAN_Constants_module
  use Petsc_Utility_module

  implicit none

  private

  type, public :: debug_type
    PetscBool :: vecview_residual
    PetscBool :: vecview_solution
    PetscBool :: matview_Matrix
    PetscBool :: matview_Matrix_detailed
    PetscBool :: norm_Matrix

    PetscInt  :: output_format
    PetscBool :: verbose_filename

    PetscBool :: print_couplers
    PetscBool :: print_regions_hdf5
    PetscBool :: print_regions_tec
    character(len=MAXSTRINGLENGTH) :: coupler_string
    PetscBool :: print_waypoints
  end type debug_type

  interface DebugMatView
    module procedure DebugMatView1
    module procedure DebugMatView2
  end interface

  interface DebugVecView
    module procedure DebugVecView1
    module procedure DebugVecView2
  end interface

  public :: DebugCreate, &
            DebugRead, &
            DebugWriteFilename, &
            DebugMatView, &
            DebugVecView, &
            DebugDestroy

contains

! ************************************************************************** !

function DebugCreate()
  !
  ! Create object that stores debugging options for PFLOW
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !

  implicit none

  type(debug_type), pointer :: DebugCreate

  type(debug_type), pointer :: debug

  allocate(debug)

  debug%vecview_residual = PETSC_FALSE
  debug%vecview_solution = PETSC_FALSE
  debug%matview_Matrix = PETSC_FALSE
  debug%matview_Matrix_detailed = PETSC_FALSE
  debug%norm_Matrix = PETSC_FALSE

  debug%output_format = VIEWER_ASCII_FORMAT
  debug%verbose_filename = PETSC_FALSE

  debug%print_couplers = PETSC_FALSE
  debug%print_regions_hdf5 = PETSC_FALSE
  debug%print_regions_tec = PETSC_FALSE
  debug%coupler_string = ''
  debug%print_waypoints = PETSC_FALSE

  DebugCreate => debug

end function DebugCreate

! ************************************************************************** !

subroutine DebugRead(debug,input,option)
  !
  ! Reads debugging data from the input file
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(debug_type) :: debug
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: keyword

  input%ierr = INPUT_ERROR_NONE
  call InputPushBlock(input,option)
  do

    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,keyword)
    call InputErrorMsg(input,option,'keyword','DEBUG')
    call StringToUpper(keyword)

    select case(trim(keyword))

      case('PRINT_SOLUTION','VECVIEW_SOLUTION','VIEW_SOLUTION')
        debug%vecview_solution = PETSC_TRUE
      case('PRINT_RESIDUAL','VECVIEW_RESIDUAL','VIEW_RESIDUAL')
        debug%vecview_residual = PETSC_TRUE
      case('PRINT_JACOBIAN','matview_Matrix','VIEW_JACOBIAN')
        debug%matview_Matrix = PETSC_TRUE
      case('PRINT_JACOBIAN_NORM','norm_Matrix')
        debug%norm_Matrix = PETSC_TRUE
      case('PRINT_MATRIX','MATVIEW_MATRIX','VIEW_MATRIX')
        debug%matview_Matrix = PETSC_TRUE
      case('PRINT_REGIONS')
        debug%print_regions_tec = PETSC_TRUE
        debug%print_regions_hdf5 = PETSC_TRUE
      case('PRINT_REGIONS_TECPLOT')
        debug%print_regions_tec = PETSC_TRUE
      case('PRINT_REGIONS_HDF5')
        debug%print_regions_hdf5 = PETSC_TRUE
      case('PRINT_COUPLERS','PRINT_COUPLER')
        debug%print_couplers = PETSC_TRUE
        debug%coupler_string = trim(adjustl(input%buf))
      case('PRINT_JACOBIAN_DETAILED','matview_Matrix_DETAILED', &
           'VIEW_JACOBIAN_DETAILED')
        debug%matview_Matrix_detailed = PETSC_TRUE
      case('PRINT_WAYPOINTS')
        debug%print_waypoints = PETSC_TRUE
      case('APPEND_COUNTS_TO_FILENAME','APPEND_COUNTS_TO_FILENAMES')
        debug%verbose_filename = PETSC_TRUE
      case('FORMAT')
        call InputReadCard(input,option,keyword)
        call InputErrorMsg(input,option,'keyword','DEBUG,FORMAT')
        call StringToUpper(keyword)
        select case(keyword)
          case('ASCII')
            debug%output_format = VIEWER_ASCII_FORMAT
          case('BINARY')
            debug%output_format = VIEWER_BINARY_FORMAT
          case('MATLAB')
            debug%output_format = VIEWER_MATLAB_FORMAT
          case('NATIVE','PARALLEL')
            debug%output_format = VIEWER_NATIVE_FORMAT
        end select
      case default
        call InputKeywordUnrecognized(input,keyword,'DEBUG',option)
    end select

  enddo
  call InputPopBlock(input,option)

end subroutine DebugRead

! ************************************************************************** !

subroutine DebugWriteFilename(debug,filename,prefix,suffix,ts,ts_cut,ni)
  !
  ! Appends timestep, timestep cut, and Newton iteration counts to a
  ! filename.
  !
  ! Author: Glenn Hammond
  ! Date: 10/23/18
  !

  implicit none

  type(debug_type) :: debug
  character(len=*) :: filename
  character(len=*) :: prefix
  character(len=*) :: suffix
  PetscInt :: ts
  PetscInt :: ts_cut
  PetscInt :: ni

  character(len=MAXWORDLENGTH) :: word

  filename = adjustl(prefix)
  if (debug%verbose_filename) then
    write(word,*) ts
    filename = trim(filename) // '_ts' // adjustl(word)
    write(word,*) ts_cut
    filename = trim(filename) // '_tc' // adjustl(word)
    write(word,*) ni
    filename = trim(filename) // '_ni' // adjustl(word)
  endif
  if (len_trim(suffix) > 0) then
    filename = trim(filename) // '.' // adjustl(suffix)
  endif

end subroutine DebugWriteFilename

! ************************************************************************** !

subroutine DebugMatView1(debug,A,prefix,suffix,ts,ts_cut,ni,option)
  !
  ! Dumps a PETSc Mat through a viewer object
  !
  ! Author: Glenn Hammond
  ! Date: 10/10/25
  !
  use Option_module

  implicit none

  type(debug_type) :: debug
  Mat :: A
  character(len=*) :: prefix
  character(len=*) :: suffix
  PetscInt :: ts
  PetscInt :: ts_cut
  PetscInt :: ni
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  call DebugWriteFilename(debug,string,prefix,suffix,ts,ts_cut,ni)
  call PUMatView(A,string,debug%output_format,option)

end subroutine DebugMatView1

! ************************************************************************** !

subroutine DebugMatView2(debug,A,filename,option)
  !
  ! Dumps a PETSc Mat through a viewer object
  !
  ! Author: Glenn Hammond
  ! Date: 10/10/25
  !
  use Option_module

  implicit none

  type(debug_type) :: debug
  Mat :: A
  character(len=*) :: filename
  type(option_type) :: option

  call PUMatView(A,filename,debug%output_format,option)

end subroutine DebugMatView2

! ************************************************************************** !

subroutine DebugVecView1(debug,v,prefix,suffix,ts,ts_cut,ni,option)
  !
  ! Dumps a PETSc Vec through a viewer object
  !
  ! Author: Glenn Hammond
  ! Date: 10/10/25
  !
  use Option_module

  implicit none

  type(debug_type) :: debug
  Vec :: v
  character(len=*) :: prefix
  character(len=*) :: suffix
  PetscInt :: ts
  PetscInt :: ts_cut
  PetscInt :: ni
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  call DebugWriteFilename(debug,string,prefix,suffix,ts,ts_cut,ni)
  call PUVecView(v,string,debug%output_format,option)

end subroutine DebugVecView1

! ************************************************************************** !

subroutine DebugVecView2(debug,v,filename,option)
  !
  ! Dumps a PETSc Vec through a viewer object
  !
  ! Author: Glenn Hammond
  ! Date: 10/10/25
  !
  use Option_module

  implicit none

  type(debug_type) :: debug
  Vec :: v
  character(len=*) :: filename
  type(option_type) :: option

  call PUVecView(v,filename,debug%output_format,option)

end subroutine DebugVecView2

! ************************************************************************** !

subroutine DebugDestroy(debug)
  !
  ! Deallocates memory associated with debug object
  !
  ! Author: Glenn Hammond
  ! Date: 12/21/07
  !
  implicit none

  type(debug_type), pointer :: debug

  if (.not.associated(debug)) return

  deallocate(debug)
  nullify(debug)

end subroutine DebugDestroy

end module Debug_module
