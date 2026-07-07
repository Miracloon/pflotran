module Anisotropy_Geom_Data_module

! Module for managing additional data needed in flux calculations to
! deal with full (off-diagonal) anisotropy terms
!
! Author: Steve Benbow
! Date: 01/06/20

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  implicit none

  ! --
  ! -- Data structures for calculating anisotropic gradient contributions
  ! --

  !  When calculating the flux across interface marked '---->' in a structured grid setting, if data
  !  is only used from cells 'up' and 'dn' then this only allows fluxes based on, e.g., permeability
  !  components in the direction of the flux to be calculated.  Anisotropic permeablity components will
  !  be ignored.
  !  To include anisotropic components, the gradient in the plane of interface must
  !  also be computed.  In a 2-D structured grid, this can be done by approximating the
  !  gradient using the up/downstream adjacent cells marked A,B.
  !
  !  .......................
  !  :          :          :
  !  :   up_B   :   dn_B   :   <- y-adjacent up/dn cells (+ve dir)
  !  :          :          :
  !  -----------------------
  !  |          |          |
  !  |    up  --|--> dn    |
  !  |          |          |
  !  -----------------------
  !  :          :          :
  !  :   up_A   :   dn_A   :   <- y-adjacent up/dn cells (-ve dir)
  !  :          :          :
  !  .......................
  !
  !
  !  The approach implemented here a general approach for unstructured grids that can
  !  equally be applied in a structured grid setting.  It is described in
  !    Benbow et al. (2026) "PFLOTRAN – Updates for Fully Anisotropic Permeability".

  ! --
  ! -- Data structures for calculating anisotropic gradient contributions
  ! --

  ! represents adjacency data for one of the cells (up or dn) neighbouring an interface
  type, public :: aniso_cell_data_type
    PetscInt :: cell_id        ! id of cell (up or dn)
    PetscReal :: cell_pos(3)   ! (x,y,z coord of cell centre)
    ! g.c, where c is cell centre (this are differenced to get g.d, where d is distance between points)
    PetscReal :: loc_gravity
    ! adjacency info
    PetscInt :: num_adj ! number of adjacent cells (<=4)
    PetscInt, pointer :: cell_id_adj(:) ! ids of adjacent cells that share an edge with the connection face. -ve if ghost cell
    PetscInt, pointer :: face_id_adj(:) ! ids of faces of cell that share an edge with the connection face.  -ve if ghost-ghost face
    PetscReal, pointer :: loc_gravity_adj(:)

    ! (u,v,w) is coordinate system orthogonal to the common face.  u is normal to the face and (v,w) is in the plane of the face.
    ! Interpolation machinery - 2x2 matrix inverse is required.
    PetscReal, pointer :: dists_u(:)    ! num_adj
    PetscReal, pointer :: dists_vw(:,:) ! 2 x num_adj
    ! LSMatInv is the inverse of the Gram matrix in the normal equations solution for the (v,w) components of grad(f)
    ! where f is the scalar of interest (e.g. h for Richards' equation)
    PetscReal :: LSMatInv(2,2)

    contains
      procedure :: InitialiseAdj
      procedure :: DestroyAdj
  end type aniso_cell_data_type

  type, public :: aniso_geom_data_type
    ! not strictly necessary here - could be inferred from parent connection_set, but convenient to pass with single ptr
    PetscInt :: conn_face_id      ! id of connection face

    type(aniso_cell_data_type), pointer :: up => null()
    type(aniso_cell_data_type), pointer :: dn => null()

    ! (u,v,w) is a coord system relative to face.
    PetscReal, pointer :: uvec(:)      ! u-direction unit vector (normal to face, +ve dnstream) - in (x,y,z) coords
    PetscReal, pointer :: vvec(:)      ! v-direction unit vector (in plane of face) - in (x,y,z) coords
    PetscReal, pointer :: wvec(:)      ! w-direction unit vector (in plane of face) - in (x,y,z) coords

    PetscBool :: is_initialised = PETSC_FALSE

    contains
      procedure :: InitialiseMembers
      procedure :: DestroyMembers
      procedure :: MakeInterpolationParams
  end type aniso_geom_data_type

  ! --
  ! -- Methods
  ! --

  CONTAINS

    subroutine InitialiseMembers(self,num_adj)
      class(aniso_geom_data_type) :: self
      PetscInt, intent(in) :: num_adj

      allocate(self%uvec(1:3))
      allocate(self%vvec(1:3))
      allocate(self%wvec(1:3))
      self%uvec = 0.d0
      self%vvec = 0.d0
      self%wvec = 0.d0

      allocate(self%up)
      allocate(self%dn)
      call self%up%InitialiseAdj(num_adj)
      call self%dn%InitialiseAdj(num_adj)

      self%is_initialised = PETSC_TRUE
    end subroutine InitialiseMembers

    subroutine DestroyMembers(self)
      class(aniso_geom_data_type) :: self

      if( self%is_initialised) then
        deallocate(self%uvec)
        deallocate(self%vvec)
        deallocate(self%wvec)

        call self%up%DestroyAdj()
        call self%dn%DestroyAdj()

        deallocate(self%up)
        deallocate(self%dn)

        self%is_initialised = PETSC_FALSE
      endif
    end subroutine DestroyMembers

    subroutine InitialiseAdj(self,num_adj)
      class(aniso_cell_data_type) :: self
      PetscInt, intent(in) :: num_adj

      allocate(self%cell_id_adj(num_adj))
      allocate(self%face_id_adj(num_adj))
      allocate(self%loc_gravity_adj(num_adj))

      allocate(self%dists_u(1:num_adj))
      allocate(self%dists_vw(1:2,1:num_adj))

      self%cell_id_adj = 0
      self%face_id_adj = 0
      self%loc_gravity_adj = 0
      self%num_adj = 0
    end subroutine InitialiseAdj

    subroutine DestroyAdj(self)
      class(aniso_cell_data_type) :: self

      deallocate(self%cell_id_adj)
      deallocate(self%face_id_adj)
      deallocate(self%loc_gravity_adj)

      deallocate(self%dists_u)
      deallocate(self%dists_vw)
    end subroutine DestroyAdj

    subroutine MakeInterpolationParams(self, xc, yc, zc )

      use Geometry_module
      use Utility_module, only : DotProduct

      class(aniso_geom_data_type) :: self
      PetscReal, intent(in) :: xc(*), yc(*), zc(*)

      PetscInt :: iupdn, iadj
      type(aniso_cell_data_type), pointer :: aniso_updn
      type(point3d_type) :: point1, point2
      PetscReal :: v1(3)
      PetscReal :: dot, det
      PetscReal :: MLS(2,2) ! Gram matrix in the normal equations solution for the (v,w) components of grad(f)
      PetscBool :: no_v_contrib, no_w_contrib

      ! Calculate distances up/dn->adj in (u,v,w) coords
      do iupdn = 1,2
        if (iupdn==1) then
          ! upstream
          aniso_updn => self%up
        else
          ! downstream
          aniso_updn => self%dn
        endif

        aniso_updn%dists_u = 0.d0 ! )
        aniso_updn%dists_vw = 0.d0 ! )

        point1%x = xc(aniso_updn%cell_id)
        point1%y = yc(aniso_updn%cell_id)
        point1%z = zc(aniso_updn%cell_id)

        no_v_contrib = PETSC_TRUE ! reset below if neighbour in v direction is found
        no_w_contrib = PETSC_TRUE ! reset below if neighbour in w direction is found

        ! loop over adjacent cells that will contribute to anisotropy, calculate (v,w) dist
        do iadj = 1,aniso_updn%num_adj
          ! adj cells have -ve cell id if they are ghost cells, hence abs(..) here
          point2%x = xc(abs(aniso_updn%cell_id_adj(iadj)))
          point2%y = yc(abs(aniso_updn%cell_id_adj(iadj)))
          point2%z = zc(abs(aniso_updn%cell_id_adj(iadj)))

          ! vec from up/dn cell centre to adj cell centre
          v1(1) = point2%x - point1%x
          v1(2) = point2%y - point1%y
          v1(3) = point2%z - point1%z

          ! distances (vector components) in (u,v,w) coords
          aniso_updn%dists_u(iadj)    = DotProduct(v1,self%uvec)
          aniso_updn%dists_vw(1,iadj) = DotProduct(v1,self%vvec)
          aniso_updn%dists_vw(2,iadj) = DotProduct(v1,self%wvec)

          if(aniso_updn%dists_vw(1,iadj) .ne. 0.d0) no_v_contrib = PETSC_FALSE
          if(aniso_updn%dists_vw(2,iadj) .ne. 0.d0) no_w_contrib = PETSC_FALSE
        enddo

        if(aniso_updn%num_adj>1) then
          ! compute least-squares mtx (sum of outer product of dists_vw vectors)
          MLS = 0.d0 ! Initialise 2x2 least-squares mtx
          do iadj = 1,aniso_updn%num_adj
            dot = aniso_updn%dists_vw(1,iadj)*aniso_updn%dists_vw(1,iadj) &
                  + aniso_updn%dists_vw(2,iadj)*aniso_updn%dists_vw(2,iadj)
            MLS(1,1) = MLS(1,1) + aniso_updn%dists_vw(1,iadj)*aniso_updn%dists_vw(1,iadj)/dot
            MLS(1,2) = MLS(1,2) + aniso_updn%dists_vw(1,iadj)*aniso_updn%dists_vw(2,iadj)/dot
            MLS(2,1) = MLS(2,1) + aniso_updn%dists_vw(2,iadj)*aniso_updn%dists_vw(1,iadj)/dot
            MLS(2,2) = MLS(2,2) + aniso_updn%dists_vw(2,iadj)*aniso_updn%dists_vw(2,iadj)/dot
          enddo

          ! invert least-squares mtx (handling singular cases - ok for LS solution)
          if(no_v_contrib) then      ! MLS only has non-zero (2,2) entry
            aniso_updn%LSMatInv = 0
            aniso_updn%LSMatInv(2,2) = 1.d0 / MLS(2,2)
          elseif(no_w_contrib) then  ! MLS only has non-zero (1,1) entry
            aniso_updn%LSMatInv = 0
            aniso_updn%LSMatInv(1,1) = 1.d0 / MLS(1,1)
          else
            ! 2x2 case
            det = MLS(1,1)*MLS(2,2) - MLS(1,2)*MLS(2,1)
            if(det .eq. 0.d0) then
              ! singular case - cell->neighbour vectors are parallel.
              ! introduce a small perturbation
              ! (could improve - e.g. calculate grad in cell->neighbour direction and assume zero orthogonal to this)
              MLS(1,1) = MLS(1,1) * 1.001d0
              det = MLS(1,1)*MLS(2,2) - MLS(1,2)*MLS(2,1)
            endif
            aniso_updn%LSMatInv(1,1) =  MLS(2,2) / det
            aniso_updn%LSMatInv(1,2) = -MLS(1,2) / det
            aniso_updn%LSMatInv(2,1) = -MLS(2,1) / det
            aniso_updn%LSMatInv(2,2) =  MLS(1,1) / det
          endif
        else
          ! under-determined.
          ! Taking MLS=I will provide an estimate of the gradient in the same ratio as components of dists_vw.
          aniso_updn%LSMatInv = 0
          aniso_updn%LSMatInv(1,1) = 1.d0
          aniso_updn%LSMatInv(2,2) = 1.d0
        endif
      enddo

    end subroutine MakeInterpolationParams

end module Anisotropy_Geom_Data_module
