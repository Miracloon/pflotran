module Richards_Aniso_module

  ! Module for managing additional data needed in Richards' flux calculations to
  ! deal with full (off-diagonal) anisotropy terms
  !
  ! The approach implemented here is described in
  !    Benbow et al. (2026) "PFLOTRAN – Updates for Fully Anisotropic Permeability".
  !
  ! Author: Steve Benbow
  ! Date: 01/06/20

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module

  use Anisotropy_Geom_Data_module


  implicit none

  public :: PermInnerProd

  ! --
  ! -- Data structures for calculating anisotropic gradient contributions
  ! --

  type, public :: aniso_richards_cell_type
    ! Permeability tensor w.r.t. (u,v,w) coords in cells - working data, overwritten throughout simulation (in case timedep)
    PetscReal :: perm_uu
    PetscReal :: perm_vv
    PetscReal :: perm_ww
    PetscReal :: perm_uv
    PetscReal :: perm_uw
    PetscReal :: perm_vw

    ! Adjacent pressure data - working data, overwritten throughout simulation
    PetscReal, pointer :: pres_adj(:) ! pressure adjacent to up/dn cell
  end type aniso_richards_cell_type

  type, public :: aniso_richards_data_type
    type(aniso_geom_data_type), pointer :: aniso_geom
    type(aniso_richards_cell_type), pointer :: up => null()
    type(aniso_richards_cell_type), pointer :: dn => null()

    PetscBool :: is_allocated = PETSC_FALSE

    ! derivative data
    PetscReal :: dv_darcy_aniso_dp_up
    PetscReal :: dv_darcy_aniso_dp_dn
    PetscReal, pointer :: dv_darcy_aniso_dp_up_adj(:)
    PetscReal, pointer :: dv_darcy_aniso_dp_dn_adj(:)
    PetscReal, pointer :: Jup_adj(:)
    PetscReal, pointer :: Jdn_adj(:)

    contains
      procedure, pass :: Initialise
      procedure, pass :: Destroy

  end type aniso_richards_data_type

  type, public :: aniso_richards_work_type
    ! structure for storing work data for anisotropic derivative calculations
    PetscReal :: d2phi_dv_up_dp_up
    PetscReal :: d2phi_dv_up_dp_dn
    PetscReal, pointer :: d2phi_dv_up_dp_up_adj(:)
    PetscReal :: d2phi_dw_up_dp_up
    PetscReal :: d2phi_dw_up_dp_dn
    PetscReal, pointer :: d2phi_dw_up_dp_up_adj(:)
    PetscReal :: d2phi_dv_dn_dp_dn
    PetscReal :: d2phi_dv_dn_dp_up
    PetscReal, pointer :: d2phi_dv_dn_dp_dn_adj(:)
    PetscReal :: d2phi_dw_dn_dp_dn
    PetscReal :: d2phi_dw_dn_dp_up
    PetscReal, pointer :: d2phi_dw_dn_dp_dn_adj(:)
  end type aniso_richards_work_type


  ! --
  ! -- Methods
  ! --


  CONTAINS
    function PermInnerProd( u, v, K_xx, K_yy, K_zz, K_xy, K_xz, K_yz )
      ! compute u^T K v, with K symmetric

      PetscReal :: u(3)
      PetscReal :: v(3)
      PetscReal :: K_xx,K_yy,K_zz,K_xy,K_xz,K_yz
      PetscReal :: PermInnerProd

      PermInnerProd =  u(1) * (K_xx*v(1) + K_xy*v(2) + K_xz*v(3)) &
                     + u(2) * (K_xy*v(1) + K_yy*v(2) + K_yz*v(3)) &
                     + u(3) * (K_xz*v(1) + K_yz*v(2) + K_zz*v(3))
    end function PermInnerProd

    subroutine Initialise(self)
      class(aniso_richards_data_type), intent(inout) :: self
      PetscInt :: sz_up
      PetscInt :: sz_dn

      sz_up = self%aniso_geom%up%num_adj
      sz_dn = self%aniso_geom%dn%num_adj

      allocate(self%dv_darcy_aniso_dp_up_adj(sz_up))
      allocate(self%dv_darcy_aniso_dp_dn_adj(sz_dn))
      allocate(self%Jup_adj(sz_up))
      allocate(self%Jdn_adj(sz_dn))

      self%dv_darcy_aniso_dp_up_adj = 0
      self%dv_darcy_aniso_dp_dn_adj = 0
      self%Jup_adj = 0
      self%Jdn_adj = 0

      ! allocate pressure adjacency data
      allocate(self%up)
      allocate(self%dn)
      allocate(self%up%pres_adj(sz_up))
      allocate(self%dn%pres_adj(sz_dn))

      self%up%pres_adj = 0.d0
      self%dn%pres_adj = 0.d0

      self%is_allocated = PETSC_TRUE
    end subroutine Initialise

    subroutine Destroy(self)
      class(aniso_richards_data_type), intent(inout) :: self

      if (self%is_allocated) then
        deallocate(self%dv_darcy_aniso_dp_up_adj)
        deallocate(self%dv_darcy_aniso_dp_dn_adj)

        deallocate(self%Jup_adj)
        deallocate(self%Jdn_adj)

        ! deallocate pressure adjacency data
        deallocate(self%up%pres_adj)
        deallocate(self%dn%pres_adj)

        deallocate(self%up)
        deallocate(self%dn)

        self%is_allocated = PETSC_FALSE
      endif
    end subroutine Destroy

    subroutine RichardsFluxAniso(aniso_rich, &
                                 material_auxvar_up,material_auxvar_dn, &
                                 global_auxvar_up,global_auxvar_dn, &
                                 option, &
                                 dist, &
                                 ukvr, &                     ! kvr = krel/visc
                                 dukvr_dp_up, dukvr_dp_dn, & ! derivs only needed if calculating derivatives, o/w set zero
                                 v_darcy_iso, &
                                 dv_darcy_iso_dp_up, dv_darcy_iso_dp_dn, & ! derivs only needed if calculating derivatives, o/w set zero
                                 calc_derivs, & ! controls whether derivatives are calculated
                                 v_darcy_aniso)
      !
      ! Analogue of RichardsFlux and RichardsFluxDerivative in richards_common
      !
      ! Computes the contributions to internal flux terms and Jacobian (when calc_derivs
      ! is true) arising from pressure gradients perpendicular to cell-centre connections
      ! due to anisotropic permeability terms
      !
      ! Author: Steve Benbow
      ! Date: 01/06/20

      use Option_module
      use Material_Aux_module
      use Global_Aux_module

      ! inputs
      type(aniso_richards_data_type), pointer :: aniso_rich

      class(material_auxvar_type) :: material_auxvar_up, material_auxvar_dn
      type(global_auxvar_type) :: global_auxvar_up, global_auxvar_dn
      type(option_type) :: option
      PetscReal :: dist(-1:3)
      PetscReal :: ukvr

      PetscReal :: dukvr_dp_up, dukvr_dp_dn
      PetscBool :: calc_derivs ! terms are stored in aniso_rich%dv_darcy_aniso_...

      PetscReal :: v_darcy_iso ! v_darcy = Dq * ukvr * dphi, as calculated by RichardsFlux()
      PetscReal :: dv_darcy_iso_dp_up,dv_darcy_iso_dp_dn ! derivs of above

      ! outputs
      PetscReal :: v_darcy_aniso

      ! local
      PetscReal :: mu_up, mu_dn, rho_up, rho_dn
      PetscReal :: v_darcy_aniso_up, v_darcy_aniso_dn

      ! anisotropy work variables
      type(aniso_geom_data_type), pointer :: aniso_geom
      type(aniso_cell_data_type), pointer :: aniso_geomupdn
      type(aniso_richards_work_type), pointer :: aniso_richw
      type(aniso_richards_cell_type), pointer :: aniso_richupdn
      PetscReal :: pres_updn, rho_updn
      PetscReal :: dphi_dv_up, dphi_dv_dn, dphi_dw_up, dphi_dw_dn
      PetscInt :: iupdn, iadj
      PetscReal :: LSrhs(2), dot, res, dres_dp_up, dres_dp_dn, dres_dp_adj
      PetscReal :: dphi_du, d2phi_du_dp_up, d2phi_du_dp_dn

      nullify(aniso_geom)
      nullify(aniso_geomupdn)
      nullify(aniso_richw)
      nullify(aniso_richupdn)

      ! Initialise output
      v_darcy_aniso = 0.d0

      ! Ideally would have these for adjacent cells too
      rho_up = global_auxvar_up%den(1) * FMWH2O
      rho_dn = global_auxvar_dn%den(1) * FMWH2O

      ! set ptrs
      aniso_geom => aniso_rich%aniso_geom

      ! v,w are the coordinate vectors in the plane of the interface over which flux is being calculated
      dphi_dv_up = 0
      dphi_dw_up = 0
      dphi_dv_dn = 0
      dphi_dw_dn = 0

      if (calc_derivs) then
        allocate(aniso_richw)
        allocate(aniso_richw%d2phi_dv_up_dp_up_adj(aniso_geom%up%num_adj))
        allocate(aniso_richw%d2phi_dw_up_dp_up_adj(aniso_geom%up%num_adj))
        allocate(aniso_richw%d2phi_dv_dn_dp_dn_adj(aniso_geom%dn%num_adj))
        allocate(aniso_richw%d2phi_dw_dn_dp_dn_adj(aniso_geom%dn%num_adj))
        aniso_richw%d2phi_dv_up_dp_up     = 0 ! scalar
        aniso_richw%d2phi_dv_up_dp_dn     = 0 ! scalar
        aniso_richw%d2phi_dv_up_dp_up_adj = 0 ! array
        aniso_richw%d2phi_dw_up_dp_up     = 0 ! scalar
        aniso_richw%d2phi_dw_up_dp_dn     = 0 ! scalar
        aniso_richw%d2phi_dw_up_dp_up_adj = 0 ! array
        aniso_richw%d2phi_dv_dn_dp_dn     = 0 ! scalar
        aniso_richw%d2phi_dv_dn_dp_up     = 0 ! scalar
        aniso_richw%d2phi_dv_dn_dp_dn_adj = 0 ! array
        aniso_richw%d2phi_dw_dn_dp_dn     = 0 ! scalar
        aniso_richw%d2phi_dw_dn_dp_up     = 0 ! scalar
        aniso_richw%d2phi_dw_dn_dp_dn_adj = 0 ! array
      endif

      do iupdn = 1, 2

        if (iupdn == 1) then
          ! upstream
          aniso_geomupdn => aniso_geom%up
          aniso_richupdn => aniso_rich%up
          pres_updn = global_auxvar_up%pres(1)
          rho_updn = rho_up
        else
          ! downstream
          aniso_geomupdn => aniso_geom%dn
          aniso_richupdn => aniso_rich%dn
          pres_updn = global_auxvar_dn%pres(1)
          rho_updn = rho_dn
        endif

        ! Calculate dphi_du from v_darcy (= Dq * ukvr * dphi) supplied by RichardsFlux(),
        ! where Dq = (perm_up * perm_dn)/(dd_up*perm_dn + dd_dn*perm_up)
        ! (NB. dphi_du is basically -(P_up-P_dn)/dist) + gravity terms)
        dphi_du = -1.d0* v_darcy_iso / (ukvr*aniso_richupdn%perm_uu)

        if (calc_derivs) then
          d2phi_du_dp_up = -1.d0* dv_darcy_iso_dp_up / (ukvr*aniso_richupdn%perm_uu) &
                            + v_darcy_iso * dukvr_dp_up / (ukvr*ukvr*aniso_richupdn%perm_uu)
          d2phi_du_dp_dn = -1.d0* dv_darcy_iso_dp_dn / (ukvr*aniso_richupdn%perm_uu) &
                            + v_darcy_iso * dukvr_dp_dn / (ukvr*ukvr*aniso_richupdn%perm_uu)
        endif

        ! Compute grad(h) = [alpha_u,alpha_v,alpha_w]' in the up/dn cell.
        ! alpha_u is known from the isotropic calculation.
        ! [alpha_v,alpha_w]' can be obtained by solving the normal equations for the (v,w) components of grad(h),
        !    M'DM [alpha_v,alpha_w]' = Sum_adj (1/dist_adj)^2 * (h_adj - h_up/dn - alpha_u*dist_adj_u) [dist_adj_v,dist_adj_w]'
        ! where adj indicates cells adjacent to up/dn that are off-diagonal anisotropy neighbours,
        ! and dist_adj = [dist_adj_u,dist_adj_v,dist_adj_w]' is the distance vector from the centre of the up/dn
        ! cell to the centre of the adjacent cell in (u,v,w) coordinates.
        ! A least-squares solution is formed since there may be more than 2 adjacent cells, in which case the
        ! approximation to grad(h) can be over-determiend.
        ! The inverse of the 2x2 matrix M'DM is already stored in aniso_geomupdn%LSMatInv (it only depends on
        ! geometry so does not vary throughout the simulation).
        !
        ! Instead forming inv(M'DM)*sum(rhs_terms), forming sum( inv(M'DM)*rhs_term ) allows the derivatives to
        ! be more easily accumulated when calculating Jacobian entries.

        ! Loop over cells adjacent to up/dn
        do iadj = 1, aniso_geomupdn%num_adj
          dot = aniso_geomupdn%dists_vw(1,iadj)*aniso_geomupdn%dists_vw(1,iadj) &
                + aniso_geomupdn%dists_vw(2,iadj)*aniso_geomupdn%dists_vw(2,iadj)
          res = (aniso_richupdn%pres_adj(iadj) - rho_updn*aniso_geomupdn%loc_gravity_adj(iadj) ) &
                - (pres_updn - rho_updn*aniso_geomupdn%loc_gravity) &
                - dphi_du*aniso_geomupdn%dists_u(iadj)

          LSrhs(1) = (1.d0/dot) * aniso_geomupdn%dists_vw(1,iadj)
          LSrhs(2) = (1.d0/dot) * aniso_geomupdn%dists_vw(2,iadj)

          if (iupdn == 1) then
            ! upstream
            ! apply inverse of LS matrix-residual
            dphi_dv_up = dphi_dv_up + aniso_geomupdn%LSMatInv(1,1)*LSrhs(1)*res + aniso_geomupdn%LSMatInv(1,2)*LSrhs(2)*res
            dphi_dw_up = dphi_dw_up + aniso_geomupdn%LSMatInv(2,1)*LSrhs(1)*res + aniso_geomupdn%LSMatInv(2,2)*LSrhs(2)*res

            if (calc_derivs) then
              dres_dp_up  = -1.d0 - d2phi_du_dp_up*aniso_geomupdn%dists_u(iadj)
              dres_dp_dn  =       - d2phi_du_dp_dn*aniso_geomupdn%dists_u(iadj)
              dres_dp_adj =  1.d0

              aniso_richw%d2phi_dv_up_dp_up = aniso_richw%d2phi_dv_up_dp_up                         &
                                            + aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_up &
                                            + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_up
              aniso_richw%d2phi_dv_up_dp_dn = aniso_richw%d2phi_dv_up_dp_dn                         &
                                            + aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_dn &
                                            + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_dn

              aniso_richw%d2phi_dw_up_dp_up = aniso_richw%d2phi_dw_up_dp_up                         &
                                            + aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_up &
                                            + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_up
              aniso_richw%d2phi_dw_up_dp_dn = aniso_richw%d2phi_dw_up_dp_dn                         &
                                            + aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_dn &
                                            + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_dn

              aniso_richw%d2phi_dv_up_dp_up_adj(iadj) = aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_adj &
                                                    + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_adj
              aniso_richw%d2phi_dw_up_dp_up_adj(iadj) = aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_adj &
                                                    + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_adj
            endif
          else
            ! downstream
            ! apply inverse of LS matrix-residual
            dphi_dv_dn = dphi_dv_dn + aniso_geomupdn%LSMatInv(1,1)*LSrhs(1)*res + aniso_geomupdn%LSMatInv(1,2)*LSrhs(2)*res
            dphi_dw_dn = dphi_dw_dn + aniso_geomupdn%LSMatInv(2,1)*LSrhs(1)*res + aniso_geomupdn%LSMatInv(2,2)*LSrhs(2)*res

            if (calc_derivs) then
              dres_dp_up  =       - d2phi_du_dp_up*aniso_geomupdn%dists_u(iadj)
              dres_dp_dn  = -1.d0 - d2phi_du_dp_dn*aniso_geomupdn%dists_u(iadj)
              dres_dp_adj =  1.d0

              aniso_richw%d2phi_dv_dn_dp_dn = aniso_richw%d2phi_dv_dn_dp_dn                         &
                                            + aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_dn &
                                            + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_dn
              aniso_richw%d2phi_dv_dn_dp_up = aniso_richw%d2phi_dv_dn_dp_up                         &
                                            + aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_up &
                                            + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_up

              aniso_richw%d2phi_dw_dn_dp_dn = aniso_richw%d2phi_dw_dn_dp_dn                         &
                                            + aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_dn &
                                            + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_dn
              aniso_richw%d2phi_dw_dn_dp_up = aniso_richw%d2phi_dw_dn_dp_up                         &
                                            + aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_up &
                                            + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_up

              aniso_richw%d2phi_dv_dn_dp_dn_adj(iadj) = aniso_geomupdn%LSMatInv(1,1) * LSrhs(1) * dres_dp_adj &
                                                    + aniso_geomupdn%LSMatInv(1,2) * LSrhs(2) * dres_dp_adj
              aniso_richw%d2phi_dw_dn_dp_dn_adj(iadj) = aniso_geomupdn%LSMatInv(2,1) * LSrhs(1) * dres_dp_adj &
                                                    + aniso_geomupdn%LSMatInv(2,2) * LSrhs(2) * dres_dp_adj
            endif
          endif
        enddo

      enddo

      ! Calculate anisotropic darcy contribution in up/dn cells
      v_darcy_aniso_up = -1.d0 * (aniso_rich%up%perm_uv * dphi_dv_up + aniso_rich%up%perm_uw * dphi_dw_up) * ukvr
      v_darcy_aniso_dn = -1.d0 * (aniso_rich%dn%perm_uv * dphi_dv_dn + aniso_rich%dn%perm_uw * dphi_dw_dn) * ukvr
      mu_up =       dist(-1) *dist(0)/aniso_rich%up%perm_uu  ! weight - up (dist to interface/Kuu)
      mu_dn = (1.d0-dist(-1))*dist(0)/aniso_rich%dn%perm_uu  ! weight - dn (dist to interface/Kuu)
      ! average up/dn components to obtain a consistent interface midpoint pressure
      v_darcy_aniso = (mu_up*v_darcy_aniso_up + mu_dn*v_darcy_aniso_dn) / (mu_up + mu_dn)

      if (calc_derivs) then
        aniso_rich%dv_darcy_aniso_dp_up = (mu_up/(mu_up+mu_dn)) * &
                                        (-1.0d0 * (  aniso_rich%up%perm_uv * aniso_richw%d2phi_dv_up_dp_up &
                                                    + aniso_rich%up%perm_uw * aniso_richw%d2phi_dw_up_dp_up) * ukvr &
                                          -1.0d0 * (aniso_rich%up%perm_uv * dphi_dv_up &
                                                    + aniso_rich%up%perm_uw * dphi_dw_up) * dukvr_dp_up) &
                                          + (mu_dn/(mu_up+mu_dn)) * &
                                          (-1.0d0 * (  aniso_rich%dn%perm_uv * aniso_richw%d2phi_dv_dn_dp_up &
                                                    + aniso_rich%dn%perm_uw * aniso_richw%d2phi_dw_dn_dp_up) * ukvr &
                                          -1.0d0 * (aniso_rich%dn%perm_uv * dphi_dv_dn &
                                                    + aniso_rich%dn%perm_uw * dphi_dw_dn) * dukvr_dp_up)
        aniso_rich%dv_darcy_aniso_dp_dn = (mu_dn/(mu_up+mu_dn)) * &
                                        (-1.0d0 * (  aniso_rich%dn%perm_uv * aniso_richw%d2phi_dv_dn_dp_dn &
                                                    + aniso_rich%dn%perm_uw * aniso_richw%d2phi_dw_dn_dp_dn) * ukvr &
                                          -1.0d0 * (aniso_rich%dn%perm_uv * dphi_dv_dn &
                                                    + aniso_rich%dn%perm_uw * dphi_dw_dn) * dukvr_dp_dn) &
                                          + (mu_up/(mu_up+mu_dn)) * &
                                          (-1.0d0 * (  aniso_rich%up%perm_uv * aniso_richw%d2phi_dv_up_dp_dn &
                                                    + aniso_rich%up%perm_uw * aniso_richw%d2phi_dw_up_dp_dn) * ukvr &
                                          -1.0d0 * (aniso_rich%up%perm_uv * dphi_dv_up &
                                                    + aniso_rich%up%perm_uw * dphi_dw_up) * dukvr_dp_dn)
        do iadj = 1, aniso_geom%up%num_adj
          aniso_rich%dv_darcy_aniso_dp_up_adj(iadj) = (mu_up/(mu_up+mu_dn)) * &
                                    -1.0d0 * (  aniso_rich%up%perm_uv * aniso_richw%d2phi_dv_up_dp_up_adj(iadj) &
                                              + aniso_rich%up%perm_uw * aniso_richw%d2phi_dw_up_dp_up_adj(iadj)) * ukvr
        enddo
        do iadj = 1, aniso_geom%dn%num_adj
          aniso_rich%dv_darcy_aniso_dp_dn_adj(iadj) = (mu_dn/(mu_up+mu_dn)) * &
                                    -1.0d0 * (  aniso_rich%dn%perm_uv * aniso_richw%d2phi_dv_dn_dp_dn_adj(iadj) &
                                              + aniso_rich%dn%perm_uw * aniso_richw%d2phi_dw_dn_dp_dn_adj(iadj)) * ukvr
        enddo
      endif

      if (calc_derivs) then
        deallocate(aniso_richw%d2phi_dv_up_dp_up_adj)
        deallocate(aniso_richw%d2phi_dw_up_dp_up_adj)
        deallocate(aniso_richw%d2phi_dv_dn_dp_dn_adj)
        deallocate(aniso_richw%d2phi_dw_dn_dp_dn_adj)
        deallocate(aniso_richw)
      endif

    end subroutine RichardsFluxAniso

end module Richards_Aniso_module

