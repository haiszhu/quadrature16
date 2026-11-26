subroutine quadSn_stok_off_w16(m,t0,h,C_Xi,C_Xij,R,E,density_map,tnx,W)
  ! QUADSN_STOK_OFF_W16  Stokes SLPn (traction of the SLP at the target normal tnx(:,itarg)) off surface (h~=0), the part
  ! ((x-y).tnx)*Dn_ab: the DLPn pipeline with psi = ((x-y).tnx)*mu formed at the radial nodes
  ! (compute_singular_correction_stok_off, kind 'Sn'). W(m,256,9), (a,b) at a+3*(b-1); quadSD_stok_w16 forms
  ! (1/2)*[(tnx.Sn)*I - W].
  ! Same text as quadDn_off_w16 except the lines marked "! <-".
  use surfquadr
  implicit none
  integer(8), parameter :: p=2, order=2, ncomp=9                              ! <- kernel constants
  logical, parameter :: if_N=.false.                                                             ! <- DLP reads N, the others J
  integer(8), intent(in) :: m
  real(8), intent(in) :: tnx(3,m)                                                              ! <- target normals
  real(8), intent(in) :: t0(2,m), h(m), C_Xi(16,16,3,2), C_Xij(16,16,3,3), density_map(16,16)
  type(radial_tables), intent(in) :: R
  type(edge_tables), intent(in) :: E
  real(8), intent(inout) :: W(m,256,ncomp)
  ! Per-call state: allocated on first use, reused for every target and edge.
  type(source_geometry) :: src
  type(target_geometry) :: geo
  type(radial_expansion_tables) :: T_radial
  type(displacement_expansion) :: X
  type(jacobian_expansion) :: Jx
  type(distance_factor) :: gp
  type(kernel_weights) :: wk
  real(8), allocatable :: d(:,:), w_s(:), G(:,:,:), G1(:,:,:), Gr(:,:,:)
  real(8) :: G0(1,6), W_leg(16,16,ncomp), weight_0(ncomp), weight_u(ncomp), weight_v(ncomp)
  integer(8), parameter :: edge_order(4)=[1_8,3_8,2_8,4_8]   ! opposite edges in pairs: edges 3, 4 reuse the step-3 tables of edges 1, 2
  integer(8) :: itarg, k, iedge, n_E

  call compute_source_geometry(C_Xi,src)                                                      ! 0
  do itarg=1,m
    call compute_target_geometry(t0(:,itarg),h(itarg),C_Xi,C_Xij,geo)                         ! 1
    W_leg=0.0d0; weight_0=0.0d0; weight_u=0.0d0; weight_v=0.0d0   ! weights on mu(t0), d_u mu(t0), d_v mu(t0)
    do k=1,4
      iedge=edge_order(k)
      call compute_edge_rule(iedge,geo,E,n_E,d,w_s)                                           ! 2
      call compute_radial_expansion(order,iedge,geo,src,R,n_E,d,G,G0,G1,Gr,T_radial)          ! 3
      call compute_displacement_expansion(order,src,R,n_E,d,G,G0,G1,Gr,X)   ! C0,q0,c1r | V1r | c1,cr  ! 4
      call compute_jacobian_expansion(order,if_N,R,n_E,G,G0,G1,Gr,Jx)    ! N,J | N1,Nr,J1,J1r,Jr  ! 5
      call compute_distance_factor(p,order,R,geo,n_E,X,gp)          ! g0,g | g1r | g1,gr,beta1  ! 6
      call compute_kernelradial_Dn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)  ! 7 <-
      call compute_singular_correction_stok_off('Sn',p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u, &
        weight_v,tnx(:,itarg))  ! 8 <-
      call compute_muhat_weights(order,iedge,R,T_radial,n_E,d,wk,W_leg)                              ! 9
    enddo
    call compute_final_weights(density_map,geo,weight_0,weight_u,weight_v,W_leg,W(itarg,:,:))                            ! 10
  enddo
end subroutine quadSn_stok_off_w16
