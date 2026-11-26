! Layer potential quadrature for one GL16 source patch.

subroutine quadSD_w16(m,tx,tnx,n,sx,snx,WS,WD,WSn,WDn,t0,h,sw,center, &
    C_X,C_Xi,C_Xij,C_X_fit,cols,target_idx,other,onsurface_comp,ib,nb,xb,cols_typeb, &
    cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,ifS,ifD,ifSn,ifDn,lwork,work,ier)
  ! QUADSD_W16  Near/self and mid weights of one source patch (n=256 GL16 nodes) for m targets.
  ! Near/self (max|t0|<=1.1): h=0 rows use quad{S,D,Sn,Dn}_on_w16, h~=0 rows use quad{S,D,Sn,Dn}_off_w16.
  ! Mid: quadSD_mid_w16 evaluates SLP and DLP, plus SLPn and DLPn when requested.
  ! lwork and work are unused.
  use surfquadr, only: radial_tables, edge_tables, compute_gauss_legendre
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  integer(8), intent(in) :: m,n,nb,cols(n),target_idx(m),other(m),onsurface_comp,ib
  integer(8), intent(in) :: cols_typeb(2304,2),if_h0(m),ifS,ifD,ifSn,ifDn,lwork
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),center(3)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3),C_X_fit(16,16,3)
  real(8), intent(in) :: xb(3,nb),cxb(16,16,3,2),cxpb(16,16,3,2,2),cxppb(16,16,3,3,2)
  real(8), intent(in) :: tcen(2),tsca(2),tau0(2,m)
  real(8), intent(inout) :: WS(m,n*ifS),WD(m,n*ifD),WSn(m,n*ifSn,3),WDn(m,n*ifDn,3,3)
  real(8), intent(inout) :: t0(2,m),h(m),work(lwork)
  integer(8), intent(inout) :: ier(m)
  type(radial_tables) :: R
  type(edge_tables) :: E
  integer(8) :: i,j,node,m_near,m_on,m_off,m_mid,near_rows(m),on_rows(m),off_rows(m),mid_rows(m)
  real(8) :: g16(16),wg16(16),st(2,256),density_map(16,16),nan
  logical(c_bool) :: root_ok(m)
  external :: quadS_on_w16,quadD_on_w16,quadSn_on_w16,quadDn_on_w16,quadS_off_w16,quadD_off_w16,quadSn_off_w16,quadDn_off_w16, &
    quadSD_mid_w16

  ! 1. Quadrature preparation
  call compute_gauss_legendre(16_8,g16,wg16)
  call prepare_tables_w16(g16,R,E,density_map)
  do j=1,16
    do i=1,16
      st(:,i+16*(j-1))=[g16(i),g16(j)]
    enddo
  enddo

  ! 2. Projection
  root_ok=all(ieee_is_finite(t0),dim=1) .and. ieee_is_finite(h)
  call project_t0_h_w16(m,n,nb,n,sx,st,cols,C_X_fit,C_Xi,C_Xij,tx,target_idx,other,onsurface_comp,ib, &
      xb,g16,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,sx,st,[1.1d0,1.1d0],if_h0,t0,h,root_ok)

  ! 3. Classification
  ier=2_8
  m_near=0_8
  m_mid=0_8
  do i=1,m
    if (root_ok(i) .and. all(ieee_is_finite(t0(:,i))) .and. ieee_is_finite(h(i)) &
        .and. maxval(abs(t0(:,i)))<=1.1d0) then
      m_near=m_near+1_8
      near_rows(m_near)=i
    else
      m_mid=m_mid+1_8
      mid_rows(m_mid)=i
    endif
  enddo

  ! 4. Near/self. The off-surface drivers take the negative of the projected signed height h.
  m_on=0_8
  m_off=0_8
  do i=1,m_near
    node=near_rows(i)
    if (h(node)==0.0d0) then
      m_on=m_on+1_8
      on_rows(m_on)=node
    else
      m_off=m_off+1_8
      off_rows(m_off)=node
    endif
  enddo
  if (m_off>0_8) then
    block
      ! W_off(:,:,1:ncomp) is each driver's W(m_off,256,ncomp); allocated, since m_off*256*9 can exceed the stack.
      real(8) :: t0_off(2,m_off),h_off(m_off)
      real(8), allocatable :: W_off(:,:,:)
      integer(8) :: ia,ik
      allocate(W_off(m_off,256,9))
      t0_off=t0(:,off_rows(1:m_off))
      h_off=-h(off_rows(1:m_off))
      if (ifS/=0_8) then
        call quadS_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off)
        WS(off_rows(1:m_off),:)=W_off(:,:,1)
      endif
      if (ifD/=0_8) then
        call quadD_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off)
        WD(off_rows(1:m_off),:)=W_off(:,:,1)
      endif
      if (ifSn/=0_8) then
        call quadSn_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off)
        WSn(off_rows(1:m_off),:,:)=W_off(:,:,1:3)
      endif
      if (ifDn/=0_8) then
        ! 9 components, (a,k) at a+3*(k-1); off surface the discrete tensor is not symmetric
        call quadDn_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off)
        do ik=1,3
          do ia=1,3
            WDn(off_rows(1:m_off),:,ia,ik)=W_off(:,:,ia+3*(ik-1))
          enddo
        enddo
      endif
      ier(off_rows(1:m_off))=0_8
    end block
  endif
  if (m_on>0_8) then
    block
      ! W_on(:,:,1:ncomp) is each driver's W(m_on,256,ncomp); allocated, since m_on*256*6 can exceed the stack.
      real(8) :: t0_on(2,m_on),h_on(m_on)
      real(8), allocatable :: W_on(:,:,:)
      integer(8) :: k,ic,jc
      allocate(W_on(m_on,256,6))
      t0_on=t0(:,on_rows(1:m_on))
      h_on=-h(on_rows(1:m_on))
      if (ifS/=0_8) then
        call quadS_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        WS(on_rows(1:m_on),:)=W_on(:,:,1)
      endif
      if (ifD/=0_8) then
        call quadD_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        WD(on_rows(1:m_on),:)=W_on(:,:,1)
      endif
      if (ifSn/=0_8) then
        call quadSn_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        WSn(on_rows(1:m_on),:,:)=W_on(:,:,1:3)
      endif
      if (ifDn/=0_8) then
        ! 6 components (11,12,13,22,23,33) -> symmetric 3x3
        call quadDn_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        k=0_8
        do ic=1,3
          do jc=ic,3
            k=k+1_8
            WDn(on_rows(1:m_on),:,ic,jc)=W_on(:,:,k)
            WDn(on_rows(1:m_on),:,jc,ic)=W_on(:,:,k)
          enddo
        enddo
      endif
      ier(on_rows(1:m_on))=0_8
    end block
  endif
  ! A non-finite weight (e.g. J0=0 at t0) sets ier=6 for that row.
  if (m_near>0_8) then
    block
      integer(8) :: ic,jc
      do j=1,n*ifS
        do i=1,m_near
          if (.not.ieee_is_finite(WS(near_rows(i),j))) ier(near_rows(i))=6_8
        enddo
      enddo
      do j=1,n*ifD
        do i=1,m_near
          if (.not.ieee_is_finite(WD(near_rows(i),j))) ier(near_rows(i))=6_8
        enddo
      enddo
      do ic=1,3*ifSn
        do j=1,n
          do i=1,m_near
            if (.not.ieee_is_finite(WSn(near_rows(i),j,ic))) ier(near_rows(i))=6_8
          enddo
        enddo
      enddo
      do jc=1,3*ifDn
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.ieee_is_finite(WDn(near_rows(i),j,ic,jc))) ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
    end block
  endif

  ! 5. Mid: SLP and DLP always; SLPn, DLPn when requested.
  if (m_mid>0_8) then
    block
      real(8) :: tx_mid(3,m_mid),tnx_mid(3,m_mid),sx_mid(3,n),t0_mid(2,m_mid),h_mid(m_mid)
      real(8) :: WS_mid(m_mid,n),WD_mid(m_mid,n),geometry_center(2),geometry_scale(2)
      real(8) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3),geometry_center_typeb(2),geometry_scale_typeb(2)
      real(8), allocatable :: WSn_mid(:,:,:),WDn_mid(:,:,:,:)   ! m_mid*256*9 can exceed the stack
      integer(8) :: classification_mid(m_mid),ier_mid(m_mid),ic,jc
      tx_mid=tx(:,mid_rows(1:m_mid))-spread(center,2,m_mid)
      sx_mid=sx-spread(center,2,n)
      tnx_mid=tnx(:,mid_rows(1:m_mid))
      t0_mid=t0(:,mid_rows(1:m_mid))
      h_mid=h(mid_rows(1:m_mid))
      if (onsurface_comp/=0_8) then
        where (other(mid_rows(1:m_mid))==0_8 .and. .not.ieee_is_finite(h_mid)) h_mid=0.0d0
      endif
      classification_mid=0_8; ier_mid=0_8
      nan=ieee_value(0.0d0,ieee_quiet_nan)
      allocate(WSn_mid(m_mid,n*ifSn,3),WDn_mid(m_mid,n*ifDn,3,3))
      WS_mid=nan; WD_mid=nan; WSn_mid=nan; WDn_mid=nan
      geometry_center=0.0d0; geometry_scale=1.0d0
      C_Xi_typeb=0.0d0; C_Xij_typeb=0.0d0
      geometry_center_typeb=0.0d0; geometry_scale_typeb=0.0d0
      if (ib/=0_8) then
        C_Xi_typeb(:,:,:,1)=cxpb(:,:,:,1,ib)/tsca(1)
        C_Xi_typeb(:,:,:,2)=cxpb(:,:,:,2,ib)/tsca(2)
        C_Xij_typeb(:,:,:,1)=cxppb(:,:,:,1,ib)/tsca(1)**2
        C_Xij_typeb(:,:,:,2)=cxppb(:,:,:,2,ib)/(tsca(1)*tsca(2))
        C_Xij_typeb(:,:,:,3)=cxppb(:,:,:,3,ib)/tsca(2)**2
        geometry_center_typeb=tcen
        geometry_scale_typeb=tsca
      endif
      call quadSD_mid_w16(m_mid,tx_mid,tnx_mid,n,sx_mid,snx,sw,WS_mid,WD_mid,ifSn,ifDn,WSn_mid,WDn_mid, &
          t0_mid,h_mid,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification_mid,C_Xi_typeb,C_Xij_typeb, &
          geometry_center_typeb,geometry_scale_typeb,ier_mid,root_ok(mid_rows(1:m_mid)))
      ! Each density column is contiguous in every mid output matrix.
      where (ier_mid==1_8) ier_mid=2_8
      do j=1,n
        do i=1,m_mid
          if (ier_mid(i)==0_8 .and. (.not.ieee_is_finite(WS_mid(i,j)) .or. &
                                      .not.ieee_is_finite(WD_mid(i,j)))) ier_mid(i)=6_8
        enddo
      enddo
      do ic=1,3*ifSn
        do j=1,n
          do i=1,m_mid
            if (ier_mid(i)==0_8 .and. .not.ieee_is_finite(WSn_mid(i,j,ic))) ier_mid(i)=6_8
          enddo
        enddo
      enddo
      do jc=1,3*ifDn
        do ic=1,3
          do j=1,n
            do i=1,m_mid
              if (ier_mid(i)==0_8 .and. .not.ieee_is_finite(WDn_mid(i,j,ic,jc))) ier_mid(i)=6_8
            enddo
          enddo
        enddo
      enddo
      ier(mid_rows(1:m_mid))=ier_mid
      if (all(ier_mid==0_8)) then
        if (ifS/=0_8) WS(mid_rows(1:m_mid),:)=WS_mid
        if (ifD/=0_8) WD(mid_rows(1:m_mid),:)=WD_mid
        if (ifSn/=0_8) WSn(mid_rows(1:m_mid),:,:)=WSn_mid
        if (ifDn/=0_8) WDn(mid_rows(1:m_mid),:,:,:)=WDn_mid
      else
        do i=1,m_mid
          if (ier_mid(i)/=0_8) cycle
          node=mid_rows(i)
          if (ifS/=0_8) WS(node,:)=WS_mid(i,:)
          if (ifD/=0_8) WD(node,:)=WD_mid(i,:)
          if (ifSn/=0_8) WSn(node,:,:)=WSn_mid(i,:,:)
          if (ifDn/=0_8) WDn(node,:,:,:)=WDn_mid(i,:,:,:)
        enddo
      endif
    end block
  endif

  ! 6. Unknown rows retained
  nan=ieee_value(0.0d0,ieee_quiet_nan)
  do i=1,m
    if (root_ok(i)) cycle
    t0(:,i)=nan
    h(i)=nan
  enddo
end subroutine quadSD_w16

subroutine quadSD_helm_w16(m,tx,tnx,n,sx,snx,WS,WD,WSn,WDn,t0,h,sw,center, &
    C_X,C_Xi,C_Xij,C_X_fit,cols,target_idx,other,onsurface_comp,ib,nb,xb,cols_typeb, &
    cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,ifS,ifD,ifSn,ifDn,zk,lwork,work,ier)
  ! QUADSD_HELM_W16  Helmholtz near/self and mid weights of one source patch (n=256 GL16 nodes) for m targets,
  ! wavenumber zk (complex allowed). quadSD_w16 with zk after ifDn; WS, WD, WSn, WDn complex, in its 1/(2*pi)
  ! normalization and component layout (zk = 0 gives its weights).
  ! Near/self (max|t0|<=1.1): h=0 rows use quad{S,D,Sn,Dn}_helm_on_w16, h~=0 rows quad{S,D,Sn,Dn}_helm_off_w16 (real,
  ! imaginary parts as components; off surface accurate while |Im zk| times the patch size is small: the split of step 7
  ! cancels as exp(2*|Im zk|*dist)). Mid: quadSD_helm_mid_w16.
  ! lwork and work are unused.
  use surfquadr, only: radial_tables, edge_tables, compute_gauss_legendre
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  integer(8), intent(in) :: m,n,nb,cols(n),target_idx(m),other(m),onsurface_comp,ib
  integer(8), intent(in) :: cols_typeb(2304,2),if_h0(m),ifS,ifD,ifSn,ifDn,lwork
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),center(3)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3),C_X_fit(16,16,3)
  real(8), intent(in) :: xb(3,nb),cxb(16,16,3,2),cxpb(16,16,3,2,2),cxppb(16,16,3,3,2)
  real(8), intent(in) :: tcen(2),tsca(2),tau0(2,m)
  complex(8), intent(in) :: zk
  complex(8), intent(inout) :: WS(m,n*ifS),WD(m,n*ifD),WSn(m,n*ifSn,3),WDn(m,n*ifDn,3,3)
  real(8), intent(inout) :: t0(2,m),h(m),work(lwork)
  integer(8), intent(inout) :: ier(m)
  type(radial_tables) :: R
  type(edge_tables) :: E
  integer(8) :: i,j,node,m_near,m_on,m_off,m_mid,near_rows(m),on_rows(m),off_rows(m),mid_rows(m)
  real(8) :: g16(16),wg16(16),st(2,256),density_map(16,16),nan
  logical(c_bool) :: root_ok(m)
  external :: quadS_helm_on_w16,quadD_helm_on_w16,quadSn_helm_on_w16,quadDn_helm_on_w16,quadSD_helm_mid_w16
  external :: quadS_helm_off_w16,quadD_helm_off_w16,quadSn_helm_off_w16,quadDn_helm_off_w16

  ! 1. Quadrature preparation
  call compute_gauss_legendre(16_8,g16,wg16)
  call prepare_tables_w16(g16,R,E,density_map)
  do j=1,16
    do i=1,16
      st(:,i+16*(j-1))=[g16(i),g16(j)]
    enddo
  enddo

  ! 2. Projection
  root_ok=all(ieee_is_finite(t0),dim=1) .and. ieee_is_finite(h)
  call project_t0_h_w16(m,n,nb,n,sx,st,cols,C_X_fit,C_Xi,C_Xij,tx,target_idx,other,onsurface_comp,ib, &
      xb,g16,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,sx,st,[1.1d0,1.1d0],if_h0,t0,h,root_ok)

  ! 3. Classification
  ier=2_8
  m_near=0_8
  m_mid=0_8
  do i=1,m
    if (root_ok(i) .and. all(ieee_is_finite(t0(:,i))) .and. ieee_is_finite(h(i)) &
        .and. maxval(abs(t0(:,i)))<=1.1d0) then
      m_near=m_near+1_8
      near_rows(m_near)=i
    else
      m_mid=m_mid+1_8
      mid_rows(m_mid)=i
    endif
  enddo

  ! 4. Near/self. The off-surface drivers take the negative of the projected signed height h. Driver components:
  ! real parts 1:ncomp_L, imaginary parts ncomp_L+1:2*ncomp_L.
  m_on=0_8
  m_off=0_8
  do i=1,m_near
    node=near_rows(i)
    if (h(node)==0.0d0) then
      m_on=m_on+1_8
      on_rows(m_on)=node
    else
      m_off=m_off+1_8
      off_rows(m_off)=node
    endif
  enddo
  if (m_off>0_8) then
    block
      ! W_off(:,:,1:ncomp) is each driver's W(m_off,256,ncomp); allocated, since m_off*256*18 can exceed the stack.
      real(8) :: t0_off(2,m_off),h_off(m_off)
      real(8), allocatable :: W_off(:,:,:)
      integer(8) :: ia,ik
      allocate(W_off(m_off,256,18))
      t0_off=t0(:,off_rows(1:m_off))
      h_off=-h(off_rows(1:m_off))
      if (ifS/=0_8) then
        call quadS_helm_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,zk,W_off)
        WS(off_rows(1:m_off),:)=cmplx(W_off(:,:,1),W_off(:,:,2),kind=8)
      endif
      if (ifD/=0_8) then
        call quadD_helm_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,zk,W_off)
        WD(off_rows(1:m_off),:)=cmplx(W_off(:,:,1),W_off(:,:,2),kind=8)
      endif
      if (ifSn/=0_8) then
        call quadSn_helm_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,zk,W_off)
        WSn(off_rows(1:m_off),:,:)=cmplx(W_off(:,:,1:3),W_off(:,:,4:6),kind=8)
      endif
      if (ifDn/=0_8) then
        ! 9 components, (a,k) at a+3*(k-1); off surface the discrete tensor is not symmetric
        call quadDn_helm_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,zk,W_off)
        do ik=1,3
          do ia=1,3
            WDn(off_rows(1:m_off),:,ia,ik)=cmplx(W_off(:,:,ia+3*(ik-1)),W_off(:,:,9+ia+3*(ik-1)),kind=8)
          enddo
        enddo
      endif
      ier(off_rows(1:m_off))=0_8
    end block
  endif
  if (m_on>0_8) then
    block
      ! W_on(:,:,1:ncomp) is each driver's W(m_on,256,ncomp); allocated, since m_on*256*12 can exceed the stack.
      real(8) :: t0_on(2,m_on),h_on(m_on)
      real(8), allocatable :: W_on(:,:,:)
      integer(8) :: k,ic,jc
      allocate(W_on(m_on,256,12))
      t0_on=t0(:,on_rows(1:m_on))
      h_on=-h(on_rows(1:m_on))
      if (ifS/=0_8) then
        call quadS_helm_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,zk,W_on)
        WS(on_rows(1:m_on),:)=cmplx(W_on(:,:,1),W_on(:,:,2),kind=8)
      endif
      if (ifD/=0_8) then
        call quadD_helm_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,zk,W_on)
        WD(on_rows(1:m_on),:)=cmplx(W_on(:,:,1),W_on(:,:,2),kind=8)
      endif
      if (ifSn/=0_8) then
        call quadSn_helm_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,zk,W_on)
        WSn(on_rows(1:m_on),:,:)=cmplx(W_on(:,:,1:3),W_on(:,:,4:6),kind=8)
      endif
      if (ifDn/=0_8) then
        ! 6 components (11,12,13,22,23,33) -> symmetric 3x3
        call quadDn_helm_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,zk,W_on)
        k=0_8
        do ic=1,3
          do jc=ic,3
            k=k+1_8
            WDn(on_rows(1:m_on),:,ic,jc)=cmplx(W_on(:,:,k),W_on(:,:,6+k),kind=8)
            WDn(on_rows(1:m_on),:,jc,ic)=WDn(on_rows(1:m_on),:,ic,jc)
          enddo
        enddo
      endif
      ier(on_rows(1:m_on))=0_8
    end block
  endif
  ! A non-finite weight (e.g. J0=0 at t0) sets ier=6 for that row.
  if (m_near>0_8) then
    block
      integer(8) :: ic,jc
      do j=1,n*ifS
        do i=1,m_near
          if (.not.(ieee_is_finite(real(WS(near_rows(i),j),8)) .and. ieee_is_finite(aimag(WS(near_rows(i),j))))) &
            ier(near_rows(i))=6_8
        enddo
      enddo
      do j=1,n*ifD
        do i=1,m_near
          if (.not.(ieee_is_finite(real(WD(near_rows(i),j),8)) .and. ieee_is_finite(aimag(WD(near_rows(i),j))))) &
            ier(near_rows(i))=6_8
        enddo
      enddo
      do ic=1,3*ifSn
        do j=1,n
          do i=1,m_near
            if (.not.(ieee_is_finite(real(WSn(near_rows(i),j,ic),8)) .and. ieee_is_finite(aimag(WSn(near_rows(i),j,ic))))) &
              ier(near_rows(i))=6_8
          enddo
        enddo
      enddo
      do jc=1,3*ifDn
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.(ieee_is_finite(real(WDn(near_rows(i),j,ic,jc),8)) .and. &
                  ieee_is_finite(aimag(WDn(near_rows(i),j,ic,jc))))) &
                ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
    end block
  endif

  ! 5. Mid: SLP and DLP always; SLPn, DLPn when requested.
  if (m_mid>0_8) then
    block
      real(8) :: tx_mid(3,m_mid),tnx_mid(3,m_mid),sx_mid(3,n),t0_mid(2,m_mid),h_mid(m_mid)
      real(8) :: geometry_center(2),geometry_scale(2)
      real(8) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3),geometry_center_typeb(2),geometry_scale_typeb(2)
      complex(8), allocatable :: WS_mid(:,:),WD_mid(:,:),WSn_mid(:,:,:),WDn_mid(:,:,:,:)   ! can exceed the stack
      integer(8) :: classification_mid(m_mid),ier_mid(m_mid),ic,jc
      tx_mid=tx(:,mid_rows(1:m_mid))-spread(center,2,m_mid)
      sx_mid=sx-spread(center,2,n)
      tnx_mid=tnx(:,mid_rows(1:m_mid))
      t0_mid=t0(:,mid_rows(1:m_mid))
      h_mid=h(mid_rows(1:m_mid))
      if (onsurface_comp/=0_8) then
        where (other(mid_rows(1:m_mid))==0_8 .and. .not.ieee_is_finite(h_mid)) h_mid=0.0d0
      endif
      classification_mid=0_8; ier_mid=0_8
      nan=ieee_value(0.0d0,ieee_quiet_nan)
      allocate(WS_mid(m_mid,n),WD_mid(m_mid,n),WSn_mid(m_mid,n*ifSn,3),WDn_mid(m_mid,n*ifDn,3,3))
      WS_mid=cmplx(nan,nan,kind=8); WD_mid=WS_mid(1,1); WSn_mid=WS_mid(1,1); WDn_mid=WS_mid(1,1)
      geometry_center=0.0d0; geometry_scale=1.0d0
      C_Xi_typeb=0.0d0; C_Xij_typeb=0.0d0
      geometry_center_typeb=0.0d0; geometry_scale_typeb=0.0d0
      if (ib/=0_8) then
        C_Xi_typeb(:,:,:,1)=cxpb(:,:,:,1,ib)/tsca(1)
        C_Xi_typeb(:,:,:,2)=cxpb(:,:,:,2,ib)/tsca(2)
        C_Xij_typeb(:,:,:,1)=cxppb(:,:,:,1,ib)/tsca(1)**2
        C_Xij_typeb(:,:,:,2)=cxppb(:,:,:,2,ib)/(tsca(1)*tsca(2))
        C_Xij_typeb(:,:,:,3)=cxppb(:,:,:,3,ib)/tsca(2)**2
        geometry_center_typeb=tcen
        geometry_scale_typeb=tsca
      endif
      call quadSD_helm_mid_w16(m_mid,tx_mid,tnx_mid,n,sx_mid,snx,sw,WS_mid,WD_mid,ifSn,ifDn,WSn_mid,WDn_mid, &
          t0_mid,h_mid,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification_mid,C_Xi_typeb,C_Xij_typeb, &
          geometry_center_typeb,geometry_scale_typeb,ier_mid,root_ok(mid_rows(1:m_mid)),zk)
      where (ier_mid==1_8) ier_mid=2_8
      do j=1,n
        do i=1,m_mid
          if (ier_mid(i)==0_8 .and. (.not.(ieee_is_finite(real(WS_mid(i,j),8)) .and. ieee_is_finite(aimag(WS_mid(i,j)))) .or. &
              .not.(ieee_is_finite(real(WD_mid(i,j),8)) .and. ieee_is_finite(aimag(WD_mid(i,j)))))) &
              ier_mid(i)=6_8
        enddo
      enddo
      do ic=1,3*ifSn
        do j=1,n
          do i=1,m_mid
            if (ier_mid(i)==0_8 .and. .not.(ieee_is_finite(real(WSn_mid(i,j,ic),8)) .and. ieee_is_finite(aimag(WSn_mid(i,j,ic))))) &
              ier_mid(i)=6_8
          enddo
        enddo
      enddo
      do jc=1,3*ifDn
        do ic=1,3
          do j=1,n
            do i=1,m_mid
              if (ier_mid(i)==0_8 .and. .not.(ieee_is_finite(real(WDn_mid(i,j,ic,jc),8)) .and. &
                  ieee_is_finite(aimag(WDn_mid(i,j,ic,jc))))) ier_mid(i)=6_8
            enddo
          enddo
        enddo
      enddo
      ier(mid_rows(1:m_mid))=ier_mid
      do i=1,m_mid
        if (ier_mid(i)/=0_8) cycle
        node=mid_rows(i)
        if (ifS/=0_8) WS(node,:)=WS_mid(i,:)
        if (ifD/=0_8) WD(node,:)=WD_mid(i,:)
        if (ifSn/=0_8) WSn(node,:,:)=WSn_mid(i,:,:)
        if (ifDn/=0_8) WDn(node,:,:,:)=WDn_mid(i,:,:,:)
      enddo
    end block
  endif

  ! 6. Unknown rows retained
  nan=ieee_value(0.0d0,ieee_quiet_nan)
  do i=1,m
    if (root_ok(i)) cycle
    t0(:,i)=nan
    h(i)=nan
  enddo

end subroutine quadSD_helm_w16

subroutine quadSD_stok_w16(m,tx,tnx,n,sx,snx,WS,WD,WSn,WDn,t0,h,sw,center, &
    C_X,C_Xi,C_Xij,C_X_fit,cols,target_idx,other,onsurface_comp,ib,nb,xb,cols_typeb, &
    cxb,cxpb,cxppb,tcen,tsca,tau0,if_h0,ifS,ifD,ifSn,ifDn,lwork,work,ier)
  ! QUADSD_STOK_W16  Stokes near/self and mid weights of one source patch (n=256 GL16 nodes) for m targets.
  ! quadSD_w16's arguments; WS, WD, WSn, WDn (m,n,3,3), weights on mu_j. With r = x-y: SLP (I/r+r*r'/r^3)/(8*pi), DLP
  ! -3*((x-y).snx)*r*r'/(4*pi*r^5) (the negative of Sto3dDLPmat);
  ! SLPn and DLPn are the tractions of the SLP and DLP at the target normal tnx. DLPn is -An of Sto3dDLPmat with source
  ! normal snx, all 9 components (ic,jc) without symmetry; on surface its finite part, no jump term.
  ! Near/self on surface (h=0): quad{S,D,Sn,Dn}_stok_on_w16 (Sn and Dn at the chart normal n_r0, scaled by tnx.n_r0).
  ! Near/self off surface: SLP, DLP, SLPn from the Laplace off drivers S, D, Sn and quad{S,D,Sn}_stok_off_w16, through
  !   (I/r+r*r'/r^3)/(8*pi) = (1/4)*[S*I + Sn_a*(y-x)_b],
  !   -3*((x-y).n)*r*r'/(4*pi*r^5) = (1/2)*[D*I - ((x-y).n)*Dn_ab],   n = snx (WD's normal) for DLP, tnx for SLPn,
  !   -3*((x-y).tnx)*r*r'/(4*pi*r^5) = (1/2)*[(tnx.Sn)*I - ((x-y).tnx)*Dn_ab],   r = x-y,
  ! S, D, Sn, Dn the Laplace kernels of quadSD_w16 (1/(2*pi); Sn the target gradient of S, Dn the target Hessian of S),
  ! the smooth factors (y-x)_b, (x-y).n applied at the radial nodes (inside quad{S,D,Sn}_stok_off_w16: the weights integrate
  ! kernel*factor*I16[mu], not kernel*I16[factor*mu]); DLPn from quadDn_stok_off_w16 at tnx.
  ! Mid rows: quadSD_stok_mid_w16, the four Stokes kernels at the integration nodes of each region.
  ! lwork and work are unused.
  use surfquadr, only: radial_tables, edge_tables, compute_gauss_legendre, compute_legendre_values
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  integer(8), intent(in) :: m,n,nb,cols(n),target_idx(m),other(m),onsurface_comp,ib
  integer(8), intent(in) :: cols_typeb(2304,2),if_h0(m),ifS,ifD,ifSn,ifDn,lwork
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),center(3)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3),C_X_fit(16,16,3)
  real(8), intent(in) :: xb(3,nb),cxb(16,16,3,2),cxpb(16,16,3,2,2),cxppb(16,16,3,3,2)
  real(8), intent(in) :: tcen(2),tsca(2),tau0(2,m)
  real(8), intent(inout) :: WS(m,n*ifS,3,3),WD(m,n*ifD,3,3),WSn(m,n*ifSn,3,3),WDn(m,n*ifDn,3,3)
  real(8), intent(inout) :: t0(2,m),h(m),work(lwork)
  integer(8), intent(inout) :: ier(m)
  type(radial_tables) :: R
  type(edge_tables) :: E
  integer(8) :: i,j,node,m_near,m_on,m_off,m_mid,near_rows(m),on_rows(m),off_rows(m),mid_rows(m)
  real(8) :: g16(16),wg16(16),st(2,256),density_map(16,16),nan
  logical(c_bool) :: root_ok(m)
  external :: quadS_stok_on_w16,quadD_stok_on_w16,quadSn_stok_on_w16,quadS_off_w16,quadD_off_w16,quadSn_off_w16, &
    quadDn_stok_on_w16,quadS_stok_off_w16,quadD_stok_off_w16,quadSn_stok_off_w16,quadDn_stok_off_w16,quadSD_stok_mid_w16

  ! 1. Quadrature preparation
  call compute_gauss_legendre(16_8,g16,wg16)
  call prepare_tables_w16(g16,R,E,density_map)
  do j=1,16
    do i=1,16
      st(:,i+16*(j-1))=[g16(i),g16(j)]
    enddo
  enddo

  ! 2. Projection
  root_ok=all(ieee_is_finite(t0),dim=1) .and. ieee_is_finite(h)
  call project_t0_h_w16(m,n,nb,n,sx,st,cols,C_X_fit,C_Xi,C_Xij,tx,target_idx,other,onsurface_comp,ib, &
      xb,g16,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,sx,st,[1.1d0,1.1d0],if_h0,t0,h,root_ok)

  ! 3. Classification
  ier=2_8
  m_near=0_8
  m_mid=0_8
  do i=1,m
    if (root_ok(i) .and. all(ieee_is_finite(t0(:,i))) .and. ieee_is_finite(h(i)) &
        .and. maxval(abs(t0(:,i)))<=1.1d0) then
      m_near=m_near+1_8
      near_rows(m_near)=i
    else
      m_mid=m_mid+1_8
      mid_rows(m_mid)=i
    endif
  enddo

  ! 4. Near/self. The off-surface drivers take the negative of the projected signed height h.
  m_on=0_8
  m_off=0_8
  do i=1,m_near
    node=near_rows(i)
    if (h(node)==0.0d0) then
      m_on=m_on+1_8
      on_rows(m_on)=node
    else
      m_off=m_off+1_8
      off_rows(m_off)=node
    endif
  enddo
  if (m_off>0_8 .and. (ifS/=0_8 .or. ifD/=0_8 .or. ifSn/=0_8 .or. ifDn/=0_8)) then
    block
      ! W_off(:,:,k): Laplace S (k=1), D (2), Sn (3:5); quadS_stok_off_w16 (6:14), quadD_stok_off_w16 (15:23),
      ! quadSn_stok_off_w16 (24:32), each (a,b) at a+3*(b-1): the smooth factors applied at the radial nodes.
      ! W_Dn(:,:,1:9): (ic,jc) at ic+3*(jc-1), no symmetry. Allocated, since these arrays can exceed the stack.
      real(8) :: t0_off(2,m_off),h_off(m_off),tnx_off(3,m_off),tnx_W_Sn
      real(8), allocatable :: W_off(:,:,:),W_Dn(:,:,:)
      integer(8) :: ic,jc,l,k
      if (ifS/=0_8 .or. ifD/=0_8 .or. ifSn/=0_8) allocate(W_off(m_off,256,32))
      t0_off=t0(:,off_rows(1:m_off))
      h_off=-h(off_rows(1:m_off))
      tnx_off=tnx(:,off_rows(1:m_off))
      if (ifS/=0_8) then
        call quadS_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off(:,:,1:1))
        call quadS_stok_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off(:,:,6:14))
      endif
      if (ifD/=0_8) then
        call quadD_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off(:,:,2:2))
        call quadD_stok_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off(:,:,15:23))
      endif
      if (ifSn/=0_8) then
        call quadSn_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,W_off(:,:,3:5))
        call quadSn_stok_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,tnx_off,W_off(:,:,24:32))
      endif
      if (ifDn/=0_8) then
        allocate(W_Dn(m_off,256,9))
        call quadDn_stok_off_w16(m_off,t0_off,h_off,C_Xi,C_Xij,R,E,density_map,tnx_off,W_Dn)
      endif
      do l=1,n
        do i=1,m_off
          node=off_rows(i)
          if (ifSn/=0_8) tnx_W_Sn=(W_off(i,l,3)*tnx(1,node)+W_off(i,l,4)*tnx(2,node))+W_off(i,l,5)*tnx(3,node)
          do jc=1,3
            do ic=1,3
              k=ic+3*(jc-1)
              if (ifS/=0_8) WS(node,l,ic,jc)=0.25d0*(merge(W_off(i,l,1),0.0d0,ic==jc)+W_off(i,l,5+k))
              if (ifD/=0_8) WD(node,l,ic,jc)=0.5d0*(merge(W_off(i,l,2),0.0d0,ic==jc)-W_off(i,l,14+k))
              if (ifSn/=0_8) WSn(node,l,ic,jc)=0.5d0*(merge(tnx_W_Sn,0.0d0,ic==jc)-W_off(i,l,23+k))
              if (ifDn/=0_8) WDn(node,l,ic,jc)=W_Dn(i,l,k)
            enddo
          enddo
        enddo
      enddo
      ier(off_rows(1:m_off))=0_8
    end block
  endif
  if (m_on>0_8 .and. (ifS/=0_8 .or. ifD/=0_8 .or. ifSn/=0_8 .or. ifDn/=0_8)) then
    block
      ! W_on(:,:,1:6): components 11,12,13,22,23,33 -> symmetric 3x3; W_Dn(:,:,1:9): (ic,jc) at ic+3*(jc-1).
      ! Allocated, since m_on*256*9 can exceed the stack.
      real(8) :: t0_on(2,m_on),h_on(m_on),u_on(m_on),v_on(m_on),P_u(m_on,0:15),P_v(m_on,0:15)
      real(8) :: r_u0(3),r_v0(3),N0(3),tnx_n_r0(m_on)
      real(8), allocatable :: W_on(:,:,:), W_Dn(:,:,:)
      integer(8) :: k,ic,jc,l
      allocate(W_on(m_on,256,6))
      t0_on=t0(:,on_rows(1:m_on))
      h_on=-h(on_rows(1:m_on))
      if (ifS/=0_8) then
        call quadS_stok_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        k=0_8
        do ic=1,3
          do jc=ic,3
            k=k+1_8
            WS(on_rows(1:m_on),:,ic,jc)=W_on(:,:,k)
            WS(on_rows(1:m_on),:,jc,ic)=W_on(:,:,k)
          enddo
        enddo
      endif
      if (ifD/=0_8) then
        call quadD_stok_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        k=0_8
        do ic=1,3
          do jc=ic,3
            k=k+1_8
            WD(on_rows(1:m_on),:,ic,jc)=W_on(:,:,k)
            WD(on_rows(1:m_on),:,jc,ic)=W_on(:,:,k)
          enddo
        enddo
      endif
      if (ifSn/=0_8 .or. ifDn/=0_8) then   ! the chart normal n_r0 = r_u0 x r_v0/|.| at t0; target normal tnx = +-n_r0
        u_on=t0_on(1,:); v_on=t0_on(2,:)
        call compute_legendre_values(15_8,u_on,p=P_u)
        call compute_legendre_values(15_8,v_on,p=P_v)
        do i=1,m_on
          do ic=1,3
            r_u0(ic)=dot_product(matmul(P_u(i,:),C_Xi(:,:,ic,1)),P_v(i,:))
            r_v0(ic)=dot_product(matmul(P_u(i,:),C_Xi(:,:,ic,2)),P_v(i,:))
          enddo
          N0=[r_u0(2)*r_v0(3)-r_u0(3)*r_v0(2),r_u0(3)*r_v0(1)-r_u0(1)*r_v0(3),r_u0(1)*r_v0(2)-r_u0(2)*r_v0(1)]
          tnx_n_r0(i)=dot_product(tnx(:,on_rows(i)),N0)/norm2(N0)
        enddo
      endif
      if (ifSn/=0_8) then
        call quadSn_stok_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_on)
        k=0_8
        do ic=1,3
          do jc=ic,3
            k=k+1_8
            do l=1,n
              WSn(on_rows(1:m_on),l,ic,jc)=tnx_n_r0*W_on(:,l,k)
              WSn(on_rows(1:m_on),l,jc,ic)=tnx_n_r0*W_on(:,l,k)
            enddo
          enddo
        enddo
      endif
      if (ifDn/=0_8) then
        allocate(W_Dn(m_on,256,9))
        call quadDn_stok_on_w16(m_on,t0_on,h_on,C_Xi,C_Xij,R,E,density_map,W_Dn)
        do jc=1,3
          do ic=1,3
            do l=1,n
              WDn(on_rows(1:m_on),l,ic,jc)=tnx_n_r0*W_Dn(:,l,ic+3*(jc-1))
            enddo
          enddo
        enddo
      endif
      ier(on_rows(1:m_on))=0_8
    end block
  endif
  ! A non-finite weight (e.g. J0=0 at t0) sets ier=6 for that row.
  if (m_near>0_8) then
    block
      integer(8) :: ic,jc
      do jc=1,3*ifS
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.ieee_is_finite(WS(near_rows(i),j,ic,jc))) ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
      do jc=1,3*ifD
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.ieee_is_finite(WD(near_rows(i),j,ic,jc))) ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
      do jc=1,3*ifSn
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.ieee_is_finite(WSn(near_rows(i),j,ic,jc))) ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
      do jc=1,3*ifDn
        do ic=1,3
          do j=1,n
            do i=1,m_near
              if (.not.ieee_is_finite(WDn(near_rows(i),j,ic,jc))) ier(near_rows(i))=6_8
            enddo
          enddo
        enddo
      enddo
    end block
  endif

  ! 5. Mid: quadSD_stok_mid_w16, the Stokes kernels at the integration nodes of each region.
  if (m_mid>0_8 .and. (ifS/=0_8 .or. ifD/=0_8 .or. ifSn/=0_8 .or. ifDn/=0_8)) then
    block
      real(8) :: tx_mid(3,m_mid),tnx_mid(3,m_mid),sx_mid(3,n),t0_mid(2,m_mid),h_mid(m_mid)
      real(8) :: geometry_center(2),geometry_scale(2)
      real(8) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3),geometry_center_typeb(2),geometry_scale_typeb(2)
      real(8), allocatable :: WS_mid(:,:,:,:),WD_mid(:,:,:,:),WSn_mid(:,:,:,:),WDn_mid(:,:,:,:)   ! can exceed the stack
      integer(8) :: classification_mid(m_mid),ier_mid(m_mid),l
      tx_mid=tx(:,mid_rows(1:m_mid))-spread(center,2,m_mid)
      sx_mid=sx-spread(center,2,n)
      tnx_mid=tnx(:,mid_rows(1:m_mid))
      t0_mid=t0(:,mid_rows(1:m_mid))
      h_mid=h(mid_rows(1:m_mid))
      if (onsurface_comp/=0_8) then
        where (other(mid_rows(1:m_mid))==0_8 .and. .not.ieee_is_finite(h_mid)) h_mid=0.0d0
      endif
      classification_mid=0_8; ier_mid=0_8
      nan=ieee_value(0.0d0,ieee_quiet_nan)
      allocate(WS_mid(m_mid,n*ifS,3,3),WD_mid(m_mid,n*ifD,3,3),WSn_mid(m_mid,n*ifSn,3,3),WDn_mid(m_mid,n*ifDn,3,3))
      WS_mid=nan; WD_mid=nan; WSn_mid=nan; WDn_mid=nan
      geometry_center=0.0d0; geometry_scale=1.0d0
      C_Xi_typeb=0.0d0; C_Xij_typeb=0.0d0
      geometry_center_typeb=0.0d0; geometry_scale_typeb=0.0d0
      if (ib/=0_8) then
        C_Xi_typeb(:,:,:,1)=cxpb(:,:,:,1,ib)/tsca(1)
        C_Xi_typeb(:,:,:,2)=cxpb(:,:,:,2,ib)/tsca(2)
        C_Xij_typeb(:,:,:,1)=cxppb(:,:,:,1,ib)/tsca(1)**2
        C_Xij_typeb(:,:,:,2)=cxppb(:,:,:,2,ib)/(tsca(1)*tsca(2))
        C_Xij_typeb(:,:,:,3)=cxppb(:,:,:,3,ib)/tsca(2)**2
        geometry_center_typeb=tcen
        geometry_scale_typeb=tsca
      endif
      call quadSD_stok_mid_w16(m_mid,tx_mid,tnx_mid,n,sx_mid,snx,sw,WS_mid,WD_mid,WSn_mid,WDn_mid,ifS,ifD,ifSn,ifDn, &
          t0_mid,h_mid,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification_mid,C_Xi_typeb,C_Xij_typeb, &
          geometry_center_typeb,geometry_scale_typeb,ier_mid,root_ok(mid_rows(1:m_mid)))
      where (ier_mid==1_8) ier_mid=2_8
      do l=1,n
        do i=1,m_mid
          if (ier_mid(i)/=0_8) cycle
          if (ifS/=0_8) then
            if (.not.all(ieee_is_finite(WS_mid(i,l,:,:)))) ier_mid(i)=6_8
          endif
          if (ifD/=0_8) then
            if (.not.all(ieee_is_finite(WD_mid(i,l,:,:)))) ier_mid(i)=6_8
          endif
          if (ifSn/=0_8) then
            if (.not.all(ieee_is_finite(WSn_mid(i,l,:,:)))) ier_mid(i)=6_8
          endif
          if (ifDn/=0_8) then
            if (.not.all(ieee_is_finite(WDn_mid(i,l,:,:)))) ier_mid(i)=6_8
          endif
        enddo
      enddo
      ier(mid_rows(1:m_mid))=ier_mid
      do i=1,m_mid
        if (ier_mid(i)/=0_8) cycle
        node=mid_rows(i)
        if (ifS/=0_8) WS(node,:,:,:)=WS_mid(i,:,:,:)
        if (ifD/=0_8) WD(node,:,:,:)=WD_mid(i,:,:,:)
        if (ifSn/=0_8) WSn(node,:,:,:)=WSn_mid(i,:,:,:)
        if (ifDn/=0_8) WDn(node,:,:,:)=WDn_mid(i,:,:,:)
      enddo
    end block
  endif

  ! 6. Unknown rows retained
  nan=ieee_value(0.0d0,ieee_quiet_nan)
  do i=1,m
    if (root_ok(i)) cycle
    t0(:,i)=nan
    h(i)=nan
  enddo
end subroutine quadSD_stok_w16

subroutine prepare_tables_w16(g16,R,E,density_map)
  ! Radial and edge quadrature tables for the near/self drivers.
  ! density_map = (P(g16)^T)^-1: GL16 node values from Legendre coefficients (step 10).
  ! R: GL24 on [0,1], 1/lambda, T1 (c1r = T1*V1r), T2 (cr = T2*(d.Gr)); off surface (step 8): Lj, the SVD
  ! U_V, S_inv, W_V of the Vandermonde matrix, ell0, ell0_coeff, ell0_over_lambda, radial_inv_sums,
  ! ell0_inv_sums, binom, c_y.
  ! E: Gauss rules of order 8r, r=1..32, on [-1,1]; budget 40.
  use surfquadr, only: radial_tables, edge_tables, compute_legendre_values, compute_gauss_legendre
  implicit none
  real(8), intent(in) :: g16(16)
  type(radial_tables), intent(inout) :: R
  type(edge_tables), intent(inout) :: E
  real(8), intent(inout) :: density_map(16,16)
  integer(8) :: i,j,k,m,ip,irule,n,info,piv(16)
  real(8) :: P(16,16),bary(24),L(24,24),V(24,24),VT(24,24),sigma(24),lapack_work(8192),excluded(24,3),s
  external :: dgesv,dgesvd

  call compute_legendre_values(15_8,g16,p=P)
  P=transpose(P)
  density_map=0.0d0
  do i=1,16
    density_map(i,i)=1.0d0
  enddo
  call dgesv(16_8,16_8,P,16_8,piv,density_map,16_8,info)

  call compute_gauss_legendre(24_8,R%lambda,R%w)
  R%lambda=(R%lambda+1.0d0)/2.0d0
  R%w=R%w/2.0d0
  R%lambda_inv=1.0d0/R%lambda
  bary=1.0d0
  do j=1,24
    do i=1,24
      if (i/=j) bary(i)=bary(i)*(R%lambda(i)-R%lambda(j))
    enddo
  enddo
  bary=1.0d0/bary
  ! T1(i,:): weights w*lambda against the interpolant at lambda(i)*lambda
  do i=1,24
    do j=1,24
      L(:,j)=bary(j)/(R%lambda(i)*R%lambda-R%lambda(j))
    enddo
    do k=1,24
      L(k,:)=L(k,:)/sum(L(k,:))
    enddo
    do j=1,24
      R%T1(i,j)=sum((R%w*R%lambda)*L(:,j))
    enddo
  enddo
  ! T2(i,:): weights w*lambda^2, same interpolant
  do i=1,24
    do j=1,24
      L(j,:)=bary/(R%lambda(i)*R%lambda(j)-R%lambda)
      L(j,:)=L(j,:)/sum(L(j,:))
    enddo
    R%T2(i,:)=matmul(R%w*R%lambda**2,L)
  enddo
  ! Off surface (step 8)
  R%ell0=bary/(-R%lambda)
  R%ell0=R%ell0/sum(R%ell0)
  do j=0,35
    R%Lj(:,j+1)=R%lambda**real(j,8)
  enddo
  V=transpose(R%Lj(:,1:24))
  call dgesvd('A','A',24_8,24_8,V,24_8,sigma,R%U_V,24_8,VT,24_8,lapack_work,8192_8,info)
  R%S_inv=1.0d0/sigma
  R%W_V=transpose(VT)
  excluded(:,1)=sum(R%lambda_inv)-R%lambda_inv   ! excluded(i,k) = sum_{j/=i} lambda_j^-k
  excluded(:,2)=sum(R%lambda_inv**2)-R%lambda_inv**2
  excluded(:,3)=sum(R%lambda_inv**3)-R%lambda_inv**3
  R%ell0_coeff(:,1)=R%ell0
  R%ell0_coeff(:,2)=-R%ell0*excluded(:,1)
  R%ell0_coeff(:,3)=R%ell0*(excluded(:,1)**2-excluded(:,2))/2.0d0
  R%ell0_coeff(:,4)=-R%ell0*(excluded(:,1)**3-3.0d0*excluded(:,1)*excluded(:,2)+2.0d0*excluded(:,3))/6.0d0
  R%ell0_over_lambda=R%ell0*R%lambda_inv
  do j=1,6
    R%radial_inv_sums(j)=sum(R%w*R%lambda_inv**j)
  enddo
  do j=1,4
    do k=1,2
      R%ell0_inv_sums(k,j)=sum(R%lambda_inv**k*R%ell0_coeff(:,j))
    enddo
  enddo
  R%binom=0.0d0
  do j=0,6
    s=1.0d0
    R%binom(1,j+1)=s
    do i=1,j
      s=s*(j-i+1)/i
      R%binom(i+1,j+1)=s
    enddo
  enddo
  R%c_y(1,:)=1.0d0
  do ip=0,3   ! p = ip (P is the Legendre matrix above)
    do m=1,6
      R%c_y(m+1,ip+1)=R%c_y(m,ip+1)*(-(real(ip,8)+0.5d0)-real(m,8)+1.0d0)/real(m,8)
    enddo
  enddo

  E%budget=40.0d0
  E%nodes=0.0d0
  E%weights=0.0d0
  do irule=1,32
    n=8_8*irule
    E%order(irule)=n
    call compute_gauss_legendre(n,E%nodes(1:n,irule),E%weights(1:n,irule))
  enddo
end subroutine prepare_tables_w16

subroutine project_t0_h_w16(m,n,nb,nr,sx,st,cols,C_X,C_Xi,C_Xij,tx,target_idx,other,onsurface_comp,ib, &
    xb,g16,cols_typeb,cxb,cxpb,cxppb,tcen,tsca,tau0,root_x,root_t,root_bound,if_h0,t0,h,root_ok)
  use, intrinsic :: iso_c_binding, only: c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  implicit none
  integer(8), intent(in) :: m,n,nb,nr,cols(n),target_idx(m),other(m),onsurface_comp,ib
  integer(8), intent(in) :: cols_typeb(2304,2),if_h0(m)
  real(8), intent(in) :: sx(3,n),st(2,n),C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3)
  real(8), intent(in) :: tx(3,m),xb(3,nb),g16(16),cxb(16,16,3,2)
  real(8), intent(in) :: cxpb(16,16,3,2,2),cxppb(16,16,3,3,2),tcen(2),tsca(2)
  real(8), intent(in) :: tau0(2,m)
  real(8), intent(in) :: root_x(3,nr),root_t(2,nr),root_bound(2)
  real(8), intent(inout) :: t0(2,m),h(m)
  logical(c_bool), intent(inout) :: root_ok(m)
  ! Complete finite pairs
  integer(8) :: itarg,node,seed,nroot,iu,iv,inode,jj,branch,degree,coordinate,iteration,info,piv(2)
  integer(8) :: pair_index(m),branches(2)
  real(8) :: tb1(48),tau1(48,2),tau2(48),C_root(16,16,3),C_root_i(16,16,3,2),C_root_ij(16,16,3,3)
  real(8) :: root_limit(2),root_center(2),root_scale(2),anchor(2),tt(2),tau(2),ends(2,2),edge(2)
  real(8) :: Pa(2,16),P(2,16),Pedge(2,16),Pe(2,16),D(2,16),De(2,16),deltae(2),delta(2)
  real(8) :: edge_to_end(3),edge_to_anchor(3),edge_to_target(3),y(3),X(3),J(3,2),D2(3,3)
  real(8) :: r(3),g(2),Hessian(2,2),step(2),trial(2),normal(3),distance2,best_distance2,nan
  real(8) :: step_size,residual_norm,gradient_tol(2),normal_length,h_try,tol_x
  real(8), allocatable :: root_x_local(:,:),root_t_local(:,:)
  logical :: converged,known(m)
  external :: dgesv

  if (m==0_8) return
  do itarg=1,m
    root_ok(itarg)=root_ok(itarg) .and. all(ieee_is_finite(t0(:,itarg))) .and. ieee_is_finite(h(itarg))
    known(itarg)=root_ok(itarg)
  enddo
  if (onsurface_comp/=0_8) then
    do itarg=1,m
      if (known(itarg) .or. target_idx(itarg)<=0_8) cycle
      do node=1,n
        if (cols(node)/=target_idx(itarg)) cycle
        t0(:,itarg)=st(:,node)
        h(itarg)=0.0d0
        root_ok(itarg)=.true.
        known(itarg)=.true.
        exit
      enddo
    enddo
  endif
  if (all(known)) return

  pair_index=0_8
  if (ib==0_8) then
    nroot=nr
    allocate(root_x_local(3,nroot),root_t_local(2,nroot))
    root_x_local=root_x; root_t_local=root_t
    C_root=C_X; C_root_i=C_Xi; C_root_ij=C_Xij
    root_limit=root_bound; root_center=0.0d0; root_scale=1.0d0
  else
    nroot=2304_8
    allocate(root_x_local(3,nroot),root_t_local(2,nroot))
    tb1(1:16)=g16/4.0d0-0.75d0
    tb1(17:32)=g16/4.0d0-0.25d0
    tb1(33:48)=g16/2.0d0+0.5d0
    tau1(:,1)=-tb1(48:1:-1); tau1(:,2)=tb1
    tau2(1:16)=(g16-2.0d0)/3.0d0
    tau2(17:32)=g16/3.0d0
    tau2(33:48)=(g16+2.0d0)/3.0d0
    do iv=1,48
      do iu=1,48
        node=iu+48_8*(iv-1_8)
        root_x_local(:,node)=xb(:,(ib-1_8)*48_8+iu+96_8*(iv-1_8))
        root_t_local(:,node)=[tau1(iu,ib),tau2(iv)]
      enddo
    enddo
    C_root=cxb(:,:,:,ib); C_root_i=cxpb(:,:,:,:,ib); C_root_ij=cxppb(:,:,:,:,ib)
    root_limit=1.0d0; root_center=tcen; root_scale=tsca
    if (onsurface_comp/=0_8) then
      do itarg=1,m
        if (other(itarg)==0_8 .or. target_idx(itarg)<=0_8) cycle
        do node=1,nroot
          if (cols_typeb(node,3_8-ib)/=target_idx(itarg)) cycle
          pair_index(itarg)=node
          exit
        enddo
      enddo
    endif
  endif

  nan=ieee_value(0.0d0,ieee_quiet_nan)
  do itarg=1,m
    if (known(itarg)) cycle
    seed=1_8; best_distance2=huge(1.0d0)
    do node=1,nroot
      distance2=sum((root_x_local(:,node)-tx(:,itarg))**2)
      if (distance2>=best_distance2) cycle
      seed=node; best_distance2=distance2
    enddo
    anchor=root_t_local(:,seed); tt=anchor
    Pa(:,1)=1.0d0; Pa(:,2)=anchor
    do degree=2,15
      Pa(:,degree+1)=((2*degree-1)*anchor*Pa(:,degree)-(degree-1)*Pa(:,degree-1))/degree
    enddo

    tau=tau0(:,itarg)
    if (pair_index(itarg)>0_8) then
      inode=pair_index(itarg)
      tau=[tau1(1_8+mod(inode-1_8,48_8),3_8-ib),tau2(1_8+(inode-1_8)/48_8)]
    endif

    if (ib>0_8 .and. onsurface_comp/=0_8 .and. other(itarg)/=0_8 .and. &
        all(ieee_is_finite(tau)) .and. all(abs(tau)<=1.0d0)) then
      branches=[ib,3_8-ib]; ends(:,1)=anchor; ends(:,2)=tau
      do jj=1,2
        branch=branches(jj)
        edge=[real(3_8-2_8*branch,8),tau(2)]
        Pedge(:,1)=1.0d0; Pedge(:,2)=edge
        Pe(:,1)=1.0d0; Pe(:,2)=ends(:,jj)
        De=0.0d0; De(:,2)=1.0d0
        do degree=2,15
          Pedge(:,degree+1)=((2*degree-1)*edge*Pedge(:,degree)-(degree-1)*Pedge(:,degree-1))/degree
          Pe(:,degree+1)=((2*degree-1)*ends(:,jj)*Pe(:,degree)-(degree-1)*Pe(:,degree-1))/degree
          De(:,degree+1)=((2*degree-1)*(Pe(:,degree)+edge*De(:,degree))-(degree-1)*De(:,degree-1))/degree
        enddo
        deltae=ends(:,jj)-edge
        do coordinate=1,3
          edge_to_end(coordinate)= &
            deltae(1)*dot_product(matmul(De(1,:),cxb(:,:,coordinate,branch)),Pe(2,:)) &
            +deltae(2)*dot_product(matmul(Pedge(1,:),cxb(:,:,coordinate,branch)),De(2,:))
        enddo
        if (jj==1_8) then
          edge_to_anchor=edge_to_end
        else
          edge_to_target=edge_to_end
        endif
      enddo
      y=edge_to_target-edge_to_anchor
    else
      y=tx(:,itarg)-root_x_local(:,seed)
    endif

    converged=.false.
    do iteration=1,30
      P(:,1)=1.0d0; P(:,2)=tt
      D=0.0d0; D(:,2)=1.0d0
      do degree=2,15
        P(:,degree+1)=((2*degree-1)*tt*P(:,degree)-(degree-1)*P(:,degree-1))/degree
        D(:,degree+1)=((2*degree-1)*(P(:,degree)+anchor*D(:,degree))-(degree-1)*D(:,degree-1))/degree
      enddo
      delta=tt-anchor
      do coordinate=1,3
        X(coordinate)=delta(1)*dot_product(matmul(D(1,:),C_root(:,:,coordinate)),P(2,:)) &
          +delta(2)*dot_product(matmul(Pa(1,:),C_root(:,:,coordinate)),D(2,:))
        J(coordinate,1)=dot_product(matmul(P(1,:),C_root_i(:,:,coordinate,1)),P(2,:))
        J(coordinate,2)=dot_product(matmul(P(1,:),C_root_i(:,:,coordinate,2)),P(2,:))
        D2(coordinate,1)=dot_product(matmul(P(1,:),C_root_ij(:,:,coordinate,1)),P(2,:))
        D2(coordinate,2)=dot_product(matmul(P(1,:),C_root_ij(:,:,coordinate,2)),P(2,:))
        D2(coordinate,3)=dot_product(matmul(P(1,:),C_root_ij(:,:,coordinate,3)),P(2,:))
      enddo
      r=X-y
      if (.not.all(ieee_is_finite(r)) .or. .not.all(ieee_is_finite(J)) .or. &
          .not.all(ieee_is_finite(D2))) exit
      g=matmul(transpose(J),r)
      Hessian=matmul(transpose(J),J)
      Hessian(1,1)=Hessian(1,1)+dot_product(r,D2(:,1))
      Hessian(1,2)=Hessian(1,2)+dot_product(r,D2(:,2))
      Hessian(2,1)=Hessian(2,1)+dot_product(r,D2(:,2))
      Hessian(2,2)=Hessian(2,2)+dot_product(r,D2(:,3))
      step=g
      call dgesv(2_8,1_8,Hessian,2_8,piv,step,2_8,info)
      if (info/=0_8 .or. .not.all(ieee_is_finite(step))) exit
      step_size=maxval(abs(step))
      residual_norm=norm2(r)
      gradient_tol=2.0d0*epsilon(1.0d0)*[norm2(J(:,1)),norm2(J(:,2))]*residual_norm
      if (step_size<=epsilon(1.0d0)*max(1.0d0,maxval(abs(tt))) .or. all(abs(g)<=gradient_tol)) then
        converged=.true.
        exit
      endif
      trial=max(-root_limit,min(root_limit,tt-step/max(1.0d0,step_size)))
      if (all(trial==tt) .or. iteration==30_8) exit
      tt=trial
    enddo

    normal=[J(2,2)*J(3,1)-J(3,2)*J(2,1),J(3,2)*J(1,1)-J(1,2)*J(3,1), &
            J(1,2)*J(2,1)-J(2,2)*J(1,1)]
    normal_length=norm2(normal)
    if (.not.ieee_is_finite(normal_length) .or. normal_length==0.0d0 .or. .not.all(ieee_is_finite(r))) then
      t0(:,itarg)=nan
      h(itarg)=nan
      cycle
    endif
    normal=normal/normal_length
    if (if_h0(itarg)/=0_8) then
      h_try=0.0d0
    else
      h_try=-dot_product(normal,r)
    endif
    tol_x=64.0d0*spacing(max(1.0d0,maxval(abs(tx(:,itarg)))))
    if (all(abs(tt)<=root_limit) .and. ieee_is_finite(h_try)) then
      t0(:,itarg)=root_center+root_scale*tt
      h(itarg)=h_try
      root_ok(itarg)=converged .and. maxval(abs(r+h_try*normal))<=tol_x
      if (.not.root_ok(itarg) .and. if_h0(itarg)/=0_8) h(itarg)=-dot_product(normal,r)
    else
      t0(:,itarg)=nan
      h(itarg)=nan
    endif
  enddo
end subroutine project_t0_h_w16

subroutine prepare_quadSD_targets_w16(m,tx,n,sx,C_X,C_Xi,geometry_center,geometry_scale, &
    inside16,outside16,xc,rt0,rp0,e1,e2,e3,r,xyz,d16,d128,lambda1)
  use surfquadr_mid, only: legendrevec
  use, intrinsic :: iso_c_binding, only: c_bool
  implicit none
  integer(8), intent(in) :: m,n
  real(8), intent(in) :: tx(3,m),sx(3,n),C_X(16,16,3),C_Xi(16,16,3,2)
  real(8), intent(in) :: geometry_center(2),geometry_scale(2)
  logical(c_bool), intent(inout) :: inside16(m),outside16(m)
  real(8), intent(inout) :: xc(3),rt0(3),rp0(3),e1(3),e2(3),e3(3),r(3),xyz(3,m)
  real(8), intent(inout) :: d16,d128,lambda1
  integer(8) :: k
  real(8) :: P(2,16),frame(3,3),q(3),x_near(3,9),abc16(3)
  real(8), allocatable :: Pu(:,:),Pv(:,:),X(:,:)

  call legendrevec(15_8,-geometry_center/geometry_scale,P)
  do k=1,3
    xc(k)=dot_product(matmul(P(1,:),C_X(:,:,k)),P(2,:))
    rt0(k)=dot_product(matmul(P(1,:),C_Xi(:,:,k,1)),P(2,:))
    rp0(k)=dot_product(matmul(P(1,:),C_Xi(:,:,k,2)),P(2,:))
  enddo
  e1=rp0/norm2(rp0)
  e2=rt0-e1*dot_product(e1,rt0); e2=e2/norm2(e2)
  e3=[e1(2)*e2(3)-e1(3)*e2(2),e1(3)*e2(1)-e1(1)*e2(3),e1(1)*e2(2)-e1(2)*e2(1)]
  frame(:,1)=e1; frame(:,2)=e2; frame(:,3)=e3

  ! Extrapolate geometry to the four corners, edge midpoints and center.
  q=[-1.1d0,0.0d0,1.1d0]
  allocate(Pu(3,16),Pv(3,16),X(3,3))
  call legendrevec(15_8,(q-geometry_center(1))/geometry_scale(1),Pu)
  call legendrevec(15_8,(q-geometry_center(2))/geometry_scale(2),Pv)
  do k=1,3
    X=matmul(matmul(Pu,C_X(:,:,k)),transpose(Pv))
    x_near(k,:)=reshape(X,[9])
  enddo
  deallocate(Pu,Pv,X)
  r=max(maxval(abs(matmul(transpose(frame),sx-spread(xc,2,n))),dim=2), &
        maxval(abs(matmul(transpose(frame),x_near-spread(xc,2,9))),dim=2))
  xyz=matmul(transpose(frame),tx-spread(xc,2,m))
  d16=1.25d0*maxval(r(1:2)); d128=0.20d0*d16
  ! One common factor keeps the layers nested and encloses the sampled bounds.
  lambda1=(1.0d0+1.0d-12)*max(1.0d0,sqrt(sum(r*r/(r*r+d128*d128))))
  abc16=lambda1*sqrt(r*r+d16*d16)
  outside16=sum((xyz/spread(abc16,2,m))**2,dim=1)>1.0d0
  inside16=.not.outside16
end subroutine prepare_quadSD_targets_w16
