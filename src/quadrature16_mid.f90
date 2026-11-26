module surfquadr_mid
  ! Legendre evaluation and Gauss rules for mid-range quadrature.
  ! The mid-range routines use -ffp-contract=fast.
  implicit none
  private
  public :: legendrevec, compute_gaussleg

contains

  subroutine legendrevec(n,x,P)
    ! Input: degree n, points x(:). Inout: P(size(x),n+1), P(:,l+1)=P_l(x).
    integer(8), intent(in) :: n
    real(8), intent(in) :: x(:)
    real(8), intent(inout) :: P(:,:)
    integer(8) :: l

    P(:,1)=1; P(:,2)=x
    do l=1,n-1
      P(:,l+2)=((2*l+1)*x*P(:,l+1)-l*P(:,l))/(l+1)
    end do
  end subroutine legendrevec

  subroutine compute_gaussleg(order,nodes,weights)
    ! Gauss-Legendre rule on [-1,1]. Input: order. Inout: nodes(order), weights(order).
    ! Newton iteration finds the roots of P_order on one half; symmetry supplies the other.
    ! Christoffel weights: weights(m)=2/sum_{ell=0}^{order-1}(2*ell+1)*P_ell(nodes(m))^2.
    integer(8), intent(in) :: order
    real(8), intent(inout) :: nodes(order), weights(order)
    real(8) :: root((order+1)/2), previous((order+1)/2), current((order+1)/2), next((order+1)/2)
    real(8) :: derivative((order+1)/2), candidate((order+1)/2), total((order+1)/2), pi
    integer(8) :: half, k, iteration, degree

    half=(order+1)/2; pi=acos(-1.0d0)
    do k=1,half
      root(k)=cos(pi*(k-0.25d0)/(order+0.5d0))
    end do
    if (mod(order,2_8)==1) root(half)=0
    do iteration=1,10
      previous=1; current=root
      do degree=1,order-1
        next=((2*degree+1)*root*current-degree*previous)/(degree+1)
        previous=current; current=next
      end do
      derivative=order*(previous-root*current)/((1-root)*(1+root))
      candidate=root-current/derivative
      if (all(candidate==root)) exit
      root=candidate
    end do
    previous=1; current=root; total=1
    do degree=1,order-1
      total=total+(2*degree+1)*current*current
      next=((2*degree+1)*root*current-degree*previous)/(degree+1)
      previous=current; current=next
    end do
    nodes(1:half)=-root; weights(1:half)=2/total
    do k=1,order/2
      nodes(half+k)=root(order/2-k+1); weights(half+k)=weights(order/2-k+1)
    end do
  end subroutine compute_gaussleg

end module surfquadr_mid

! SLP, DLP and the requested SLPn and DLPn components share quadrature nodes within each region.
subroutine quadSD_mid_w16(m,tx,tnx,n,sx,snx,sw,WS,WD,ifSn,ifDn,WSn,WDn, &
    t0,h,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification, &
    C_Xi_typeb,C_Xij_typeb,geometry_center_typeb,geometry_scale_typeb,ier,root_ok)
  use surfquadr_mid, only: legendrevec, gauss => compute_gaussleg
  use, intrinsic :: iso_c_binding, only: c_double,c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  implicit none
  integer(8), intent(in) :: m,n
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),t0(2,m),h(m)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3)
  real(8), intent(in) :: geometry_center(2),geometry_scale(2)
  ! *_typeb: geometry of the type-b chart the source sits in; geometry_scale_typeb = 0: no type-b chart.
  real(8), intent(in) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3)
  real(8), intent(in) :: geometry_center_typeb(2),geometry_scale_typeb(2)
  integer(8), intent(in) :: ifSn,ifDn
  real(8), intent(inout) :: WS(m,n),WD(m,n),WSn(m,n*ifSn,3),WDn(m,n*ifDn,3,3)
  integer(8), intent(inout) :: classification(m),ier(m)
  logical(c_bool), intent(in) :: root_ok(m)
  real(8) :: t0_local(2,m),h_local(m)
  ! The 256 source nodes are ordered with u varying fastest.
  ! The caller excludes near/self and supplies the Region 7 projection and signed height.
  ! Regions 1:3 use 1D parity factors for every kernel.
  ! At a node with d = y-x, r = |d| and SLP value v, SLPn_c = (v/r^2)*d_c and
  ! DLPn_ab = (v/r^2)*(3*d_a*d_b/r^2-delta_ab), DLPn components (a,b) = (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
  ! Where the kernels share one contraction (regions 1:6) they are stacked by kernel index k: 1 SLP, 2 DLP,
  ! k_Sn+kn SLPn kn, k_Dn+kn DLPn kn; n_kernel = 2, 5, 8 or 11.
  integer(8), parameter :: i_Dn(6)=[1_8,1_8,1_8,2_8,2_8,3_8], j_Dn(6)=[1_8,2_8,3_8,2_8,3_8,3_8]
  integer(8) :: n_kernel,k_Sn,k_Dn,kn,a_c,b_c,q_node
  real(8) :: delta_c
  ! SLPn, DLPn: displacement y-x, 1/r^2 and v/r^2 at the nodes of one evaluation step
  real(8), allocatable :: node_d(:,:,:),node_r2(:,:),node_f(:,:),W_k(:,:)
  integer(8) :: i,j,k,m16,mr,nq,nh,ld,layer,info,piv(16)
  ! integer(8), parameter :: order_p(3)=[32,64,128]
  integer(8), parameter :: order_p(3)=[32,48,96]
  integer(8), allocatable :: idx16(:),rows(:)
  real(8) :: g16(16),wg16(16),V16(16,16),Vsolve(16,16),density_map(16,16),source_t(2,256)
  real(8) :: xc(3),rt0(3),rp0(3),e1(3),e2(3),e3(3),frame(3,3),r(3),xyz(3,m)
  real(8) :: ratio
  real(8) :: lambda1,d16,d32,d64,d128,abc32(3),abc64(3),abc128(3)
  real(8), parameter :: pi=acos(-1.0d0),inv2pi=0.5d0/pi
  logical(c_bool) :: outside16(m),inside16(m)
  logical :: inside32(m),inside64(m),inside128(m),shell32(m),shell64(m),shell128(m)
  logical :: target_groups(m,3),remaining(m),use_sinh(m),off_column4(m)
  real(8), allocatable :: Pu(:,:),Pv(:,:),X(:,:),rt(:,:,:),rp(:,:,:),Xtmp(:,:)
  real(8), allocatable :: g(:),wg(:),L(:,:,:),sq_x(:,:),sq_nx(:,:),sq_w(:),target(:,:),Y(:,:)
  real(8) :: P_half(maxval(order_p)/2,16),L_half(maxval(order_p)/2,16)
  real(8) :: dx,dy,dz,distance2_inv,distance_inv,kernel_value
  integer(8) :: iu,iv,mode,iq,is,node(4,4096),column(4,64)
  real(8) :: v1,v2,v3,v4,p1,p2,p3,p4,d1,d2,d3,d4
  real(8), allocatable :: parity_K(:,:,:),parity_W(:,:,:)
  external :: dgemm,dgesv

  ier=0_8
  n_kernel=2_8+3_8*ifSn+6_8*ifDn; k_Sn=2_8; k_Dn=2_8+3_8*ifSn
  call prepare_quadSD_targets_w16(m,tx,n,sx,C_X,C_Xi,geometry_center,geometry_scale, &
      inside16,outside16,xc,rt0,rp0,e1,e2,e3,r,xyz,d16,d128,lambda1)
  frame(:,1)=e1; frame(:,2)=e2; frame(:,3)=e3

  d32=0.50d0*d16; d64=0.25d0*d16
  abc32=lambda1*sqrt(r*r+d32*d32)
  abc64=lambda1*sqrt(r*r+d64*d64); abc128=lambda1*sqrt(r*r+d128*d128)

  ! 1st layer of guard.
  where (outside16) classification=0_8
  m16=count(outside16); idx16=pack([(k,k=1,m)],outside16)
  if (m16>0) then
    target=tx(:,idx16)
    do j=1,n
      !GCC$ ivdep
      do i=1,m16
        dx=sx(1,j)-target(1,i); dy=sx(2,j)-target(2,i); dz=sx(3,j)-target(3,i)
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz); distance_inv=sqrt(distance2_inv)
        kernel_value=(sw(j)*inv2pi)*distance_inv
        WS(idx16(i),j)=kernel_value
        WD(idx16(i),j)=kernel_value*(dx*snx(1,j)+dy*snx(2,j)+dz*snx(3,j))*distance2_inv
      enddo
    enddo
    if (n_kernel>2_8) then
      allocate(node_d(m16,3,1),node_r2(m16,1),node_f(m16,1))
      do j=1,n
        !GCC$ ivdep
        do i=1,m16
          dx=sx(1,j)-target(1,i); dy=sx(2,j)-target(2,i); dz=sx(3,j)-target(3,i)
          distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
          node_d(i,1,1)=dx; node_d(i,2,1)=dy; node_d(i,3,1)=dz; node_r2(i,1)=distance2_inv
          node_f(i,1)=((sw(j)*inv2pi)*sqrt(distance2_inv))*distance2_inv
        enddo
        do kn=1,3*ifSn
          do i=1,m16
            WSn(idx16(i),j,kn)=node_f(i,1)*node_d(i,kn,1)
          enddo
        enddo
        do kn=1,6*ifDn
          a_c=i_Dn(kn); b_c=j_Dn(kn); delta_c=merge(1.0d0,0.0d0,a_c==b_c)
          do i=1,m16
            v1=node_f(i,1)*(3.0d0*node_d(i,a_c,1)*node_d(i,b_c,1)*node_r2(i,1)-delta_c)
            WDn(idx16(i),j,a_c,b_c)=v1; WDn(idx16(i),j,b_c,a_c)=v1
          enddo
        enddo
      enddo
      deallocate(node_d,node_r2,node_f)
    endif
    deallocate(target)
  endif
  deallocate(idx16)

  inside32=inside16 .and. sum((xyz/spread(abc32,2,m))**2,dim=1)<=1.0d0
  shell32=inside16 .and. .not.inside32
  inside64=inside32 .and. sum((xyz/spread(abc64,2,m))**2,dim=1)<=1.0d0
  shell64=inside32 .and. .not.inside64
  inside128=inside64 .and. sum((xyz/spread(abc128,2,m))**2,dim=1)<=1.0d0
  shell128=inside64 .and. .not.inside128
  target_groups(:,1)=shell32; target_groups(:,2)=shell64; target_groups(:,3)=shell128

  if (any(inside16)) then
    call gauss(16_8,g16,wg16)
    call legendrevec(15_8,g16,V16)
    Vsolve=transpose(V16); density_map=0.0d0
    do j=1,16
      density_map(j,j)=1.0d0
    enddo
    call dgesv(16_8,16_8,Vsolve,16_8,piv,density_map,16_8,info)
    if (info/=0) return
  endif

  ! Common p refinement: contract v, then u for order_p.
  do layer=1,3
    where (target_groups(:,layer)) classification=layer
    mr=count(target_groups(:,layer))
    if (mr==0) cycle
    rows=pack([(k,k=1,m)],target_groups(:,layer))
    ! nq=16*2**layer; nh=nq/2
    nq=order_p(layer); nh=nq/2
    allocate(g(nq),wg(nq),L(nh,8,2))
    call gauss(nq,g,wg)
    call legendrevec(15_8,g(:nh),P_half(:nh,:))
    call dgemm('N','T',nh,16_8,16_8,1d0,P_half,size(P_half,1,kind=8), &
      density_map,16_8,0d0,L_half,size(L_half,1,kind=8))
    L(:,:,1)=0.5d0*(L_half(:nh,1:8)+L_half(:nh,16:9:-1))
    L(:,:,2)=0.5d0*(L_half(:nh,1:8)-L_half(:nh,16:9:-1))
    allocate(Pu(nq,16),Pv(nq,16),X(nq,nq),rt(nq,nq,3),rp(nq,nq,3),Xtmp(nq,16))
    allocate(sq_x(3,nq*nq),sq_nx(3,nq*nq),sq_w(nq*nq))
    call legendrevec(15_8,(g-geometry_center(1))/geometry_scale(1),Pu)
    call legendrevec(15_8,(g-geometry_center(2))/geometry_scale(2),Pv)
    do k=1,3
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_X(:,:,k),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,X,nq)
      sq_x(k,:)=reshape(X,[nq*nq])
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,1),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rt(:,:,k),nq)
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,2),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rp(:,:,k),nq)
    enddo
    sq_nx(1,:)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[nq*nq])
    sq_nx(2,:)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[nq*nq])
    sq_nx(3,:)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[nq*nq])
    sq_w=sqrt(sum(sq_nx*sq_nx,dim=1)); sq_nx=sq_nx/spread(sq_w,1,3)
    sq_w=sq_w*reshape(spread(wg,2,nq)*spread(wg,1,nq),[nq*nq])
    target=transpose(tx(:,rows))
    ! Keep targets contiguous; contract v, then u without transposing kernels.
    ! Row i+mr*(k-1): target i, kernel k.
    ld=8*((n_kernel*mr+7)/8)
    allocate(parity_K(ld,nh*nh,4),parity_W(ld,64,4),Y(ld*nh,8))
    parity_K(n_kernel*mr+1:ld,:,:)=0.0d0
    do iv=1,nh
      do iu=1,nh
        iq=iu+nh*(iv-1)
        node(:,iq)=[iu+nq*(iv-1),nq+1-iu+nq*(iv-1),iu+nq*(nq-iv),nq+1-iu+nq*(nq-iv)]
      enddo
    enddo
    do iv=1,8
      do iu=1,8
        is=iu+8*(iv-1)
        column(:,is)=[iu+16*(iv-1),17-iu+16*(iv-1),iu+16*(16-iv),17-iu+16*(16-iv)]
      enddo
    enddo
    sq_nx=sq_nx*spread(sq_w*inv2pi,1,3)
    do j=1,nh*nh
      !GCC$ ivdep
      !GCC$ unroll 1
      do i=1,mr
        dx=sq_x(1,node(1,j))-target(i,1); dy=sq_x(2,node(1,j))-target(i,2); dz=sq_x(3,node(1,j))-target(i,3)
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
        distance_inv=sqrt(distance2_inv); v1=(sq_w(node(1,j))*inv2pi)*distance_inv
        d1=((dx*sq_nx(1,node(1,j))+dy*sq_nx(2,node(1,j))+dz*sq_nx(3,node(1,j)))*distance2_inv)*distance_inv
        dx=sq_x(1,node(2,j))-target(i,1); dy=sq_x(2,node(2,j))-target(i,2); dz=sq_x(3,node(2,j))-target(i,3)
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
        distance_inv=sqrt(distance2_inv); v2=(sq_w(node(2,j))*inv2pi)*distance_inv
        d2=((dx*sq_nx(1,node(2,j))+dy*sq_nx(2,node(2,j))+dz*sq_nx(3,node(2,j)))*distance2_inv)*distance_inv
        parity_K(i,j,1)=v1+v2; parity_K(i,j,2)=v1-v2
        parity_K(mr+i,j,1)=d1+d2; parity_K(mr+i,j,2)=d1-d2
      enddo
      !GCC$ ivdep
      !GCC$ unroll 1
      do i=1,mr
        dx=sq_x(1,node(3,j))-target(i,1); dy=sq_x(2,node(3,j))-target(i,2); dz=sq_x(3,node(3,j))-target(i,3)
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
        distance_inv=sqrt(distance2_inv); v3=(sq_w(node(3,j))*inv2pi)*distance_inv
        d3=((dx*sq_nx(1,node(3,j))+dy*sq_nx(2,node(3,j))+dz*sq_nx(3,node(3,j)))*distance2_inv)*distance_inv
        dx=sq_x(1,node(4,j))-target(i,1); dy=sq_x(2,node(4,j))-target(i,2); dz=sq_x(3,node(4,j))-target(i,3)
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
        distance_inv=sqrt(distance2_inv); v4=(sq_w(node(4,j))*inv2pi)*distance_inv
        d4=((dx*sq_nx(1,node(4,j))+dy*sq_nx(2,node(4,j))+dz*sq_nx(3,node(4,j)))*distance2_inv)*distance_inv
        p1=v3+v4; p2=v3-v4; p3=parity_K(i,j,1); p4=parity_K(i,j,2)
        parity_K(i,j,1)=p3+p1; parity_K(i,j,2)=p4+p2
        parity_K(i,j,3)=p3-p1; parity_K(i,j,4)=p4-p2
        p1=d3+d4; p2=d3-d4; p3=parity_K(mr+i,j,1); p4=parity_K(mr+i,j,2)
        parity_K(mr+i,j,1)=p3+p1; parity_K(mr+i,j,2)=p4+p2
        parity_K(mr+i,j,3)=p3-p1; parity_K(mr+i,j,4)=p4-p2
      enddo
    enddo
    if (n_kernel>2_8) then   ! SLPn, DLPn rows: the same parity sums at the four mirrored nodes
      allocate(node_d(mr,3,4),node_r2(mr,4),node_f(mr,4))
      do j=1,nh*nh
        do q_node=1,4
          !GCC$ ivdep
          do i=1,mr
            dx=sq_x(1,node(q_node,j))-target(i,1); dy=sq_x(2,node(q_node,j))-target(i,2)
            dz=sq_x(3,node(q_node,j))-target(i,3)
            distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz)
            node_d(i,1,q_node)=dx; node_d(i,2,q_node)=dy; node_d(i,3,q_node)=dz; node_r2(i,q_node)=distance2_inv
            node_f(i,q_node)=((sq_w(node(q_node,j))*inv2pi)*sqrt(distance2_inv))*distance2_inv
          enddo
        enddo
        do kn=1,3*ifSn
          k=k_Sn+kn
          !GCC$ ivdep
          do i=1,mr
            v1=node_f(i,1)*node_d(i,kn,1); v2=node_f(i,2)*node_d(i,kn,2)
            v3=node_f(i,3)*node_d(i,kn,3); v4=node_f(i,4)*node_d(i,kn,4)
            p1=v1+v2; p2=v1-v2; p3=v3+v4; p4=v3-v4
            parity_K(i+mr*(k-1),j,1)=p1+p3; parity_K(i+mr*(k-1),j,2)=p2+p4
            parity_K(i+mr*(k-1),j,3)=p1-p3; parity_K(i+mr*(k-1),j,4)=p2-p4
          enddo
        enddo
        do kn=1,6*ifDn
          k=k_Dn+kn; a_c=i_Dn(kn); b_c=j_Dn(kn); delta_c=merge(1.0d0,0.0d0,a_c==b_c)
          !GCC$ ivdep
          do i=1,mr
            v1=node_f(i,1)*(3.0d0*node_d(i,a_c,1)*node_d(i,b_c,1)*node_r2(i,1)-delta_c)
            v2=node_f(i,2)*(3.0d0*node_d(i,a_c,2)*node_d(i,b_c,2)*node_r2(i,2)-delta_c)
            v3=node_f(i,3)*(3.0d0*node_d(i,a_c,3)*node_d(i,b_c,3)*node_r2(i,3)-delta_c)
            v4=node_f(i,4)*(3.0d0*node_d(i,a_c,4)*node_d(i,b_c,4)*node_r2(i,4)-delta_c)
            p1=v1+v2; p2=v1-v2; p3=v3+v4; p4=v3-v4
            parity_K(i+mr*(k-1),j,1)=p1+p3; parity_K(i+mr*(k-1),j,2)=p2+p4
            parity_K(i+mr*(k-1),j,3)=p1-p3; parity_K(i+mr*(k-1),j,4)=p2-p4
          enddo
        enddo
      enddo
      deallocate(node_d,node_r2,node_f)
    endif
    do mode=1,4
      call dgemm('N','N',ld*nh,8_8,nh,1.0d0,parity_K(1,1,mode),ld*nh, &
        L(1,1,1+(mode-1)/2),nh,0.0d0,Y,ld*nh)
      do iv=1,8
        call dgemm('N','N',n_kernel*mr,8_8,nh,1.0d0,Y(1,iv),ld,L(1,1,1+mod(mode-1,2)),nh, &
          0.0d0,parity_W(1,1+8*(iv-1),mode),ld)
      enddo
    enddo
    do j=1,64
      !GCC$ ivdep
      do i=1,mr
        v1=parity_W(i,j,1); v2=parity_W(i,j,2); v3=parity_W(i,j,3); v4=parity_W(i,j,4)
        p1=v1+v2; p2=v3+v4; p3=v1-v2; p4=v3-v4
        WS(rows(i),column(1,j))=p1+p2; WS(rows(i),column(2,j))=p3+p4
        WS(rows(i),column(3,j))=p1-p2; WS(rows(i),column(4,j))=p3-p4
        v1=parity_W(mr+i,j,1); v2=parity_W(mr+i,j,2)
        v3=parity_W(mr+i,j,3); v4=parity_W(mr+i,j,4)
        p1=v1+v2; p2=v3+v4; p3=v1-v2; p4=v3-v4
        WD(rows(i),column(1,j))=p1+p2; WD(rows(i),column(2,j))=p3+p4
        WD(rows(i),column(3,j))=p1-p2; WD(rows(i),column(4,j))=p3-p4
      enddo
    enddo
    if (n_kernel>2_8) then   ! SLPn, DLPn rows, unpacked as WS, WD
      allocate(W_k(mr,n))
      do k=3,n_kernel
        do j=1,64
          !GCC$ ivdep
          do i=1,mr
            v1=parity_W(i+mr*(k-1),j,1); v2=parity_W(i+mr*(k-1),j,2)
            v3=parity_W(i+mr*(k-1),j,3); v4=parity_W(i+mr*(k-1),j,4)
            p1=v1+v2; p2=v3+v4; p3=v1-v2; p4=v3-v4
            W_k(i,column(1,j))=p1+p2; W_k(i,column(2,j))=p3+p4
            W_k(i,column(3,j))=p1-p2; W_k(i,column(4,j))=p3-p4
          enddo
        enddo
        if (k<=k_Dn) then
          WSn(rows,:,k-k_Sn)=W_k
        else
          WDn(rows,:,i_Dn(k-k_Dn),j_Dn(k-k_Dn))=W_k; WDn(rows,:,j_Dn(k-k_Dn),i_Dn(k-k_Dn))=W_k
        endif
      enddo
      deallocate(W_k)
    endif
    deallocate(parity_K,parity_W,Y)
    deallocate(Pu,Pv,X,Xtmp,rt,rp,sq_x,sq_nx,sq_w,target)
    deallocate(g,wg,L,rows)
  enddo

  remaining=inside128
  off_column4=.false.
  t0_local=t0; h_local=h
  if (any(inside128)) then
    do iv=1,16
      do iu=1,16
        source_t(:,iu+16*(iv-1))=[g16(iu),g16(iv)]
      enddo
    enddo
  endif

  ! Regions 4:6: panel refinement for patches up to 10:1.
  ratio = norm2(rp0)/norm2(rt0)
  ratio = max(ratio,1d0/ratio)
  ! if (any(remaining) .and. ratio <= 10d0) then
  if (any(remaining) .and. ratio <= 8d0) then ! heuristic
    block
      integer(8), parameter :: order4(2)=[64,16] ! Common orders for regions 4:6.
      integer(8), parameter :: nv=order4(2),npq=product(order4)
      integer(8), parameter :: preset=2 ! 1: five rows; 2: three rows.
      integer(8) :: p,ip,ip_right,pr,direction,n_panel,np,cols,parameter_order(2),density_column(16,16)
      integer(8) :: panel_count(3),row,nparent,nchild,nsplit,nlayout,first_panel,last_panel
      integer(8), allocatable :: ncols(:),panel_id(:,:),panel_row(:),needed(:),panel_mirror(:),panel_pair(:,:,:)
      integer(8) :: node_order(npq),reverse_nodes(npq)
      real(8) :: box_center(3),box_frame(3,3),box_halfsize4(3),box_halfsize5(3),box_halfsize6(3)
      real(8) :: scale_box4,scale_box5,scale_box6,d_box4,d_box5,d_box6
      real(8), allocatable :: row_ratio(:),split_ratio(:),v_edges(:),panels_all(:,:)
      real(8) :: width,top_height
      real(8) :: b(4),c(2),half(2),rotation(2,2,4)
      real(8) :: tu(order4(1)),tv(nv),g_h(order4(1)),wg_h(order4(1)),g_v(nv),wg_v(nv),wq(npq)
      real(8) :: C_Xp(16,16,3),C_Xip(16,16,3,2)
      real(8) :: normal_h(npq,3),weight_h(npq)
      real(8), allocatable :: Lu(:,:,:),Lv(:,:,:),W_mid(:,:)
      logical :: inside_box4(m),inside_box5(m),inside_box6(m),region(m,3),direction_groups(m,4),active(3)
      logical, allocatable :: use_panel(:)
      ! Cylindrical bounds around the self patch for regions 4:6.
      integer(8), parameter :: p_R=3,n_R=(p_R+1)*(p_R+2)/2,n_b=68
      integer(8) :: i_R(n_R),j_R(n_R),iR,jR,kR,itarg
      real(8) :: m_box(3),H_box,gram(3),t_b(n_b,2),P_bu(n_b,16),P_bv(n_b,16),x_b(3,n_b),xyz_s(3,n+n_b),s_R(2)
      real(8) :: B_s(n+n_b,n_R),B_fit(n,n_R),c_R(n,2),r_s(2),zeta_lo,zeta_hi,ell_R(2),b_box(3,2),lapack_work(4096)
      real(8), allocatable :: P_a(:,:),P_b(:,:),B_t(:,:),R_t(:,:)
      real(8) :: P_est(2,16),X_est(3),J_est(3,2),n_est(3),t_direction(2,m)
      logical :: in_slab(m),beside(m)
      external :: dgels

      box_center=xc; box_frame=frame
      t_direction=t0_local
      if (.false.) then ! old rule, turns out to be not valid for sphere, near triple junction
        ! scale_box4=0.95d0; d_box4=0.15d0*maxval(r(1:2))
        scale_box4=0.95d0; d_box4=0.0425d0*maxval(r(1:2))
        ! scale_box5=0.95d0; d_box5=0.080d0*maxval(r(1:2))
        scale_box5=0.65d0; d_box5=0.0d0*maxval(r(1:2))
        ! scale_box6=0.95d0; d_box6=0.0d0*maxval(r(1:2))
        ! scale_box6=0.90d0; d_box6=0.0d0*maxval(r(1:2))
        scale_box6=0.60d0; d_box6=0.0d0*maxval(r(1:2))
        box_halfsize4=scale_box4*r+d_box4
        box_halfsize5=scale_box5*r+d_box5
        box_halfsize6=scale_box6*r+d_box6
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        inside_box4=all(abs(xyz)<=spread(box_halfsize4,2,m),dim=1)
        inside_box5=all(abs(xyz)<=spread(box_halfsize5,2,m),dim=1)
        inside_box6=all(abs(xyz)<=spread(box_halfsize6,2,m),dim=1)
      else
        ! Boxes 4:6 as cylinders around the self patch.
        ! (xi,eta,zeta) = box_frame^T (x-box_center); (u,v) = R(xi,eta): least-squares inverse fit of degree p_R
        ! on the source nodes. box_k = {|R_i| <= r_s,i + m_k*l_other/l_i, zeta_lo-H_box <= zeta <= zeta_hi+H_box};
        ! r_s = extent of R over self, so self is inside every box.
        ! Region 7 = box6; a target with trusted (t0,h) must also satisfy |t0_i| <= b6_i
        ! (keeps targets beside a curved edge out of Region 7).
        !
        ! things to tweak inside, heuristic: p_R (parameter above), m_box, H_box
        m_box=[0.10d0,0.06d0,0.02d0]
        gram=[dot_product(rt0,rt0),dot_product(rt0,rp0),dot_product(rp0,rp0)]
        H_box=0.6d0*sqrt(0.5d0*(gram(1)+gram(3))+sqrt(0.25d0*(gram(1)-gram(3))**2+gram(2)**2)) ! 0.6*||[rt0,rp0]||_2
        !
        kR=0 ! total degree <= p_R
        do jR=1,p_R+1
          do iR=1,p_R+2-jR
            kR=kR+1; i_R(kR)=iR; j_R(kR)=jR
          enddo
        enddo
        t_b(1:16,1)=g16; t_b(1:16,2)=-1d0 ! self boundary: GL16 on each edge + corners
        t_b(17:32,1)=1d0; t_b(17:32,2)=g16
        t_b(33:48,1)=g16; t_b(33:48,2)=1d0
        t_b(49:64,1)=-1d0; t_b(49:64,2)=g16
        t_b(65:68,1)=[-1d0,1d0,1d0,-1d0]; t_b(65:68,2)=[-1d0,-1d0,1d0,1d0]
        call legendrevec(15_8,(t_b(:,1)-geometry_center(1))/geometry_scale(1),P_bu)
        call legendrevec(15_8,(t_b(:,2)-geometry_center(2))/geometry_scale(2),P_bv)
        do k=1,3
          x_b(k,:)=sum(matmul(P_bu,C_X(:,:,k))*P_bv,dim=2)
        enddo
        xyz_s(:,1:n)=matmul(transpose(box_frame),sx-spread(box_center,2,n)) ! self samples: n nodes, then boundary
        xyz_s(:,n+1:)=matmul(transpose(box_frame),x_b-spread(box_center,2,n_b))
        s_R=maxval(abs(xyz_s(1:2,1:n)),dim=2)
        allocate(P_a(n+n_b,p_R+1),P_b(n+n_b,p_R+1))
        call legendrevec(p_R,xyz_s(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz_s(2,:)/s_R(2),P_b)
        B_s=P_a(:,i_R)*P_b(:,j_R)
        deallocate(P_a,P_b)
        B_fit=B_s(1:n,:); c_R=transpose(source_t) ! fit on the n source nodes
        call dgels('N',n,n_R,2_8,B_fit,n,c_R,n,lapack_work,4096_8,info)
        r_s=maxval(abs(matmul(B_s,c_R(1:n_R,:))),dim=1) ! self extent of R (nodes and boundary)
        zeta_lo=minval(xyz_s(3,:)); zeta_hi=maxval(xyz_s(3,:))
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        allocate(P_a(m,p_R+1),P_b(m,p_R+1),B_t(m,n_R),R_t(m,2))
        call legendrevec(p_R,xyz(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz(2,:)/s_R(2),P_b)
        B_t=P_a(:,i_R)*P_b(:,j_R)
        R_t=matmul(B_t,c_R(1:n_R,:))
        do itarg=1,m
          if (.not.root_ok(itarg) .or. .not.all(ieee_is_finite(t0_local(:,itarg)))) t_direction(:,itarg)=R_t(itarg,:)
        enddo
        in_slab=xyz(3,:)>=zeta_lo-H_box .and. xyz(3,:)<=zeta_hi+H_box
        ! facing-edge guard: the margin across an edge scales with that edge's length (anisotropic patches)
        ell_R=[norm2(rt0),norm2(rp0)] ! physical length per unit t1, t2
        b_box(:,1)=r_s(1)+m_box*ell_R(2)/ell_R(1) ! rows box4/5/6, columns t1/t2
        b_box(:,2)=r_s(2)+m_box*ell_R(1)/ell_R(2)
        ! Off-surface targets in column4 bypass the slab and regions 4:6.
        off_column4=remaining .and. abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. h_local/=0d0
        do itarg=1,m
          if (.not.off_column4(itarg)) cycle
          if (all(ieee_is_finite(t0_local(:,itarg))) .and. ieee_is_finite(h_local(itarg))) cycle
          t0_local(:,itarg)=R_t(itarg,:)
          call legendrevec(15_8,(t0_local(:,itarg)-geometry_center)/geometry_scale,P_est)
          do k=1,3
            X_est(k)=dot_product(matmul(P_est(1,:),C_X(:,:,k)),P_est(2,:))
            J_est(k,1)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,1)),P_est(2,:))
            J_est(k,2)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,2)),P_est(2,:))
          enddo
          n_est=[J_est(2,2)*J_est(3,1)-J_est(3,2)*J_est(2,1), &
                 J_est(3,2)*J_est(1,1)-J_est(1,2)*J_est(3,1), &
                 J_est(1,2)*J_est(2,1)-J_est(2,2)*J_est(1,1)]
          n_est=n_est/norm2(n_est)
          h_local(itarg)=dot_product(tx(:,itarg)-X_est,n_est)
        enddo
        off_column4=off_column4 .and. ieee_is_finite(h_local) .and. all(ieee_is_finite(t0_local),dim=1)
        remaining=remaining .and. .not.off_column4
        inside_box4=abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. in_slab
        inside_box5=abs(R_t(:,1))<=b_box(2,1) .and. abs(R_t(:,2))<=b_box(2,2) .and. in_slab
        inside_box6=abs(R_t(:,1))<=b_box(3,1) .and. abs(R_t(:,2))<=b_box(3,2) .and. in_slab
        ! Region 7 only above the patch: a target with trusted t0,h must also have its root inside box6
        beside=root_ok .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local) .and. &
               (abs(t0_local(1,:))>b_box(3,1) .or. abs(t0_local(2,:))>b_box(3,2))
        inside_box6=inside_box6 .and. .not.beside
        deallocate(P_a,P_b,B_t,R_t)
      endif
      region(:,1)=remaining .and. .not.inside_box4; remaining=remaining .and. inside_box4
      region(:,2)=remaining .and. .not.inside_box5; remaining=remaining .and. inside_box5
      region(:,3)=remaining .and. .not.inside_box6; remaining=remaining .and. inside_box6
      where (region(:,1)) classification=4_8
      where (region(:,2)) classification=5_8
      where (region(:,3)) classification=6_8

      ! List row heights and even column counts from the target-facing edge down.
      select case (preset)
      case (1)
        row_ratio=[1d0,1d0,1d0,3d0,3d0]
        ncols=[4,4,4,2,2]
        split_ratio=[1d0,1d0,1d0]
      case (2)
        row_ratio=[1d0,2d0,6d0]
        ncols=[4,4,2]
        split_ratio=[1d0,2d0]
      end select
      nsplit=size(split_ratio)
      panel_count(1)=sum(ncols); nlayout=panel_count(1); nparent=ncols(1)
      do layer=2,3
        nchild=2*nparent*nsplit
        panel_count(layer)=panel_count(layer-1)-nparent+nchild
        nlayout=nlayout+nchild; nparent=2*nparent
      enddo
      allocate(panels_all(4,nlayout),panel_id(maxval(panel_count),3),panel_mirror(nlayout))
      allocate(panel_pair(2,maxval(panel_count)/2,3),panel_row(nlayout),needed(nlayout),use_panel(nlayout))
      v_edges=[(1d0-2d0*sum(row_ratio(:row))/sum(row_ratio),row=0,size(row_ratio))]
      last_panel=0
      do layer=1,3
        if (layer>1) then
          nparent=ncols(1); top_height=1d0-v_edges(2)
          v_edges=[(1d0-top_height*sum(split_ratio(:row))/sum(split_ratio),row=0,nsplit)]
          ncols=spread(2*nparent,1,nsplit)
        endif
        first_panel=last_panel+1
        do row=1,size(ncols)
          width=2d0/ncols(row)
          do j=1,ncols(row)
            p=last_panel+j
            panels_all(:,p)=[-1d0+(j-1)*width,-1d0+j*width,v_edges(row+1),v_edges(row)]
            panel_mirror(p)=last_panel+ncols(row)+1-j
          enddo
          last_panel=last_panel+ncols(row)
        enddo
        nchild=last_panel-first_panel+1
        panel_id(1:nchild,layer)=[(p,p=first_panel,last_panel)]
        ! Lower rows retain their IDs from the previous layout.
        if (layer>1) panel_id(nchild+1:panel_count(layer),layer)=panel_id(nparent+1:panel_count(layer-1),layer-1)
      enddo
      ! Pair canonical left/right rectangles; geometry and targets need no symmetry.
      do layer=1,3
        is=0
        do p=1,panel_count(layer)
          ip=panel_id(p,layer)
          if (panels_all(1,ip)+panels_all(2,ip)>=0d0) cycle
          is=is+1
          panel_pair(:,is,layer)=[p,findloc(panel_id(1:panel_count(layer),layer),panel_mirror(ip),dim=1)]
        enddo
      enddo
      direction_groups(:,1)=t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,2)=t_direction(1,:)>abs(t_direction(2,:))
      direction_groups(:,3)=t_direction(2,:)<0 .and. -t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,4)=-t_direction(1,:)>abs(t_direction(2,:))
      rotation(:,:,1)=reshape([1d0,0d0,0d0,1d0],[2,2])
      rotation(:,:,2)=reshape([0d0,-1d0,1d0,0d0],[2,2])
      rotation(:,:,3)=-rotation(:,:,1); rotation(:,:,4)=-rotation(:,:,2)
      nq=order4(1)
      call gauss(nq,g_h,wg_h)
      call gauss(nv,g_v,wg_v)
      wq=reshape(spread(wg_h,2,nv)*spread(wg_v,1,nq),[npq])

      reverse_nodes=[(nq-mod(j-1,nq)+nq*((j-1)/nq),j=1,npq)]

      do direction=1,4
        use_panel=.false.
        do layer=1,3
          target_groups(:,layer)=region(:,layer) .and. direction_groups(:,direction)
          active(layer)=any(target_groups(:,layer))
          if (active(layer)) use_panel(panel_id(1:panel_count(layer),layer))=.true.
        enddo
        n_panel=count(use_panel)
        if (n_panel==0) cycle
        needed(1:n_panel)=pack([(p,p=1,nlayout)],use_panel)
        panel_row=0; panel_row(needed(1:n_panel))=[(p,p=1,n_panel)]
        parameter_order=maxloc(abs(rotation(:,:,direction)),dim=1)
        C_Xp=C_X; C_Xip=C_Xi
        density_column=reshape([(j,j=1,n)],[16,16])
        if (parameter_order(1)==2) then
          do k=1,3
            C_Xp(:,:,k)=transpose(C_X(:,:,k))
            do j=1,2
              C_Xip(:,:,k,j)=transpose(C_Xi(:,:,k,j))
            enddo
          enddo
          density_column=transpose(density_column)
        endif
        allocate(sq_x(npq*n_panel,3),sq_nx(npq*n_panel,3),sq_w(npq*n_panel))
        allocate(Lu(nq,16,n_panel),Lv(nv,16,n_panel),Pu(nq,16),Pv(nv,16),X(nq,nv),Xtmp(nq,16))
        allocate(rt(nq,nv,3),rp(nq,nv,3))
        do p=1,n_panel
          b=panels_all(:,needed(p)); c=(b([1,3])+b([2,4]))/2; half=(b([2,4])-b([1,3]))/2
          cols=(p-1)*npq
          node_order=[(j,j=1,npq)]
          if (b(1)+b(2)>0d0) node_order=reverse_nodes
          tu=rotation(parameter_order(1),1,direction)*(c(1)+half(1)*g_h)
          tv=rotation(parameter_order(2),2,direction)*(c(2)+half(2)*g_v)
          call legendrevec(15_8,(tu-geometry_center(parameter_order(1))) &
            /geometry_scale(parameter_order(1)),Pu)
          call legendrevec(15_8,(tv-geometry_center(parameter_order(2))) &
            /geometry_scale(parameter_order(2)),Pv)
          do k=1,3
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xp(:,:,k),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,X,nq)
            sq_x(cols+node_order,k)=reshape(X,[npq])
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,1),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rt(:,:,k),nq)
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,2),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rp(:,:,k),nq)
          enddo
          normal_h(:,1)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[npq])
          normal_h(:,2)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[npq])
          normal_h(:,3)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[npq])
          weight_h=sqrt(sum(normal_h**2,dim=2))
          do k=1,3
            sq_nx(cols+node_order,k)=normal_h(:,k)/weight_h
          enddo
          sq_w(cols+node_order)=product(half)*weight_h*wq
          ! Only the canonical-left panel supplies the paired density operators.
          if (b(1)+b(2)<0d0) then
            call legendrevec(15_8,tu,Pu)
            call legendrevec(15_8,tv,Pv)
            call dgemm('N','T',nq,16_8,16_8,1d0,Pu,nq,density_map,16_8,0d0,Lu(:,:,p),nq)
            Xtmp=Lu(:,:,p)
            Lu(:,1:8,p)=(Xtmp(:,1:8)+Xtmp(:,16:9:-1))/2d0
            Lu(:,9:16,p)=(Xtmp(:,1:8)-Xtmp(:,16:9:-1))/2d0
            call dgemm('N','T',nv,16_8,16_8,1d0,Pv,nv,density_map,16_8,0d0,Lv(:,:,p),nv)
          endif
        enddo
        deallocate(Pu,Pv,X,Xtmp,rt,rp)

        do layer=1,3
          if (.not.active(layer)) cycle
          mr=count(target_groups(:,layer)); rows=pack([(i,i=1,m)],target_groups(:,layer))
          target=transpose(tx(:,rows)); np=panel_count(layer)
          ! Column i+mr*(k-1)+n_kernel*mr*(iv-1): target i, kernel k, short node iv.
          allocate(W_mid(n_kernel*mr,n),Y(8,n_kernel*mr*nv))
          allocate(parity_K(nq,n_kernel*mr*nv,2),parity_W(8,n_kernel*mr*16,2))
          parity_W=0d0
          do p=1,np/2
            ip=panel_row(panel_id(panel_pair(1,p,layer),layer))
            ip_right=panel_row(panel_id(panel_pair(2,p,layer),layer))
            do iv=1,nv
              do i=1,mr
                cols=i+n_kernel*mr*(iv-1)
                !GCC$ ivdep
                !GCC$ unroll 1
                do iu=1,nq
                  j=(ip-1)*npq+iu+nq*(iv-1)
                  pr=(ip_right-1)*npq+iu+nq*(iv-1)
                  dx=sq_x(j,1)-target(i,1); dy=sq_x(j,2)-target(i,2); dz=sq_x(j,3)-target(i,3)
                  distance2_inv=1d0/(dx*dx+dy*dy+dz*dz)
                  v1=(sq_w(j)*inv2pi)*sqrt(distance2_inv)
                  d1=v1*(dx*sq_nx(j,1)+dy*sq_nx(j,2)+dz*sq_nx(j,3))*distance2_inv
                  dx=sq_x(pr,1)-target(i,1); dy=sq_x(pr,2)-target(i,2); dz=sq_x(pr,3)-target(i,3)
                  distance2_inv=1d0/(dx*dx+dy*dy+dz*dz)
                  v2=(sq_w(pr)*inv2pi)*sqrt(distance2_inv)
                  d2=v2*(dx*sq_nx(pr,1)+dy*sq_nx(pr,2)+dz*sq_nx(pr,3))*distance2_inv
                  parity_K(iu,cols,1)=v1+v2; parity_K(iu,cols,2)=v1-v2
                  parity_K(iu,mr+cols,1)=d1+d2; parity_K(iu,mr+cols,2)=d1-d2
                enddo
              enddo
            enddo
            if (n_kernel>2_8) then   ! SLPn, DLPn columns: the same mirrored-pair sums
              allocate(node_d(nq,3,2),node_r2(nq,2),node_f(nq,2))
              do iv=1,nv
                do i=1,mr
                  cols=i+n_kernel*mr*(iv-1)
                  !GCC$ ivdep
                  do iu=1,nq
                    j=(ip-1)*npq+iu+nq*(iv-1)
                    pr=(ip_right-1)*npq+iu+nq*(iv-1)
                    dx=sq_x(j,1)-target(i,1); dy=sq_x(j,2)-target(i,2); dz=sq_x(j,3)-target(i,3)
                    distance2_inv=1d0/(dx*dx+dy*dy+dz*dz)
                    node_d(iu,1,1)=dx; node_d(iu,2,1)=dy; node_d(iu,3,1)=dz; node_r2(iu,1)=distance2_inv
                    node_f(iu,1)=((sq_w(j)*inv2pi)*sqrt(distance2_inv))*distance2_inv
                    dx=sq_x(pr,1)-target(i,1); dy=sq_x(pr,2)-target(i,2); dz=sq_x(pr,3)-target(i,3)
                    distance2_inv=1d0/(dx*dx+dy*dy+dz*dz)
                    node_d(iu,1,2)=dx; node_d(iu,2,2)=dy; node_d(iu,3,2)=dz; node_r2(iu,2)=distance2_inv
                    node_f(iu,2)=((sq_w(pr)*inv2pi)*sqrt(distance2_inv))*distance2_inv
                  enddo
                  do kn=1,3*ifSn
                    k=k_Sn+kn
                    !GCC$ ivdep
                    do iu=1,nq
                      v1=node_f(iu,1)*node_d(iu,kn,1); v2=node_f(iu,2)*node_d(iu,kn,2)
                      parity_K(iu,cols+mr*(k-1),1)=v1+v2; parity_K(iu,cols+mr*(k-1),2)=v1-v2
                    enddo
                  enddo
                  do kn=1,6*ifDn
                    k=k_Dn+kn; a_c=i_Dn(kn); b_c=j_Dn(kn); delta_c=merge(1d0,0d0,a_c==b_c)
                    !GCC$ ivdep
                    do iu=1,nq
                      v1=node_f(iu,1)*(3d0*node_d(iu,a_c,1)*node_d(iu,b_c,1)*node_r2(iu,1)-delta_c)
                      v2=node_f(iu,2)*(3d0*node_d(iu,a_c,2)*node_d(iu,b_c,2)*node_r2(iu,2)-delta_c)
                      parity_K(iu,cols+mr*(k-1),1)=v1+v2; parity_K(iu,cols+mr*(k-1),2)=v1-v2
                    enddo
                  enddo
                enddo
              enddo
              deallocate(node_d,node_r2,node_f)
            endif
            ! Reduce the paired coordinate to eight density nodes first.
            do mode=1,2
              call dgemm('T','N',8_8,n_kernel*mr*nv,nq,1d0,Lu(1,1+8*(mode-1),ip),nq, &
                parity_K(1,1,mode),nq,0d0,Y,8_8)
              call dgemm('N','N',8*n_kernel*mr,16_8,nv,1d0,Y,8*n_kernel*mr,Lv(:,:,ip),nv, &
                1d0,parity_W(1,1,mode),8*n_kernel*mr)
            enddo
          enddo
          do iv=1,16
            do iu=1,8
              W_mid(:,density_column(iu,iv))=parity_W(iu,1+n_kernel*mr*(iv-1):n_kernel*mr*iv,1) &
                +parity_W(iu,1+n_kernel*mr*(iv-1):n_kernel*mr*iv,2)
              W_mid(:,density_column(17-iu,iv))=parity_W(iu,1+n_kernel*mr*(iv-1):n_kernel*mr*iv,1) &
                -parity_W(iu,1+n_kernel*mr*(iv-1):n_kernel*mr*iv,2)
            enddo
          enddo
          WS(rows,:)=W_mid(1:mr,:); WD(rows,:)=W_mid(mr+1:2*mr,:)
          do kn=1,3*ifSn
            k=k_Sn+kn
            WSn(rows,:,kn)=W_mid(mr*(k-1)+1:mr*k,:)
          enddo
          do kn=1,6*ifDn
            k=k_Dn+kn
            WDn(rows,:,i_Dn(kn),j_Dn(kn))=W_mid(mr*(k-1)+1:mr*k,:)
            WDn(rows,:,j_Dn(kn),i_Dn(kn))=W_mid(mr*(k-1)+1:mr*k,:)
          enddo
          deallocate(W_mid,Y,parity_K,parity_W,target,rows)
        enddo
        deallocate(sq_x,sq_nx,sq_w,Lu,Lv)
      enddo
    end block
  endif

  ! Region 7: tensor sinh for targets with finite source-chart coordinates.
  remaining=remaining .or. off_column4
  where (remaining) classification=7_8
  use_sinh=remaining .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local)
  if (any(use_sinh)) then
    block
      integer(8), parameter :: nlong=192,nshort=32
      integer(8) :: axis(2),n_axis(2),n_panels,it,itarg,np,first,last,l,j_sigma
      real(8) :: c(2),scale(2),J_center(3,2),a2,b2,shear,area,t(2),dc(2),center(2),beta(2)
      real(8) :: C_pos(16,16,3),Xi(16,16,3,2),Xij(16,16,3,3),orientation
      real(8) :: c_s(2),scale_s(2),Xi_s(16,16,3,2),Xij_s(16,16,3,3)
      real(8) :: c_b(2),scale_b(2),Xi_b(16,16,3,2),Xij_b(16,16,3,3)
      logical :: use_typeb,inside
      real(8) :: g32(32),wg32(32),g96(96),wg96(96),g192(192),wg192(192),sigma(16),w_sigma(16)
      real(8) :: edges(3),g_axis(nlong),wg_axis(nlong),a,b,eta0,span,step,eta,z,A0,r0,p0_axis,r0_axis
      real(8) :: alpha,ell,sgn,f,q(nlong,2),d(nlong,2),w(nlong,2),P(nlong,16,2)
      real(8) :: line(nlong*16),B_line(nlong*16,16),T1(nlong,16,2),T0L(nlong,16)
      real(8) :: P0(2,16),J0(3,2),n0(3),offset(3)
      real(8) :: dL(nlong,nshort),dS(nlong,nshort),HSS(nshort),HLS(nlong,nshort),HLL(nlong,nshort)
      real(8) :: XL(nlong,nshort,3),XS(nlong,nshort,3),chord(nlong,nshort,3),R(nlong,nshort,3)
      real(8) :: Nvec(nlong,nshort,3),distance2(nlong,nshort),rinv(nlong,nshort),omega(nlong,nshort)
      real(8) :: kernel_S(nlong,nshort),kernel_D(nlong,nshort),LL(nlong,16),LS(nshort,16)
      real(8) :: Y_leg(16,nshort),W_S(16,16),W_D(16,16)
      ! SLPn, DLPn: kernel_n(:,:,kn), kn = SLPn x, y, z then DLPn (a,b); 1/r^2 and kernel_S/r^2
      real(8), allocatable :: kernel_n(:,:,:),r2_inv(:,:),factor_r2(:,:)
      logical :: naive
      ! Restructured evaluation (same nodes and formulas as the reference branches below).
      integer(8) :: il,js
      real(8) :: z_line(nlong),P_prev(nlong),P_cur(nlong),P_next(nlong),T0_line(nlong,16),w_T0,w_T1,rec1,rec2
      real(8) :: C_XLS(16,6,16),C_HLS(16,3,16),C_HLL(16,3,16)
      real(8) :: S_XLS(16,6,nshort),S_HLS(16,3,nshort),S_HLL(16,3,nshort)
      real(8) :: XLS(nlong,6,nshort),HLS3(nlong,3,nshort),HLL3(nlong,3,nshort),HSS3(nshort,3)
      real(8) :: d_L,d_S,hs1,hs2,hs3,r1,r2,r3,c1,c2,c3,Nv1,Nv2,Nv3,dist2,r_inv,om
      interface
        pure function expm1(x) bind(C,name='expm1') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function expm1
        pure function log1p(x) bind(C,name='log1p') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function log1p
      end interface

      axis=[1_8,2_8]
      if (dot_product(rp0,rp0)>dot_product(rt0,rt0)) axis=[2_8,1_8]
      orientation=real(3-2*axis(1),8)
      c=geometry_center(axis); scale=geometry_scale(axis)
      J_center(:,1)=rt0; J_center(:,2)=rp0; J_center=J_center(:,axis)
      a2=dot_product(J_center(:,1),J_center(:,1)); b2=dot_product(J_center(:,2),J_center(:,2))
      shear=dot_product(J_center(:,1),J_center(:,2))
      area=norm2([rt0(2)*rp0(3)-rt0(3)*rp0(2),rt0(3)*rp0(1)-rt0(1)*rp0(3),rt0(1)*rp0(2)-rt0(2)*rp0(1)])
      C_pos=C_X
      Xi=C_Xi(:,:,:,axis); Xij=C_Xij(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
      if (axis(1)==2) then
        do k=1,3
          C_pos(:,:,k)=transpose(C_pos(:,:,k))
          do j=1,2
            Xi(:,:,k,j)=transpose(Xi(:,:,k,j))
          enddo
          do j=1,3
            Xij(:,:,k,j)=transpose(Xij(:,:,k,j))
          enddo
        enddo
      endif
      c_s=c; scale_s=scale; Xi_s=Xi; Xij_s=Xij ! source chart
      use_typeb=all(geometry_scale_typeb>0d0)
      if (use_typeb) then ! type-b chart: the source sits inside it; used when t0 falls inside
        c_b=geometry_center_typeb(axis); scale_b=geometry_scale_typeb(axis)
        Xi_b=C_Xi_typeb(:,:,:,axis); Xij_b=C_Xij_typeb(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
        if (axis(1)==2) then
          do k=1,3
            do j=1,2
              Xi_b(:,:,k,j)=transpose(Xi_b(:,:,k,j))
            enddo
            do j=1,3
              Xij_b(:,:,k,j)=transpose(Xij_b(:,:,k,j))
            enddo
          enddo
        endif
      endif
      call gauss(32_8,g32,wg32)
      call gauss(96_8,g96,wg96)
      call gauss(192_8,g192,wg192)
      sigma=(g16+1d0)/2; w_sigma=wg16/2
      rows=pack([(j,j=1,m)],use_sinh)
      n_axis=[nlong,nshort]
      if (n_kernel>2_8) allocate(kernel_n(nlong,nshort,n_kernel-2),r2_inv(nlong,nshort),factor_r2(nlong,nshort))

      do it=1,size(rows,kind=8)
        itarg=rows(it)
        t=t0_local(axis,itarg)
        ! naive=abs(t(2))>2.5d0
        ! naive=off_column4(itarg) .or. abs(t(2))>1.5d0 ! heuristic
        inside=.false.
        if (use_typeb) inside=all(abs(t-c_b)<=scale_b*(1d0+1d-12)) ! t0 inside the type-b chart: geometry is valid at t0
        if (inside .and. .not.off_column4(itarg) .and. root_ok(itarg)) then
          c=c_b; scale=scale_b; Xi=Xi_b; Xij=Xij_b; naive=.false.
        else
          c=c_s; scale=scale_s; Xi=Xi_s; Xij=Xij_s; naive=.not.root_ok(itarg) .or. off_column4(itarg) .or. abs(t(2))>1.5d0
        endif
        dc=min(1d0,max(-1d0,t))-t
        center=t-[shear/a2*dc(2),shear/b2*dc(1)]
        beta=[hypot(h_local(itarg),area/sqrt(a2)*dc(2))/sqrt(a2),hypot(h_local(itarg),area/sqrt(b2)*dc(1))/sqrt(b2)]
        if (.not.root_ok(itarg) .and. (any(.not.ieee_is_finite(center)) .or. any(.not.ieee_is_finite(beta)) .or. &
            any(beta==0d0 .and. abs(center)<=1d0))) then
          use_sinh(itarg)=.false.
          cycle
        endif

        do k=1,2
          n_panels=1; edges(1:2)=[-1d0,1d0]
          if (k==2) then
            np=nshort; g_axis(:np)=g32; wg_axis(:np)=wg32
          else if (abs(center(k))<1d0) then
            np=nlong/2; g_axis(:np)=g96; wg_axis(:np)=wg96
            n_panels=2; edges=[-1d0,center(k),1d0]
          else
            np=nlong; g_axis=g192; wg_axis=wg192
          endif
          do j=1,n_panels
            a=edges(j); b=edges(j+1)
            first=(j-1)*np+1; last=j*np
            if (beta(k)>0d0) then
              eta0=asinh((a-center(k))/beta(k))
              span=asinh((b-center(k))/beta(k))-eta0
              if (abs(center(k))>=2d0) then
                span=2*atanh((b-a)/(hypot(a-center(k),beta(k))+hypot(b-center(k),beta(k))))
              endif
              A0=a-center(k); r0=hypot(A0,beta(k))
              if (A0>=0d0) then
                p0_axis=r0+A0; r0_axis=beta(k)*(beta(k)/p0_axis)
              else
                r0_axis=r0-A0; p0_axis=beta(k)*(beta(k)/r0_axis)
              endif
              do i=first,last
                step=(g_axis(i-first+1)+1d0)*(span/2)
                eta=eta0+step; z=beta(k)*sinh(eta)
                w(i,k)=wg_axis(i-first+1)*((span/2)*beta(k)*cosh(eta))
                if (abs(center(k))<=1d0) then
                  q(i,k)=center(k)+z
                else
                  q(i,k)=a+(p0_axis*expm1(step)-r0_axis*expm1(-step))/2
                endif
                d(i,k)=(center(k)-t(k))+z
              enddo
            else
              alpha=abs(center(k)); sgn=sign(1d0,center(k))
              ell=-log1p(2/(alpha-1d0))
              do i=first,last
                f=(sgn*g_axis(i-first+1)+1d0)/2
                z=-sgn*(alpha+1d0)*exp(f*ell)
                q(i,k)=sgn*(-1d0-(alpha+1d0)*expm1(f*ell))
                w(i,k)=wg_axis(i-first+1)*(-abs(z)*ell/2)
                d(i,k)=(center(k)-t(k))+z
              enddo
            endif
          enddo
          np=n_axis(k)
          call legendrevec(15_8,(q(:np,k)-c(k))/scale(k),P(:np,:,k))
          if (.not.naive) then
            ! Change .false. to .true. (here and below) to force the original evaluation for comparison.
            if (.false.) then
            do j_sigma=1,16
              line((j_sigma-1)*np+1:j_sigma*np)=(t(k)+d(:np,k)*sigma(j_sigma)-c(k))/scale(k)
            enddo
            call legendrevec(15_8,line(:np*16),B_line(:np*16,:))
            do l=1,16
              T1(:np,l,k)=matmul(reshape(B_line(:np*16,l),[np,16_8]),w_sigma*sigma)
              if (k==1) T0L(:,l)=matmul(reshape(B_line(:,l),[nlong,16_8]),w_sigma)
            enddo
            else
              ! Same GL16 sums over sigma, accumulated inside the Legendre recurrence:
              ! T1(i,l) = sum_j w_sigma(j)*sigma(j)*P_{l-1}(z_ij), T0(i,l) = sum_j w_sigma(j)*P_{l-1}(z_ij).
              T1(:np,:,k)=0d0; T0_line(:np,:)=0d0
              do j_sigma=1,16
                w_T1=w_sigma(j_sigma)*sigma(j_sigma); w_T0=w_sigma(j_sigma)
                do i=1,np
                  z_line(i)=(t(k)+d(i,k)*sigma(j_sigma)-c(k))/scale(k)
                  P_prev(i)=1d0; P_cur(i)=z_line(i)
                  T1(i,1,k)=T1(i,1,k)+w_T1; T1(i,2,k)=T1(i,2,k)+P_cur(i)*w_T1
                  T0_line(i,1)=T0_line(i,1)+w_T0; T0_line(i,2)=T0_line(i,2)+P_cur(i)*w_T0
                enddo
                do l=1,14
                  rec1=real(2*l+1,8)/real(l+1,8); rec2=real(l,8)/real(l+1,8)
                  !GCC$ ivdep
                  do i=1,np
                    P_next(i)=rec1*z_line(i)*P_cur(i)-rec2*P_prev(i)
                    T1(i,l+2,k)=T1(i,l+2,k)+P_next(i)*w_T1
                    T0_line(i,l+2)=T0_line(i,l+2)+P_next(i)*w_T0
                    P_prev(i)=P_cur(i); P_cur(i)=P_next(i)
                  enddo
                enddo
              enddo
              if (k==1) T0L=T0_line
            endif
          endif
        enddo

        ! Naive targets keep the original evaluation; .false. -> .true. forces it for every target.
        if (.false. .or. naive) then
        do k=1,3
          XL(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,1)),transpose(P(:nshort,:,2)))
          XS(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,2)),transpose(P(:nshort,:,2)))
        enddo
        if (naive) then
          do k=1,3
            chord(:,:,k)=matmul(matmul(P(:,:,1),C_pos(:,:,k)),transpose(P(:nshort,:,2)))-tx(k,itarg)
          enddo
          R=chord
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          dL=spread(d(:,1),2,nshort); dS=spread(d(:nshort,2),1,nlong)
          do k=1,3
            HSS=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
            HLS=matmul(matmul(T0L,Xij(:,:,k,2)),transpose(P(:nshort,:,2)))
            HLL=matmul(matmul(T1(:,:,1),Xij(:,:,k,1)),transpose(P(:nshort,:,2)))
            R(:,:,k)=offset(k)-dS**2*spread(HSS,1,nlong)-dL*dS*HLS-dL**2*HLL
            chord(:,:,k)=dL*XL(:,:,k)+dS*XS(:,:,k)+R(:,:,k)
          enddo
        endif

        Nvec(:,:,1)=orientation*(XS(:,:,2)*XL(:,:,3)-XS(:,:,3)*XL(:,:,2))
        Nvec(:,:,2)=orientation*(XS(:,:,3)*XL(:,:,1)-XS(:,:,1)*XL(:,:,3))
        Nvec(:,:,3)=orientation*(XS(:,:,1)*XL(:,:,2)-XS(:,:,2)*XL(:,:,1))
        distance2=sum(chord**2,dim=3)
        rinv=1d0/sqrt(distance2)
        omega=inv2pi*spread(w(:,1),2,nshort)*spread(w(:nshort,2),1,nlong)
        kernel_S=omega*sqrt(sum(Nvec**2,dim=3))*rinv
        kernel_D=omega*sum(Nvec*R,dim=3)*rinv/distance2
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          ! Coefficients stacked as (long mode, field, short mode): XL (fields 1:3) and XS (4:6); HLS; HLL.
          do l=1,16
            do k=1,3
              C_XLS(:,k,l)=Xi(:,l,k,1); C_XLS(:,3+k,l)=Xi(:,l,k,2)
              C_HLS(:,k,l)=Xij(:,l,k,2); C_HLL(:,k,l)=Xij(:,l,k,1)
            enddo
          enddo
          ! X = P_long*(C*P_short^T): the short direction first, then one GEMM per group along the long direction.
          call dgemm('N','T',96_8,nshort,16_8,1d0,C_XLS,96_8,P(1,1,2),nlong,0d0,S_XLS,96_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLS,48_8,P(1,1,2),nlong,0d0,S_HLS,48_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLL,48_8,P(1,1,2),nlong,0d0,S_HLL,48_8)
          call dgemm('N','N',nlong,6*nshort,16_8,1d0,P(1,1,1),nlong,S_XLS,16_8,0d0,XLS,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T0L,nlong,S_HLS,16_8,0d0,HLS3,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T1(1,1,1),nlong,S_HLL,16_8,0d0,HLL3,nlong)
          do k=1,3
            HSS3(:,k)=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
          enddo
          ! R, chord, normal and both kernels at each node, same expressions as the reference.
          do js=1,nshort
            d_S=d(js,2)
            hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
            !GCC$ ivdep
            do il=1,nlong
              d_L=d(il,1)
              r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
              r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
              r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
              c1=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
              c2=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
              c3=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
              Nv1=orientation*(XLS(il,5,js)*XLS(il,3,js)-XLS(il,6,js)*XLS(il,2,js))
              Nv2=orientation*(XLS(il,6,js)*XLS(il,1,js)-XLS(il,4,js)*XLS(il,3,js))
              Nv3=orientation*(XLS(il,4,js)*XLS(il,2,js)-XLS(il,5,js)*XLS(il,1,js))
              dist2=c1**2+c2**2+c3**2
              r_inv=1d0/sqrt(dist2)
              om=inv2pi*w(il,1)*w(js,2)
              kernel_S(il,js)=om*sqrt(Nv1**2+Nv2**2+Nv3**2)*r_inv
              kernel_D(il,js)=om*(Nv1*r1+Nv2*r2+Nv3*r3)*r_inv/dist2
            enddo
          enddo
          if (n_kernel>2_8) then   ! chord and distance2 for SLPn, DLPn, as in the loop above
            do js=1,nshort
              d_S=d(js,2)
              hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
              do il=1,nlong
                d_L=d(il,1)
                r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
                r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
                r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
                chord(il,js,1)=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
                chord(il,js,2)=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
                chord(il,js,3)=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
                distance2(il,js)=chord(il,js,1)**2+chord(il,js,2)**2+chord(il,js,3)**2
              enddo
            enddo
          endif
        endif

        if (n_kernel>2_8) then   ! SLPn, DLPn from the chord
          r2_inv=1d0/distance2
          factor_r2=kernel_S*r2_inv
          do kn=1,3*ifSn
            kernel_n(:,:,kn)=factor_r2*chord(:,:,kn)
          enddo
          do kn=1,6*ifDn
            kernel_n(:,:,3*ifSn+kn)=factor_r2*(3d0*chord(:,:,i_Dn(kn))*chord(:,:,j_Dn(kn))*r2_inv &
              -merge(1d0,0d0,i_Dn(kn)==j_Dn(kn)))
          enddo
        endif

        call legendrevec(15_8,q(:,1),P(:,:,1))
        call legendrevec(15_8,q(:nshort,2),P(:nshort,:,2))
        LL=matmul(P(:,:,1),transpose(density_map)); LS=matmul(P(:nshort,:,2),transpose(density_map))
        call dgemm('T','N',16_8,nshort,nlong,1d0,LL,nlong,kernel_S,nlong,0d0,Y_leg,16_8)
        call dgemm('N','N',16_8,16_8,nshort,1d0,Y_leg,16_8,LS,nshort,0d0,W_S,16_8)
        call dgemm('T','N',16_8,nshort,nlong,1d0,LL,nlong,kernel_D,nlong,0d0,Y_leg,16_8)
        call dgemm('N','N',16_8,16_8,nshort,1d0,Y_leg,16_8,LS,nshort,0d0,W_D,16_8)
        if (axis(1)==2) then
          W_S=transpose(W_S); W_D=transpose(W_D)
        endif
        WS(itarg,:)=reshape(W_S,[256]); WD(itarg,:)=reshape(W_D,[256])
        do kn=1,n_kernel-2   ! SLPn, DLPn, contracted as SLP, DLP
          call dgemm('T','N',16_8,nshort,nlong,1d0,LL,nlong,kernel_n(1,1,kn),nlong,0d0,Y_leg,16_8)
          call dgemm('N','N',16_8,16_8,nshort,1d0,Y_leg,16_8,LS,nshort,0d0,W_S,16_8)
          if (axis(1)==2) W_S=transpose(W_S)
          if (kn<=3*ifSn) then
            WSn(itarg,:,kn)=reshape(W_S,[256])
          else
            WDn(itarg,:,i_Dn(kn-3*ifSn),j_Dn(kn-3*ifSn))=reshape(W_S,[256])
            WDn(itarg,:,j_Dn(kn-3*ifSn),i_Dn(kn-3*ifSn))=reshape(W_S,[256])
          endif
        enddo
      enddo
      deallocate(rows)
      if (n_kernel>2_8) deallocate(kernel_n,r2_inv,factor_r2)
    end block
  endif
  remaining=remaining .and. .not.use_sinh

  ! Previous Region 7: numerical area-to-line quadrature and radial H0 geometry.
  if (.false. .and. any(remaining)) then
    block
      integer(8), parameter :: n_base=32
      integer(8) :: ninner,it,itarg,iteration,correction,ir,iedge,n_radial,n_edge,nquad
      real(8) :: B2(6,256),Bfit(256,6),rhs(256,3),C2(3,6)
      real(8) :: linear_fit(3,2),linear_rhs(3,3),initial_map(2,3),lapack_work(4096)
      real(8) :: X2(3),J1(3),J2(3),residual(3),g11,g12,g22,rhs1,rhs2,detg,step(2)
      real(8) :: g_base(32),w_base(32),lambda_base(32),P_base(32,30),projection_base(30,32)
      real(8) :: sigma(16),w_sigma(16),P_sigma(16,30),corners(2,4),Pc(2,16),Pa1(16,16),Pa2(16,16)
      real(8) :: tc_a2l(2),Jc(3,2),da(2),tc_to_target(3),delta,radius,beta_radial,eta_max
      real(8) :: C_offset(16,16,3),t0_num(2),P0(2,16),J0(3,2),offset(3),n0(3),h_use
      real(8) :: edge_a(2,6),edge_b(2,6),edge_mid(2,6),edge_dt(2,6),d_mid(2),edge_factor
      real(8) :: v_edge(3),b_edge(3),cross_edge(3),s0,beta_edge,bounds(2),half,W_leg(16,16)
      complex(8) :: z
      real(8), allocatable :: tc(:,:),g_radial(:),wg_radial(:),lambda(:),w_radial(:),T1_base(:,:)
      real(8), allocatable :: g_edge(:),wg_edge(:),s_edge(:),w_edge(:),d(:,:),q1(:),q2(:)
      real(8), allocatable :: rt_nodes(:,:),rp_nodes(:,:),normal_nodes(:,:),chord(:,:),rinv(:)
      real(8), allocatable :: q1_base(:),q2_base(:),P1_base(:,:),P2_base(:,:),base_tmp(:,:)
      real(8), allocatable :: delta_rt(:),delta_rp(:),a0(:,:),q_lambda_base(:,:,:),H0(:,:,:)
      real(8), allocatable :: kernel_S(:),kernel_D(:),omega(:),Lu(:,:),Lv(:,:),weighted_Lv(:,:)
      external :: dgels

      do iv=1,16
        do iu=1,16
          j=iu+16*(iv-1)
          B2(:,j)=[1.0d0,source_t(:,j),g16(iu)**2,g16(iu)*g16(iv),g16(iv)**2]
        enddo
      enddo
      Bfit=transpose(B2); rhs=transpose(sx-spread(xc,2,n))
      call dgels('N',256_8,6_8,3_8,Bfit,256_8,rhs,256_8,lapack_work,4096_8,info)
      if (info/=0) return
      C2=transpose(rhs(1:6,:))
      linear_fit=C2(:,2:3); linear_rhs=0.0d0
      do j=1,3
        linear_rhs(j,j)=1.0d0
      enddo
      call dgels('N',3_8,2_8,3_8,linear_fit,3_8,linear_rhs,3_8,lapack_work,4096_8,info)
      if (info/=0) return
      initial_map=linear_rhs(1:2,:)
      ninner=count(remaining); rows=pack([(j,j=1,m)],remaining)
      allocate(tc(2,ninner))
      do it=1,ninner
        residual=tx(:,rows(it))-xc-C2(:,1)
        tc(:,it)=min(1.0d0,max(-1.0d0,matmul(initial_map,residual)))
      enddo
      do iteration=1,2
        do it=1,ninner
          dx=tc(1,it); dy=tc(2,it)
          X2=matmul(C2,[1.0d0,dx,dy,dx*dx,dx*dy,dy*dy])
          J1=C2(:,2)+2.0d0*C2(:,4)*dx+C2(:,5)*dy
          J2=C2(:,3)+C2(:,5)*dx+2.0d0*C2(:,6)*dy
          residual=tx(:,rows(it))-xc-X2
          g11=sum(J1*J1); g12=sum(J1*J2); g22=sum(J2*J2)
          rhs1=sum(J1*residual); rhs2=sum(J2*residual); detg=g11*g22-g12*g12
          step=[(g22*rhs1-g12*rhs2)/detg,(g11*rhs2-g12*rhs1)/detg]
          tc(:,it)=min(1.0d0,max(-1.0d0,tc(:,it)+step))
        enddo
      enddo

      call gauss(32_8,g_base,w_base)
      lambda_base=(g_base+1.0d0)/2.0d0; w_base=w_base/2.0d0
      call legendrevec(29_8,g_base,P_base)
      do j=1,30
        projection_base(j,:)=(2*j-1)*P_base(:,j)*w_base
      enddo
      sigma=(g16+1.0d0)/2.0d0; w_sigma=wg16/2.0d0
      corners=reshape([-1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,1.0d0],[2,4])

      do it=1,ninner
        itarg=rows(it); tc_a2l=tc(:,it)
        ! Refine C inside the source square.
        C_offset=C_X
        C_offset(1,1,:)=C_offset(1,1,:)-tx(:,itarg)
        do correction=1,3
          call legendrevec(15_8,(tc_a2l-geometry_center)/geometry_scale,Pc)
          do k=1,3
            Jc(k,1)=dot_product(matmul(Pc(1,:),C_Xi(:,:,k,1)),Pc(2,:))
            Jc(k,2)=dot_product(matmul(Pc(1,:),C_Xi(:,:,k,2)),Pc(2,:))
            tc_to_target(k)=dot_product(matmul(Pc(1,:),C_offset(:,:,k)),Pc(2,:))
          enddo
          if (correction<3) then
            linear_fit=Jc; residual=tc_to_target
            call dgels('N',3_8,2_8,1_8,linear_fit,3_8,residual,3_8,lapack_work,4096_8,info)
            if (info/=0) return
            tc_a2l=min(1.0d0,max(-1.0d0,tc_a2l-residual(1:2)))
          endif
        enddo

        ! Only refine the projection within a few short-direction half-widths.
        if (norm2(tc_to_target)<=4.0d0*min(norm2(Jc(:,1)),norm2(Jc(:,2)))) then
          ! Three Gauss-Newton updates, then evaluate the final displacement.
          t0_num=tc_a2l
          do iteration=1,4
            da=t0_num-tc_a2l
            call legendrevec(15_8,(tc_a2l(1)+sigma*da(1)-geometry_center(1))/geometry_scale(1),Pa1)
            call legendrevec(15_8,(tc_a2l(2)+sigma*da(2)-geometry_center(2))/geometry_scale(2),Pa2)
            call legendrevec(15_8,(t0_num-geometry_center)/geometry_scale,P0)
            do k=1,3
              offset(k)=dot_product(w_sigma, &
                da(1)*sum(matmul(Pa1,C_Xi(:,:,k,1))*Pa2,dim=2)+da(2)*sum(matmul(Pa1,C_Xi(:,:,k,2))*Pa2,dim=2))
              J0(k,1)=dot_product(matmul(P0(1,:),C_Xi(:,:,k,1)),P0(2,:))
              J0(k,2)=dot_product(matmul(P0(1,:),C_Xi(:,:,k,2)),P0(2,:))
            enddo
            if (iteration<4) then
              linear_fit=J0; residual=tc_to_target+offset
              call dgels('N',3_8,2_8,1_8,linear_fit,3_8,residual,3_8,lapack_work,4096_8,info)
              if (info/=0) return
              t0_num=t0_num-residual(1:2)
            endif
          enddo
          n0=[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
              J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          n0=n0/norm2(n0)
          h_use=h_local(itarg)
          if (ieee_is_nan(h_use)) h_use=-dot_product(n0,tc_to_target+offset)
          residual=-offset-h_use*n0
          if (norm2(residual-tc_to_target)<=32.0d0*epsilon(1.0d0)*max(maxval(abs(C_X)),norm2(tx(:,itarg)))) then
            tc_to_target=residual
          endif
        endif
        delta=norm2(tc_to_target)
        radius=maxval(sqrt(sum(matmul(Jc,corners-spread(tc_a2l,2,4))**2,dim=1)))
        if (delta==0.0d0 .or. .not.ieee_is_finite(delta/radius)) cycle
        beta_radial=delta/radius; eta_max=asinh(1.0d0/beta_radial)
        z=cmplx(-1.0d0,pi/eta_max,kind=8)
        n_radial=max(32_8,8_8*ceiling(40.0d0/(8.0d0*abs(real(acosh(z),8))),kind=8))
        allocate(g_radial(n_radial),wg_radial(n_radial),lambda(n_radial),w_radial(n_radial),T1_base(n_radial,32))
        call gauss(n_radial,g_radial,wg_radial)
        lambda=beta_radial*sinh(eta_max*(g_radial+1.0d0)/2.0d0)
        w_radial=wg_radial*(eta_max*beta_radial/2.0d0)*cosh(eta_max*(g_radial+1.0d0)/2.0d0)
        do ir=1,n_radial
          call legendrevec(29_8,2.0d0*lambda(ir)*sigma-1.0d0,P_sigma)
          T1_base(ir,:)=matmul(matmul(w_sigma*sigma,P_sigma),projection_base)
        enddo
        edge_a=reshape([-1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0,tc_a2l(2), &
                        1.0d0,1.0d0,-1.0d0,1.0d0,-1.0d0,tc_a2l(2)],[2,6])
        edge_a=edge_a-spread(tc_a2l,2,6)
        edge_b=edge_a(:,[2,3,4,5,6,1]); edge_mid=(edge_a+edge_b)/2.0d0; edge_dt=(edge_b-edge_a)/2.0d0
        WS(itarg,:)=0.0d0; WD(itarg,:)=0.0d0
        do iedge=1,6
          d_mid=edge_mid(:,iedge)
          edge_factor=d_mid(1)*edge_dt(2,iedge)-d_mid(2)*edge_dt(1,iedge)
          if (edge_factor==0.0d0) cycle
          v_edge=matmul(Jc,edge_dt(:,iedge)); b_edge=matmul(Jc,d_mid)
          cross_edge=[v_edge(2)*b_edge(3)-v_edge(3)*b_edge(2),v_edge(3)*b_edge(1)-v_edge(1)*b_edge(3), &
                      v_edge(1)*b_edge(2)-v_edge(2)*b_edge(1)]
          s0=-dot_product(v_edge,b_edge)/sum(v_edge*v_edge)
          beta_edge=norm2(cross_edge)/sum(v_edge*v_edge)
          bounds=asinh(([-1.0d0,1.0d0]-s0)/beta_edge); half=(bounds(2)-bounds(1))/2.0d0
          z=cmplx(-sum(bounds)/2.0d0,pi/2.0d0,kind=8)/half
          n_edge=max(32_8,8_8*ceiling(40.0d0/(8.0d0*abs(real(acosh(z),8))),kind=8))
          nquad=n_radial*n_edge
          allocate(g_edge(n_edge),wg_edge(n_edge),s_edge(n_edge),w_edge(n_edge),d(2,n_edge),q1(nquad),q2(nquad))
          allocate(Pu(nquad,16),Pv(nquad,16),Xtmp(nquad,16),rt_nodes(nquad,3),rp_nodes(nquad,3))
          allocate(normal_nodes(nquad,3),chord(nquad,3),rinv(nquad),omega(nquad),kernel_S(nquad),kernel_D(nquad))
          allocate(q1_base(32*n_edge),q2_base(32*n_edge),P1_base(32*n_edge,16),P2_base(32*n_edge,16))
          allocate(base_tmp(32*n_edge,16),delta_rt(32*n_edge),delta_rp(32*n_edge),a0(3,n_edge))
          allocate(q_lambda_base(32,n_edge,3),H0(n_radial,n_edge,3),Lu(nquad,16),Lv(nquad,16),weighted_Lv(nquad,16))
          call gauss(n_edge,g_edge,wg_edge)
          s_edge=beta_edge*sinh(bounds(1)+(g_edge+1.0d0)*half)
          w_edge=wg_edge*(half*beta_edge)*cosh(bounds(1)+(g_edge+1.0d0)*half)
          do j=1,2
            d(j,:)=(d_mid(j)+edge_dt(j,iedge)*s0)+edge_dt(j,iedge)*s_edge
          enddo
          q1=reshape(tc_a2l(1)+spread(lambda,2,n_edge)*spread(d(1,:),1,n_radial),[nquad])
          q2=reshape(tc_a2l(2)+spread(lambda,2,n_edge)*spread(d(2,:),1,n_radial),[nquad])
          call legendrevec(15_8,(q1-geometry_center(1))/geometry_scale(1),Pu)
          call legendrevec(15_8,(q2-geometry_center(2))/geometry_scale(2),Pv)
          do k=1,3
            call dgemm('N','N',nquad,16_8,16_8,1.0d0,Pu,nquad,C_Xi(:,:,k,1),16_8,0.0d0,Xtmp,nquad)
            rt_nodes(:,k)=sum(Xtmp*Pv,dim=2)
            call dgemm('N','N',nquad,16_8,16_8,1.0d0,Pu,nquad,C_Xi(:,:,k,2),16_8,0.0d0,Xtmp,nquad)
            rp_nodes(:,k)=sum(Xtmp*Pv,dim=2)
          enddo
          normal_nodes(:,1)=rp_nodes(:,2)*rt_nodes(:,3)-rp_nodes(:,3)*rt_nodes(:,2)
          normal_nodes(:,2)=rp_nodes(:,3)*rt_nodes(:,1)-rp_nodes(:,1)*rt_nodes(:,3)
          normal_nodes(:,3)=rp_nodes(:,1)*rt_nodes(:,2)-rp_nodes(:,2)*rt_nodes(:,1)

          ! Same near/self radial H0 construction, using the supplied tangent fields.
          q1_base=reshape(tc_a2l(1)+spread(lambda_base,2,n_edge)*spread(d(1,:),1,32),[32*n_edge])
          q2_base=reshape(tc_a2l(2)+spread(lambda_base,2,n_edge)*spread(d(2,:),1,32),[32*n_edge])
          call legendrevec(15_8,(q1_base-geometry_center(1))/geometry_scale(1),P1_base)
          call legendrevec(15_8,(q2_base-geometry_center(2))/geometry_scale(2),P2_base)
          a0=matmul(Jc,d)
          do k=1,3
            call dgemm('N','N',32*n_edge,16_8,16_8,1.0d0,P1_base,32*n_edge, &
              C_Xi(:,:,k,1),16_8,0.0d0,base_tmp,32*n_edge)
            delta_rt=sum(base_tmp*P2_base,dim=2)-Jc(k,1)
            call dgemm('N','N',32*n_edge,16_8,16_8,1.0d0,P1_base,32*n_edge, &
              C_Xi(:,:,k,2),16_8,0.0d0,base_tmp,32*n_edge)
            delta_rp=sum(base_tmp*P2_base,dim=2)-Jc(k,2)
            q_lambda_base(:,:,k)=(reshape(delta_rt,[32,n_edge])*spread(d(1,:),1,32) &
              +reshape(delta_rp,[32,n_edge])*spread(d(2,:),1,32))/spread(lambda_base,2,n_edge)
          enddo
          call dgemm('N','N',n_radial,3*n_edge,32_8,1.0d0,T1_base,n_radial,q_lambda_base,32_8,0.0d0,H0,n_radial)
          do k=1,3
            chord(:,k)=reshape(spread(lambda,2,n_edge)*spread(a0(k,:),1,n_radial) &
              +spread(lambda**2,2,n_edge)*H0(:,:,k),[nquad])+tc_to_target(k)
          enddo
          rinv=1.0d0/sqrt(sum(chord*chord,dim=2))
          omega=reshape(spread(lambda*w_radial,2,n_edge)*spread(edge_factor*w_edge,1,n_radial),[nquad])*inv2pi
          kernel_S=omega*sqrt(sum(normal_nodes*normal_nodes,dim=2))*rinv
          kernel_D=omega*sum(normal_nodes*chord,dim=2)*rinv**3

          ! Both kernels act on the same GL16 density interpolation basis.
          call legendrevec(15_8,q1,Pu)
          call legendrevec(15_8,q2,Pv)
          call dgemm('N','T',nquad,16_8,16_8,1.0d0,Pu,nquad,density_map,16_8,0.0d0,Lu,nquad)
          call dgemm('N','T',nquad,16_8,16_8,1.0d0,Pv,nquad,density_map,16_8,0.0d0,Lv,nquad)
          weighted_Lv=Lv*spread(kernel_S,2,16)
          call dgemm('T','N',16_8,16_8,nquad,1.0d0,Lu,nquad,weighted_Lv,nquad,0.0d0,W_leg,16_8)
          WS(itarg,:)=WS(itarg,:)+reshape(W_leg,[256])
          weighted_Lv=Lv*spread(kernel_D,2,16)
          call dgemm('T','N',16_8,16_8,nquad,1.0d0,Lu,nquad,weighted_Lv,nquad,0.0d0,W_leg,16_8)
          WD(itarg,:)=WD(itarg,:)+reshape(W_leg,[256])
          deallocate(g_edge,wg_edge,s_edge,w_edge,d,q1,q2,Pu,Pv,Xtmp,rt_nodes,rp_nodes,normal_nodes,chord,rinv)
          deallocate(omega,kernel_S,kernel_D,q1_base,q2_base,P1_base,P2_base,base_tmp,delta_rt,delta_rp,a0)
          deallocate(q_lambda_base,H0,Lu,Lv,weighted_Lv)
        enddo
        deallocate(g_radial,wg_radial,lambda,w_radial,T1_base)
      enddo
    end block
  else
    where (remaining) ier=1_8
  endif

end subroutine quadSD_mid_w16


! Helmholtz version of quadSD_mid_w16: exp(1i*zk*r)/r, any complex zk, in the 1/(2*pi) normalization, with
! the same regions, nodes and contractions. At a node with d = y-x, r = |d|, z = zk*r, e = exp(1i*z),
! c = e*(1-1i*z), SLP value v and f = v/r^2: SLP v*e, DLP f*(d.n)*c, SLPn_c f*d_c*c,
! DLPn_ab f*(c*(3*d_a*d_b/r^2-delta_ab)-zk^2*e*d_a*d_b), evaluated per batch of nodes in compute_helm_kernels.
! Real and imaginary parts are contracted as separate real rows (regions 1:6: n_row = 2*n_kernel, kernel k
! in row block k, its imaginary part in block n_kernel+k; Region 7: two contractions per kernel).
! WS, WD, WSn, WDn are complex.
subroutine quadSD_helm_mid_w16(m,tx,tnx,n,sx,snx,sw,WS,WD,ifSn,ifDn,WSn,WDn, &
    t0,h,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification, &
    C_Xi_typeb,C_Xij_typeb,geometry_center_typeb,geometry_scale_typeb,ier,root_ok,zk)
  use surfquadr_mid, only: legendrevec, gauss => compute_gaussleg
  use, intrinsic :: iso_c_binding, only: c_double,c_bool,c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  implicit none
  integer(8), intent(in) :: m,n
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),t0(2,m),h(m)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3)
  real(8), intent(in) :: geometry_center(2),geometry_scale(2)
  ! *_typeb: geometry of the type-b chart the source sits in; geometry_scale_typeb = 0: no type-b chart.
  real(8), intent(in) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3)
  real(8), intent(in) :: geometry_center_typeb(2),geometry_scale_typeb(2)
  integer(8), intent(in) :: ifSn,ifDn
  complex(8), intent(inout) :: WS(m,n),WD(m,n),WSn(m,n*ifSn,3),WDn(m,n*ifDn,3,3)
  integer(8), intent(inout) :: classification(m),ier(m)
  logical(c_bool), intent(in) :: root_ok(m)
  real(8) :: t0_local(2,m),h_local(m)
  ! The 256 source nodes are ordered with u varying fastest.
  ! The caller excludes near/self and supplies the Region 7 projection and signed height.
  ! Regions 1:3 use 1D parity factors for every kernel.
  ! At a node with d = y-x, r = |d| and SLP value v, SLPn_c = (v/r^2)*d_c and
  ! DLPn_ab = (v/r^2)*(3*d_a*d_b/r^2-delta_ab), DLPn components (a,b) = (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
  ! Where the kernels share one contraction (regions 1:6) they are stacked by kernel index k: 1 SLP, 2 DLP,
  ! k_Sn+kn SLPn kn, k_Dn+kn DLPn kn; n_kernel = 2, 5, 8 or 11.
  integer(8), parameter :: i_Dn(6)=[1_8,1_8,1_8,2_8,2_8,3_8], j_Dn(6)=[1_8,2_8,3_8,2_8,3_8,3_8]
  integer(8) :: n_kernel,k_Sn,k_Dn,kn,a_c,b_c,q_node
  real(8) :: delta_c
  ! SLPn, DLPn: displacement y-x, 1/r^2 and v/r^2 at the nodes of one evaluation step
  real(8), allocatable :: node_d(:,:,:),node_r2(:,:),node_f(:,:),W_k(:,:)
  integer(8) :: i,j,k,m16,mr,nq,nh,ld,layer,info,piv(16)
  ! integer(8), parameter :: order_p(3)=[32,64,128]
  integer(8), parameter :: order_p(3)=[32,48,96]
  integer(8), allocatable :: idx16(:),rows(:)
  real(8) :: g16(16),wg16(16),V16(16,16),Vsolve(16,16),density_map(16,16),source_t(2,256)
  real(8) :: xc(3),rt0(3),rp0(3),e1(3),e2(3),e3(3),frame(3,3),r(3),xyz(3,m)
  real(8) :: ratio
  real(8) :: lambda1,d16,d32,d64,d128,abc32(3),abc64(3),abc128(3)
  real(8), parameter :: pi=acos(-1.0d0),inv2pi=0.5d0/pi
  logical(c_bool) :: outside16(m),inside16(m)
  logical :: inside32(m),inside64(m),inside128(m),shell32(m),shell64(m),shell128(m)
  logical :: target_groups(m,3),remaining(m),use_sinh(m),off_column4(m)
  real(8), allocatable :: Pu(:,:),Pv(:,:),X(:,:),rt(:,:,:),rp(:,:,:),Xtmp(:,:)
  real(8), allocatable :: g(:),wg(:),L(:,:,:),sq_x(:,:),sq_nx(:,:),sq_w(:),target(:,:),Y(:,:)
  real(8) :: P_half(maxval(order_p)/2,16),L_half(maxval(order_p)/2,16)
  real(8) :: dx,dy,dz,distance2_inv,distance_inv,kernel_value
  integer(8) :: iu,iv,mode,iq,is,node(4,4096),column(4,64)
  real(8) :: v1,v2,v3,v4,p1,p2,p3,p4,d1,d2,d3,d4
  real(8), allocatable :: parity_K(:,:,:),parity_W(:,:,:)
  external :: dgemm,dgesv
  ! Helmholtz: wavenumber; Laplace values, chord, 1/r, r and the Helmholtz kernels (Re, Im columns) at the
  ! nodes of one evaluation step (compute_helm_kernels).
  complex(8), intent(in) :: zk
  integer(8) :: n_row,part,row0,kk,jn,ie,ln,off_node,base_node
  real(8) :: w_node
  real(8), allocatable :: v_node(:),dlp_node(:),d_node(:,:),rinv_node(:),r_node(:),K_val(:,:),W_part(:,:,:)
  complex(8), allocatable :: W_c(:,:)

#ifdef SURFQUADR_USE_VFORCE
  interface
    subroutine vvsincos(sine,cosine,phase,count) bind(C,name="vvsincos")
      import c_double,c_int
      real(c_double), intent(out) :: sine(*),cosine(*)
      real(c_double), intent(in) :: phase(*)
      integer(c_int), intent(in) :: count
    end subroutine vvsincos
    subroutine vvexp(output,input,count) bind(C,name="vvexp")
      import c_double,c_int
      real(c_double), intent(out) :: output(*)
      real(c_double), intent(in) :: input(*)
      integer(c_int), intent(in) :: count
    end subroutine vvexp
  end interface
#endif

  ier=0_8
  n_kernel=2_8+3_8*ifSn+6_8*ifDn; k_Sn=2_8; k_Dn=2_8+3_8*ifSn; n_row=2_8*n_kernel
  call prepare_quadSD_targets_w16(m,tx,n,sx,C_X,C_Xi,geometry_center,geometry_scale, &
      inside16,outside16,xc,rt0,rp0,e1,e2,e3,r,xyz,d16,d128,lambda1)
  frame(:,1)=e1; frame(:,2)=e2; frame(:,3)=e3

  d32=0.50d0*d16; d64=0.25d0*d16
  abc32=lambda1*sqrt(r*r+d32*d32)
  abc64=lambda1*sqrt(r*r+d64*d64); abc128=lambda1*sqrt(r*r+d128*d128)

  ! 1st layer of guard.
  where (outside16) classification=0_8
  m16=count(outside16); idx16=pack([(k,k=1,m)],outside16)
  if (m16>0) then
    target=tx(:,idx16)
    allocate(v_node(m16),dlp_node(m16),d_node(m16,3),rinv_node(m16),r_node(m16),K_val(m16,n_row))
    do j=1,n
      !GCC$ ivdep
      do i=1,m16
        dx=sx(1,j)-target(1,i); dy=sx(2,j)-target(2,i); dz=sx(3,j)-target(3,i)
        d_node(i,1)=dx; d_node(i,2)=dy; d_node(i,3)=dz
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz); rinv_node(i)=sqrt(distance2_inv)
        r_node(i)=(dx*dx+dy*dy+dz*dz)*rinv_node(i)
        v_node(i)=(sw(j)*inv2pi)*rinv_node(i)
        dlp_node(i)=v_node(i)*(dx*snx(1,j)+dy*snx(2,j)+dz*snx(3,j))*distance2_inv
      enddo
      call compute_helm_kernels(m16,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
      WS(idx16,j)=cmplx(K_val(:,1),K_val(:,n_kernel+1),kind=8)
      WD(idx16,j)=cmplx(K_val(:,2),K_val(:,n_kernel+2),kind=8)
      do kn=1,3*ifSn
        WSn(idx16,j,kn)=cmplx(K_val(:,k_Sn+kn),K_val(:,n_kernel+k_Sn+kn),kind=8)
      enddo
      do kn=1,6*ifDn
        WDn(idx16,j,i_Dn(kn),j_Dn(kn))=cmplx(K_val(:,k_Dn+kn),K_val(:,n_kernel+k_Dn+kn),kind=8)
        WDn(idx16,j,j_Dn(kn),i_Dn(kn))=WDn(idx16,j,i_Dn(kn),j_Dn(kn))
      enddo
    enddo
    deallocate(target,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
  endif
  deallocate(idx16)

  inside32=inside16 .and. sum((xyz/spread(abc32,2,m))**2,dim=1)<=1.0d0
  shell32=inside16 .and. .not.inside32
  inside64=inside32 .and. sum((xyz/spread(abc64,2,m))**2,dim=1)<=1.0d0
  shell64=inside32 .and. .not.inside64
  inside128=inside64 .and. sum((xyz/spread(abc128,2,m))**2,dim=1)<=1.0d0
  shell128=inside64 .and. .not.inside128
  target_groups(:,1)=shell32; target_groups(:,2)=shell64; target_groups(:,3)=shell128

  if (any(inside16)) then
    call gauss(16_8,g16,wg16)
    call legendrevec(15_8,g16,V16)
    Vsolve=transpose(V16); density_map=0.0d0
    do j=1,16
      density_map(j,j)=1.0d0
    enddo
    call dgesv(16_8,16_8,Vsolve,16_8,piv,density_map,16_8,info)
    if (info/=0) return
  endif

  ! Common p refinement: contract v, then u for order_p.
  do layer=1,3
    where (target_groups(:,layer)) classification=layer
    mr=count(target_groups(:,layer))
    if (mr==0) cycle
    rows=pack([(k,k=1,m)],target_groups(:,layer))
    ! nq=16*2**layer; nh=nq/2
    nq=order_p(layer); nh=nq/2
    allocate(g(nq),wg(nq),L(nh,8,2))
    call gauss(nq,g,wg)
    call legendrevec(15_8,g(:nh),P_half(:nh,:))
    call dgemm('N','T',nh,16_8,16_8,1d0,P_half,size(P_half,1,kind=8), &
      density_map,16_8,0d0,L_half,size(L_half,1,kind=8))
    L(:,:,1)=0.5d0*(L_half(:nh,1:8)+L_half(:nh,16:9:-1))
    L(:,:,2)=0.5d0*(L_half(:nh,1:8)-L_half(:nh,16:9:-1))
    allocate(Pu(nq,16),Pv(nq,16),X(nq,nq),rt(nq,nq,3),rp(nq,nq,3),Xtmp(nq,16))
    allocate(sq_x(3,nq*nq),sq_nx(3,nq*nq),sq_w(nq*nq))
    call legendrevec(15_8,(g-geometry_center(1))/geometry_scale(1),Pu)
    call legendrevec(15_8,(g-geometry_center(2))/geometry_scale(2),Pv)
    do k=1,3
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_X(:,:,k),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,X,nq)
      sq_x(k,:)=reshape(X,[nq*nq])
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,1),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rt(:,:,k),nq)
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,2),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rp(:,:,k),nq)
    enddo
    sq_nx(1,:)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[nq*nq])
    sq_nx(2,:)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[nq*nq])
    sq_nx(3,:)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[nq*nq])
    sq_w=sqrt(sum(sq_nx*sq_nx,dim=1)); sq_nx=sq_nx/spread(sq_w,1,3)
    sq_w=sq_w*reshape(spread(wg,2,nq)*spread(wg,1,nq),[nq*nq])
    target=transpose(tx(:,rows))
    ! Keep targets contiguous; contract v, then u without transposing kernels.
    ! Row i+mr*(k-1): target i, kernel k.
    ! Helmholtz: row i+mr*(k-1) holds the real part of kernel k, row i+mr*(n_kernel+k-1) its imaginary part.
    ld=8*((n_row*mr+7)/8)
    allocate(parity_K(ld,nh*nh,4),parity_W(ld,64,4),Y(ld*nh,8))
    allocate(v_node(4*mr),dlp_node(4*mr),d_node(4*mr,3),rinv_node(4*mr),r_node(4*mr),K_val(4*mr,n_row))
    parity_K(n_row*mr+1:ld,:,:)=0.0d0
    do iv=1,nh
      do iu=1,nh
        iq=iu+nh*(iv-1)
        node(:,iq)=[iu+nq*(iv-1),nq+1-iu+nq*(iv-1),iu+nq*(nq-iv),nq+1-iu+nq*(nq-iv)]
      enddo
    enddo
    do iv=1,8
      do iu=1,8
        is=iu+8*(iv-1)
        column(:,is)=[iu+16*(iv-1),17-iu+16*(iv-1),iu+16*(16-iv),17-iu+16*(16-iv)]
      enddo
    enddo
    do j=1,nh*nh
      ! The four mirrored nodes (entries i+mr*(q_node-1)), their kernel values, then the parity sums.
      do q_node=1,4
        jn=node(q_node,j); w_node=sq_w(jn)*inv2pi
        !GCC$ ivdep
        do i=1,mr
          dx=sq_x(1,jn)-target(i,1); dy=sq_x(2,jn)-target(i,2); dz=sq_x(3,jn)-target(i,3)
          ie=i+mr*(q_node-1)
          d_node(ie,1)=dx; d_node(ie,2)=dy; d_node(ie,3)=dz
          distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz); rinv_node(ie)=sqrt(distance2_inv)
          r_node(ie)=(dx*dx+dy*dy+dz*dz)*rinv_node(ie)
          v_node(ie)=w_node*rinv_node(ie)
          dlp_node(ie)=v_node(ie)*(dx*sq_nx(1,jn)+dy*sq_nx(2,jn)+dz*sq_nx(3,jn))*distance2_inv
        enddo
      enddo
      call compute_helm_kernels(4*mr,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
      do kk=1,n_row
        row0=mr*(kk-1)
        !GCC$ ivdep
        do i=1,mr
          p1=K_val(i,kk)+K_val(mr+i,kk); p2=K_val(i,kk)-K_val(mr+i,kk)
          p3=K_val(2*mr+i,kk)+K_val(3*mr+i,kk); p4=K_val(2*mr+i,kk)-K_val(3*mr+i,kk)
          parity_K(row0+i,j,1)=p1+p3; parity_K(row0+i,j,2)=p2+p4
          parity_K(row0+i,j,3)=p1-p3; parity_K(row0+i,j,4)=p2-p4
        enddo
      enddo
    enddo
    do mode=1,4
      call dgemm('N','N',ld*nh,8_8,nh,1.0d0,parity_K(1,1,mode),ld*nh, &
        L(1,1,1+(mode-1)/2),nh,0.0d0,Y,ld*nh)
      do iv=1,8
        call dgemm('N','N',n_row*mr,8_8,nh,1.0d0,Y(1,iv),ld,L(1,1,1+mod(mode-1,2)),nh, &
          0.0d0,parity_W(1,1+8*(iv-1),mode),ld)
      enddo
    enddo
    allocate(W_part(mr,n,2))
    do k=1,n_kernel
      do part=1,2
        row0=mr*(k-1)+mr*n_kernel*(part-1)
        do j=1,64
          !GCC$ ivdep
          do i=1,mr
            v1=parity_W(row0+i,j,1); v2=parity_W(row0+i,j,2)
            v3=parity_W(row0+i,j,3); v4=parity_W(row0+i,j,4)
            p1=v1+v2; p2=v3+v4; p3=v1-v2; p4=v3-v4
            W_part(i,column(1,j),part)=p1+p2; W_part(i,column(2,j),part)=p3+p4
            W_part(i,column(3,j),part)=p1-p2; W_part(i,column(4,j),part)=p3-p4
          enddo
        enddo
      enddo
      W_c=cmplx(W_part(:,:,1),W_part(:,:,2),kind=8)
      if (k==1) then
        WS(rows,:)=W_c
      else if (k==2) then
        WD(rows,:)=W_c
      else if (k<=k_Dn) then
        WSn(rows,:,k-k_Sn)=W_c
      else
        WDn(rows,:,i_Dn(k-k_Dn),j_Dn(k-k_Dn))=W_c; WDn(rows,:,j_Dn(k-k_Dn),i_Dn(k-k_Dn))=W_c
      endif
    enddo
    deallocate(W_part,W_c,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
    deallocate(parity_K,parity_W,Y)
    deallocate(Pu,Pv,X,Xtmp,rt,rp,sq_x,sq_nx,sq_w,target)
    deallocate(g,wg,L,rows)
  enddo

  remaining=inside128
  off_column4=.false.
  t0_local=t0; h_local=h
  if (any(inside128)) then
    do iv=1,16
      do iu=1,16
        source_t(:,iu+16*(iv-1))=[g16(iu),g16(iv)]
      enddo
    enddo
  endif

  ! Regions 4:6: panel refinement for patches up to 10:1.
  ratio = norm2(rp0)/norm2(rt0)
  ratio = max(ratio,1d0/ratio)
  ! if (any(remaining) .and. ratio <= 10d0) then
  if (any(remaining) .and. ratio <= 8d0) then ! heuristic
    block
      integer(8), parameter :: order4(2)=[64,16] ! Common orders for regions 4:6.
      integer(8), parameter :: nv=order4(2),npq=product(order4)
      integer(8), parameter :: preset=2 ! 1: five rows; 2: three rows.
      integer(8) :: p,ip,ip_right,pr,direction,n_panel,np,cols,parameter_order(2),density_column(16,16)
      integer(8) :: panel_count(3),row,nparent,nchild,nsplit,nlayout,first_panel,last_panel
      integer(8), allocatable :: ncols(:),panel_id(:,:),panel_row(:),needed(:),panel_mirror(:),panel_pair(:,:,:)
      integer(8) :: node_order(npq),reverse_nodes(npq)
      real(8) :: box_center(3),box_frame(3,3),box_halfsize4(3),box_halfsize5(3),box_halfsize6(3)
      real(8) :: scale_box4,scale_box5,scale_box6,d_box4,d_box5,d_box6
      real(8), allocatable :: row_ratio(:),split_ratio(:),v_edges(:),panels_all(:,:)
      real(8) :: width,top_height
      real(8) :: b(4),c(2),half(2),rotation(2,2,4)
      real(8) :: tu(order4(1)),tv(nv),g_h(order4(1)),wg_h(order4(1)),g_v(nv),wg_v(nv),wq(npq)
      real(8) :: C_Xp(16,16,3),C_Xip(16,16,3,2)
      real(8) :: normal_h(npq,3),weight_h(npq)
      real(8), allocatable :: Lu(:,:,:),Lv(:,:,:),W_mid(:,:)
      logical :: inside_box4(m),inside_box5(m),inside_box6(m),region(m,3),direction_groups(m,4),active(3)
      logical, allocatable :: use_panel(:)
      ! Cylindrical bounds around the self patch for regions 4:6.
      integer(8), parameter :: p_R=3,n_R=(p_R+1)*(p_R+2)/2,n_b=68
      integer(8) :: i_R(n_R),j_R(n_R),iR,jR,kR,itarg
      real(8) :: m_box(3),H_box,gram(3),t_b(n_b,2),P_bu(n_b,16),P_bv(n_b,16),x_b(3,n_b),xyz_s(3,n+n_b),s_R(2)
      real(8) :: B_s(n+n_b,n_R),B_fit(n,n_R),c_R(n,2),r_s(2),zeta_lo,zeta_hi,ell_R(2),b_box(3,2),lapack_work(4096)
      real(8), allocatable :: P_a(:,:),P_b(:,:),B_t(:,:),R_t(:,:)
      real(8) :: P_est(2,16),X_est(3),J_est(3,2),n_est(3),t_direction(2,m)
      logical :: in_slab(m),beside(m)
      external :: dgels

      box_center=xc; box_frame=frame
      t_direction=t0_local
      if (.false.) then ! old rule, turns out to be not valid for sphere, near triple junction
        ! scale_box4=0.95d0; d_box4=0.15d0*maxval(r(1:2))
        scale_box4=0.95d0; d_box4=0.0425d0*maxval(r(1:2))
        ! scale_box5=0.95d0; d_box5=0.080d0*maxval(r(1:2))
        scale_box5=0.65d0; d_box5=0.0d0*maxval(r(1:2))
        ! scale_box6=0.95d0; d_box6=0.0d0*maxval(r(1:2))
        ! scale_box6=0.90d0; d_box6=0.0d0*maxval(r(1:2))
        scale_box6=0.60d0; d_box6=0.0d0*maxval(r(1:2))
        box_halfsize4=scale_box4*r+d_box4
        box_halfsize5=scale_box5*r+d_box5
        box_halfsize6=scale_box6*r+d_box6
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        inside_box4=all(abs(xyz)<=spread(box_halfsize4,2,m),dim=1)
        inside_box5=all(abs(xyz)<=spread(box_halfsize5,2,m),dim=1)
        inside_box6=all(abs(xyz)<=spread(box_halfsize6,2,m),dim=1)
      else
        ! Boxes 4:6 as cylinders around the self patch.
        ! (xi,eta,zeta) = box_frame^T (x-box_center); (u,v) = R(xi,eta): least-squares inverse fit of degree p_R
        ! on the source nodes. box_k = {|R_i| <= r_s,i + m_k*l_other/l_i, zeta_lo-H_box <= zeta <= zeta_hi+H_box};
        ! r_s = extent of R over self, so self is inside every box.
        ! Region 7 = box6; a target with trusted (t0,h) must also satisfy |t0_i| <= b6_i
        ! (keeps targets beside a curved edge out of Region 7).
        !
        ! things to tweak inside, heuristic: p_R (parameter above), m_box, H_box
        m_box=[0.10d0,0.06d0,0.02d0]
        gram=[dot_product(rt0,rt0),dot_product(rt0,rp0),dot_product(rp0,rp0)]
        H_box=0.6d0*sqrt(0.5d0*(gram(1)+gram(3))+sqrt(0.25d0*(gram(1)-gram(3))**2+gram(2)**2)) ! 0.6*||[rt0,rp0]||_2
        !
        kR=0 ! total degree <= p_R
        do jR=1,p_R+1
          do iR=1,p_R+2-jR
            kR=kR+1; i_R(kR)=iR; j_R(kR)=jR
          enddo
        enddo
        t_b(1:16,1)=g16; t_b(1:16,2)=-1d0 ! self boundary: GL16 on each edge + corners
        t_b(17:32,1)=1d0; t_b(17:32,2)=g16
        t_b(33:48,1)=g16; t_b(33:48,2)=1d0
        t_b(49:64,1)=-1d0; t_b(49:64,2)=g16
        t_b(65:68,1)=[-1d0,1d0,1d0,-1d0]; t_b(65:68,2)=[-1d0,-1d0,1d0,1d0]
        call legendrevec(15_8,(t_b(:,1)-geometry_center(1))/geometry_scale(1),P_bu)
        call legendrevec(15_8,(t_b(:,2)-geometry_center(2))/geometry_scale(2),P_bv)
        do k=1,3
          x_b(k,:)=sum(matmul(P_bu,C_X(:,:,k))*P_bv,dim=2)
        enddo
        xyz_s(:,1:n)=matmul(transpose(box_frame),sx-spread(box_center,2,n)) ! self samples: n nodes, then boundary
        xyz_s(:,n+1:)=matmul(transpose(box_frame),x_b-spread(box_center,2,n_b))
        s_R=maxval(abs(xyz_s(1:2,1:n)),dim=2)
        allocate(P_a(n+n_b,p_R+1),P_b(n+n_b,p_R+1))
        call legendrevec(p_R,xyz_s(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz_s(2,:)/s_R(2),P_b)
        B_s=P_a(:,i_R)*P_b(:,j_R)
        deallocate(P_a,P_b)
        B_fit=B_s(1:n,:); c_R=transpose(source_t) ! fit on the n source nodes
        call dgels('N',n,n_R,2_8,B_fit,n,c_R,n,lapack_work,4096_8,info)
        r_s=maxval(abs(matmul(B_s,c_R(1:n_R,:))),dim=1) ! self extent of R (nodes and boundary)
        zeta_lo=minval(xyz_s(3,:)); zeta_hi=maxval(xyz_s(3,:))
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        allocate(P_a(m,p_R+1),P_b(m,p_R+1),B_t(m,n_R),R_t(m,2))
        call legendrevec(p_R,xyz(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz(2,:)/s_R(2),P_b)
        B_t=P_a(:,i_R)*P_b(:,j_R)
        R_t=matmul(B_t,c_R(1:n_R,:))
        do itarg=1,m
          if (.not.root_ok(itarg) .or. .not.all(ieee_is_finite(t0_local(:,itarg)))) t_direction(:,itarg)=R_t(itarg,:)
        enddo
        in_slab=xyz(3,:)>=zeta_lo-H_box .and. xyz(3,:)<=zeta_hi+H_box
        ! facing-edge guard: the margin across an edge scales with that edge's length (anisotropic patches)
        ell_R=[norm2(rt0),norm2(rp0)] ! physical length per unit t1, t2
        b_box(:,1)=r_s(1)+m_box*ell_R(2)/ell_R(1) ! rows box4/5/6, columns t1/t2
        b_box(:,2)=r_s(2)+m_box*ell_R(1)/ell_R(2)
        ! Off-surface targets in column4 bypass the slab and regions 4:6.
        off_column4=remaining .and. abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. h_local/=0d0
        do itarg=1,m
          if (.not.off_column4(itarg)) cycle
          if (all(ieee_is_finite(t0_local(:,itarg))) .and. ieee_is_finite(h_local(itarg))) cycle
          t0_local(:,itarg)=R_t(itarg,:)
          call legendrevec(15_8,(t0_local(:,itarg)-geometry_center)/geometry_scale,P_est)
          do k=1,3
            X_est(k)=dot_product(matmul(P_est(1,:),C_X(:,:,k)),P_est(2,:))
            J_est(k,1)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,1)),P_est(2,:))
            J_est(k,2)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,2)),P_est(2,:))
          enddo
          n_est=[J_est(2,2)*J_est(3,1)-J_est(3,2)*J_est(2,1), &
                 J_est(3,2)*J_est(1,1)-J_est(1,2)*J_est(3,1), &
                 J_est(1,2)*J_est(2,1)-J_est(2,2)*J_est(1,1)]
          n_est=n_est/norm2(n_est)
          h_local(itarg)=dot_product(tx(:,itarg)-X_est,n_est)
        enddo
        off_column4=off_column4 .and. ieee_is_finite(h_local) .and. all(ieee_is_finite(t0_local),dim=1)
        remaining=remaining .and. .not.off_column4
        inside_box4=abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. in_slab
        inside_box5=abs(R_t(:,1))<=b_box(2,1) .and. abs(R_t(:,2))<=b_box(2,2) .and. in_slab
        inside_box6=abs(R_t(:,1))<=b_box(3,1) .and. abs(R_t(:,2))<=b_box(3,2) .and. in_slab
        ! Region 7 only above the patch: a target with trusted t0,h must also have its root inside box6
        beside=root_ok .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local) .and. &
               (abs(t0_local(1,:))>b_box(3,1) .or. abs(t0_local(2,:))>b_box(3,2))
        inside_box6=inside_box6 .and. .not.beside
        deallocate(P_a,P_b,B_t,R_t)
      endif
      region(:,1)=remaining .and. .not.inside_box4; remaining=remaining .and. inside_box4
      region(:,2)=remaining .and. .not.inside_box5; remaining=remaining .and. inside_box5
      region(:,3)=remaining .and. .not.inside_box6; remaining=remaining .and. inside_box6
      where (region(:,1)) classification=4_8
      where (region(:,2)) classification=5_8
      where (region(:,3)) classification=6_8

      ! List row heights and even column counts from the target-facing edge down.
      select case (preset)
      case (1)
        row_ratio=[1d0,1d0,1d0,3d0,3d0]
        ncols=[4,4,4,2,2]
        split_ratio=[1d0,1d0,1d0]
      case (2)
        row_ratio=[1d0,2d0,6d0]
        ncols=[4,4,2]
        split_ratio=[1d0,2d0]
      end select
      nsplit=size(split_ratio)
      panel_count(1)=sum(ncols); nlayout=panel_count(1); nparent=ncols(1)
      do layer=2,3
        nchild=2*nparent*nsplit
        panel_count(layer)=panel_count(layer-1)-nparent+nchild
        nlayout=nlayout+nchild; nparent=2*nparent
      enddo
      allocate(panels_all(4,nlayout),panel_id(maxval(panel_count),3),panel_mirror(nlayout))
      allocate(panel_pair(2,maxval(panel_count)/2,3),panel_row(nlayout),needed(nlayout),use_panel(nlayout))
      v_edges=[(1d0-2d0*sum(row_ratio(:row))/sum(row_ratio),row=0,size(row_ratio))]
      last_panel=0
      do layer=1,3
        if (layer>1) then
          nparent=ncols(1); top_height=1d0-v_edges(2)
          v_edges=[(1d0-top_height*sum(split_ratio(:row))/sum(split_ratio),row=0,nsplit)]
          ncols=spread(2*nparent,1,nsplit)
        endif
        first_panel=last_panel+1
        do row=1,size(ncols)
          width=2d0/ncols(row)
          do j=1,ncols(row)
            p=last_panel+j
            panels_all(:,p)=[-1d0+(j-1)*width,-1d0+j*width,v_edges(row+1),v_edges(row)]
            panel_mirror(p)=last_panel+ncols(row)+1-j
          enddo
          last_panel=last_panel+ncols(row)
        enddo
        nchild=last_panel-first_panel+1
        panel_id(1:nchild,layer)=[(p,p=first_panel,last_panel)]
        ! Lower rows retain their IDs from the previous layout.
        if (layer>1) panel_id(nchild+1:panel_count(layer),layer)=panel_id(nparent+1:panel_count(layer-1),layer-1)
      enddo
      ! Pair canonical left/right rectangles; geometry and targets need no symmetry.
      do layer=1,3
        is=0
        do p=1,panel_count(layer)
          ip=panel_id(p,layer)
          if (panels_all(1,ip)+panels_all(2,ip)>=0d0) cycle
          is=is+1
          panel_pair(:,is,layer)=[p,findloc(panel_id(1:panel_count(layer),layer),panel_mirror(ip),dim=1)]
        enddo
      enddo
      direction_groups(:,1)=t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,2)=t_direction(1,:)>abs(t_direction(2,:))
      direction_groups(:,3)=t_direction(2,:)<0 .and. -t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,4)=-t_direction(1,:)>abs(t_direction(2,:))
      rotation(:,:,1)=reshape([1d0,0d0,0d0,1d0],[2,2])
      rotation(:,:,2)=reshape([0d0,-1d0,1d0,0d0],[2,2])
      rotation(:,:,3)=-rotation(:,:,1); rotation(:,:,4)=-rotation(:,:,2)
      nq=order4(1)
      call gauss(nq,g_h,wg_h)
      call gauss(nv,g_v,wg_v)
      wq=reshape(spread(wg_h,2,nv)*spread(wg_v,1,nq),[npq])

      reverse_nodes=[(nq-mod(j-1,nq)+nq*((j-1)/nq),j=1,npq)]

      do direction=1,4
        use_panel=.false.
        do layer=1,3
          target_groups(:,layer)=region(:,layer) .and. direction_groups(:,direction)
          active(layer)=any(target_groups(:,layer))
          if (active(layer)) use_panel(panel_id(1:panel_count(layer),layer))=.true.
        enddo
        n_panel=count(use_panel)
        if (n_panel==0) cycle
        needed(1:n_panel)=pack([(p,p=1,nlayout)],use_panel)
        panel_row=0; panel_row(needed(1:n_panel))=[(p,p=1,n_panel)]
        parameter_order=maxloc(abs(rotation(:,:,direction)),dim=1)
        C_Xp=C_X; C_Xip=C_Xi
        density_column=reshape([(j,j=1,n)],[16,16])
        if (parameter_order(1)==2) then
          do k=1,3
            C_Xp(:,:,k)=transpose(C_X(:,:,k))
            do j=1,2
              C_Xip(:,:,k,j)=transpose(C_Xi(:,:,k,j))
            enddo
          enddo
          density_column=transpose(density_column)
        endif
        allocate(sq_x(npq*n_panel,3),sq_nx(npq*n_panel,3),sq_w(npq*n_panel))
        allocate(Lu(nq,16,n_panel),Lv(nv,16,n_panel),Pu(nq,16),Pv(nv,16),X(nq,nv),Xtmp(nq,16))
        allocate(rt(nq,nv,3),rp(nq,nv,3))
        do p=1,n_panel
          b=panels_all(:,needed(p)); c=(b([1,3])+b([2,4]))/2; half=(b([2,4])-b([1,3]))/2
          cols=(p-1)*npq
          node_order=[(j,j=1,npq)]
          if (b(1)+b(2)>0d0) node_order=reverse_nodes
          tu=rotation(parameter_order(1),1,direction)*(c(1)+half(1)*g_h)
          tv=rotation(parameter_order(2),2,direction)*(c(2)+half(2)*g_v)
          call legendrevec(15_8,(tu-geometry_center(parameter_order(1))) &
            /geometry_scale(parameter_order(1)),Pu)
          call legendrevec(15_8,(tv-geometry_center(parameter_order(2))) &
            /geometry_scale(parameter_order(2)),Pv)
          do k=1,3
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xp(:,:,k),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,X,nq)
            sq_x(cols+node_order,k)=reshape(X,[npq])
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,1),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rt(:,:,k),nq)
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,2),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rp(:,:,k),nq)
          enddo
          normal_h(:,1)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[npq])
          normal_h(:,2)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[npq])
          normal_h(:,3)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[npq])
          weight_h=sqrt(sum(normal_h**2,dim=2))
          do k=1,3
            sq_nx(cols+node_order,k)=normal_h(:,k)/weight_h
          enddo
          sq_w(cols+node_order)=product(half)*weight_h*wq
          ! Only the canonical-left panel supplies the paired density operators.
          if (b(1)+b(2)<0d0) then
            call legendrevec(15_8,tu,Pu)
            call legendrevec(15_8,tv,Pv)
            call dgemm('N','T',nq,16_8,16_8,1d0,Pu,nq,density_map,16_8,0d0,Lu(:,:,p),nq)
            Xtmp=Lu(:,:,p)
            Lu(:,1:8,p)=(Xtmp(:,1:8)+Xtmp(:,16:9:-1))/2d0
            Lu(:,9:16,p)=(Xtmp(:,1:8)-Xtmp(:,16:9:-1))/2d0
            call dgemm('N','T',nv,16_8,16_8,1d0,Pv,nv,density_map,16_8,0d0,Lv(:,:,p),nv)
          endif
        enddo
        deallocate(Pu,Pv,X,Xtmp,rt,rp)

        do layer=1,3
          if (.not.active(layer)) cycle
          mr=count(target_groups(:,layer)); rows=pack([(i,i=1,m)],target_groups(:,layer))
          target=transpose(tx(:,rows)); np=panel_count(layer)
          ! Column i+mr*(k-1)+n_row*mr*(iv-1): target i, real part of kernel k (imaginary part: k -> n_kernel+k),
          ! short node iv.
          allocate(W_mid(n_row*mr,n),Y(8,n_row*mr*nv))
          allocate(parity_K(nq,n_row*mr*nv,2),parity_W(8,n_row*mr*16,2))
          allocate(v_node(2*npq),dlp_node(2*npq),d_node(2*npq,3),rinv_node(2*npq),r_node(2*npq),K_val(2*npq,n_row))
          parity_W=0d0
          do p=1,np/2
            ip=panel_row(panel_id(panel_pair(1,p,layer),layer))
            ip_right=panel_row(panel_id(panel_pair(2,p,layer),layer))
            do i=1,mr
              ! The mirrored pair (entries ln and npq+ln, ln = iu+nq*(iv-1)), its kernel values, then the pair sums.
              do q_node=1,2
                off_node=(merge(ip,ip_right,q_node==1)-1)*npq; base_node=npq*(q_node-1)
                !GCC$ ivdep
                do ln=1,npq
                  dx=sq_x(off_node+ln,1)-target(i,1); dy=sq_x(off_node+ln,2)-target(i,2); dz=sq_x(off_node+ln,3)-target(i,3)
                  ie=base_node+ln
                  d_node(ie,1)=dx; d_node(ie,2)=dy; d_node(ie,3)=dz
                  distance2_inv=1d0/(dx*dx+dy*dy+dz*dz); rinv_node(ie)=sqrt(distance2_inv)
                  r_node(ie)=(dx*dx+dy*dy+dz*dz)*rinv_node(ie)
                  v_node(ie)=(sq_w(off_node+ln)*inv2pi)*rinv_node(ie)
                  dlp_node(ie)=v_node(ie)*(dx*sq_nx(off_node+ln,1)+dy*sq_nx(off_node+ln,2)+dz*sq_nx(off_node+ln,3))*distance2_inv
                enddo
              enddo
              call compute_helm_kernels(2*npq,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
              do kk=1,n_row
                do iv=1,nv
                  cols=i+mr*(kk-1)+n_row*mr*(iv-1)
                  !GCC$ ivdep
                  do iu=1,nq
                    ln=iu+nq*(iv-1)
                    parity_K(iu,cols,1)=K_val(ln,kk)+K_val(npq+ln,kk)
                    parity_K(iu,cols,2)=K_val(ln,kk)-K_val(npq+ln,kk)
                  enddo
                enddo
              enddo
            enddo
            ! Reduce the paired coordinate to eight density nodes first.
            do mode=1,2
              call dgemm('T','N',8_8,n_row*mr*nv,nq,1d0,Lu(1,1+8*(mode-1),ip),nq, &
                parity_K(1,1,mode),nq,0d0,Y,8_8)
              call dgemm('N','N',8*n_row*mr,16_8,nv,1d0,Y,8*n_row*mr,Lv(:,:,ip),nv, &
                1d0,parity_W(1,1,mode),8*n_row*mr)
            enddo
          enddo
          do iv=1,16
            do iu=1,8
              W_mid(:,density_column(iu,iv))=parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,1) &
                +parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,2)
              W_mid(:,density_column(17-iu,iv))=parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,1) &
                -parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,2)
            enddo
          enddo
          ! Kernel k: real part in rows mr*(k-1)+1:mr*k, imaginary part in rows mr*(n_kernel+k-1)+1:mr*(n_kernel+k).
          do k=1,n_kernel
            row0=mr*(k-1)
            W_c=cmplx(W_mid(row0+1:row0+mr,:),W_mid(mr*n_kernel+row0+1:mr*n_kernel+row0+mr,:),kind=8)
            if (k==1) then
              WS(rows,:)=W_c
            else if (k==2) then
              WD(rows,:)=W_c
            else if (k<=k_Dn) then
              WSn(rows,:,k-k_Sn)=W_c
            else
              WDn(rows,:,i_Dn(k-k_Dn),j_Dn(k-k_Dn))=W_c; WDn(rows,:,j_Dn(k-k_Dn),i_Dn(k-k_Dn))=W_c
            endif
          enddo
          deallocate(W_c,v_node,dlp_node,d_node,rinv_node,r_node,K_val)
          deallocate(W_mid,Y,parity_K,parity_W,target,rows)
        enddo
        deallocate(sq_x,sq_nx,sq_w,Lu,Lv)
      enddo
    end block
  endif

  ! Region 7: tensor sinh for targets with finite source-chart coordinates.
  remaining=remaining .or. off_column4
  where (remaining) classification=7_8
  use_sinh=remaining .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local)
  if (any(use_sinh)) then
    block
      integer(8), parameter :: nlong=192,nshort=32
      integer(8) :: axis(2),n_axis(2),n_panels,it,itarg,np,first,last,l,j_sigma
      real(8) :: c(2),scale(2),J_center(3,2),a2,b2,shear,area,t(2),dc(2),center(2),beta(2)
      real(8) :: C_pos(16,16,3),Xi(16,16,3,2),Xij(16,16,3,3),orientation
      real(8) :: c_s(2),scale_s(2),Xi_s(16,16,3,2),Xij_s(16,16,3,3)
      real(8) :: c_b(2),scale_b(2),Xi_b(16,16,3,2),Xij_b(16,16,3,3)
      logical :: use_typeb,inside
      real(8) :: g32(32),wg32(32),g96(96),wg96(96),g192(192),wg192(192),sigma(16),w_sigma(16)
      real(8) :: edges(3),g_axis(nlong),wg_axis(nlong),a,b,eta0,span,step,eta,z,A0,r0,p0_axis,r0_axis
      real(8) :: alpha,ell,sgn,f,q(nlong,2),d(nlong,2),w(nlong,2),P(nlong,16,2)
      real(8) :: line(nlong*16),B_line(nlong*16,16),T1(nlong,16,2),T0L(nlong,16)
      real(8) :: P0(2,16),J0(3,2),n0(3),offset(3)
      real(8) :: dL(nlong,nshort),dS(nlong,nshort),HSS(nshort),HLS(nlong,nshort),HLL(nlong,nshort)
      real(8) :: XL(nlong,nshort,3),XS(nlong,nshort,3),chord(nlong,nshort,3),R(nlong,nshort,3)
      real(8) :: Nvec(nlong,nshort,3),distance2(nlong,nshort),rinv(nlong,nshort),omega(nlong,nshort)
      real(8) :: kernel_S(nlong,nshort),kernel_D(nlong,nshort),LL(nlong,16),LS(nshort,16)
      real(8) :: Y_leg(16,nshort),W_part7(16,16,2)
      ! K7(:,k): real part of Helmholtz kernel k at the nodes, K7(:,n_kernel+k) its imaginary part
      real(8), allocatable :: K7(:,:),r_dist(:,:)
      complex(8) :: W_c7(256)
      logical :: naive
      ! Restructured evaluation (same nodes and formulas as the reference branches below).
      integer(8) :: il,js
      real(8) :: z_line(nlong),P_prev(nlong),P_cur(nlong),P_next(nlong),T0_line(nlong,16),w_T0,w_T1,rec1,rec2
      real(8) :: C_XLS(16,6,16),C_HLS(16,3,16),C_HLL(16,3,16)
      real(8) :: S_XLS(16,6,nshort),S_HLS(16,3,nshort),S_HLL(16,3,nshort)
      real(8) :: XLS(nlong,6,nshort),HLS3(nlong,3,nshort),HLL3(nlong,3,nshort),HSS3(nshort,3)
      real(8) :: d_L,d_S,hs1,hs2,hs3,r1,r2,r3,c1,c2,c3,Nv1,Nv2,Nv3,dist2,r_inv,om
      interface
        pure function expm1(x) bind(C,name='expm1') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function expm1
        pure function log1p(x) bind(C,name='log1p') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function log1p
      end interface

      axis=[1_8,2_8]
      if (dot_product(rp0,rp0)>dot_product(rt0,rt0)) axis=[2_8,1_8]
      orientation=real(3-2*axis(1),8)
      c=geometry_center(axis); scale=geometry_scale(axis)
      J_center(:,1)=rt0; J_center(:,2)=rp0; J_center=J_center(:,axis)
      a2=dot_product(J_center(:,1),J_center(:,1)); b2=dot_product(J_center(:,2),J_center(:,2))
      shear=dot_product(J_center(:,1),J_center(:,2))
      area=norm2([rt0(2)*rp0(3)-rt0(3)*rp0(2),rt0(3)*rp0(1)-rt0(1)*rp0(3),rt0(1)*rp0(2)-rt0(2)*rp0(1)])
      C_pos=C_X
      Xi=C_Xi(:,:,:,axis); Xij=C_Xij(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
      if (axis(1)==2) then
        do k=1,3
          C_pos(:,:,k)=transpose(C_pos(:,:,k))
          do j=1,2
            Xi(:,:,k,j)=transpose(Xi(:,:,k,j))
          enddo
          do j=1,3
            Xij(:,:,k,j)=transpose(Xij(:,:,k,j))
          enddo
        enddo
      endif
      c_s=c; scale_s=scale; Xi_s=Xi; Xij_s=Xij ! source chart
      use_typeb=all(geometry_scale_typeb>0d0)
      if (use_typeb) then ! type-b chart: the source sits inside it; used when t0 falls inside
        c_b=geometry_center_typeb(axis); scale_b=geometry_scale_typeb(axis)
        Xi_b=C_Xi_typeb(:,:,:,axis); Xij_b=C_Xij_typeb(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
        if (axis(1)==2) then
          do k=1,3
            do j=1,2
              Xi_b(:,:,k,j)=transpose(Xi_b(:,:,k,j))
            enddo
            do j=1,3
              Xij_b(:,:,k,j)=transpose(Xij_b(:,:,k,j))
            enddo
          enddo
        endif
      endif
      call gauss(32_8,g32,wg32)
      call gauss(96_8,g96,wg96)
      call gauss(192_8,g192,wg192)
      sigma=(g16+1d0)/2; w_sigma=wg16/2
      rows=pack([(j,j=1,m)],use_sinh)
      n_axis=[nlong,nshort]
      allocate(K7(nlong*nshort,n_row),r_dist(nlong,nshort))

      do it=1,size(rows,kind=8)
        itarg=rows(it)
        t=t0_local(axis,itarg)
        ! naive=abs(t(2))>2.5d0
        ! naive=off_column4(itarg) .or. abs(t(2))>1.5d0 ! heuristic
        inside=.false.
        if (use_typeb) inside=all(abs(t-c_b)<=scale_b*(1d0+1d-12)) ! t0 inside the type-b chart: geometry is valid at t0
        if (inside .and. .not.off_column4(itarg) .and. root_ok(itarg)) then
          c=c_b; scale=scale_b; Xi=Xi_b; Xij=Xij_b; naive=.false.
        else
          c=c_s; scale=scale_s; Xi=Xi_s; Xij=Xij_s; naive=.not.root_ok(itarg) .or. off_column4(itarg) .or. abs(t(2))>1.5d0
        endif
        dc=min(1d0,max(-1d0,t))-t
        center=t-[shear/a2*dc(2),shear/b2*dc(1)]
        beta=[hypot(h_local(itarg),area/sqrt(a2)*dc(2))/sqrt(a2),hypot(h_local(itarg),area/sqrt(b2)*dc(1))/sqrt(b2)]
        if (.not.root_ok(itarg) .and. (any(.not.ieee_is_finite(center)) .or. any(.not.ieee_is_finite(beta)) .or. &
            any(beta==0d0 .and. abs(center)<=1d0))) then
          use_sinh(itarg)=.false.
          cycle
        endif

        do k=1,2
          n_panels=1; edges(1:2)=[-1d0,1d0]
          if (k==2) then
            np=nshort; g_axis(:np)=g32; wg_axis(:np)=wg32
          else if (abs(center(k))<1d0) then
            np=nlong/2; g_axis(:np)=g96; wg_axis(:np)=wg96
            n_panels=2; edges=[-1d0,center(k),1d0]
          else
            np=nlong; g_axis=g192; wg_axis=wg192
          endif
          do j=1,n_panels
            a=edges(j); b=edges(j+1)
            first=(j-1)*np+1; last=j*np
            if (beta(k)>0d0) then
              eta0=asinh((a-center(k))/beta(k))
              span=asinh((b-center(k))/beta(k))-eta0
              if (abs(center(k))>=2d0) then
                span=2*atanh((b-a)/(hypot(a-center(k),beta(k))+hypot(b-center(k),beta(k))))
              endif
              A0=a-center(k); r0=hypot(A0,beta(k))
              if (A0>=0d0) then
                p0_axis=r0+A0; r0_axis=beta(k)*(beta(k)/p0_axis)
              else
                r0_axis=r0-A0; p0_axis=beta(k)*(beta(k)/r0_axis)
              endif
              do i=first,last
                step=(g_axis(i-first+1)+1d0)*(span/2)
                eta=eta0+step; z=beta(k)*sinh(eta)
                w(i,k)=wg_axis(i-first+1)*((span/2)*beta(k)*cosh(eta))
                if (abs(center(k))<=1d0) then
                  q(i,k)=center(k)+z
                else
                  q(i,k)=a+(p0_axis*expm1(step)-r0_axis*expm1(-step))/2
                endif
                d(i,k)=(center(k)-t(k))+z
              enddo
            else
              alpha=abs(center(k)); sgn=sign(1d0,center(k))
              ell=-log1p(2/(alpha-1d0))
              do i=first,last
                f=(sgn*g_axis(i-first+1)+1d0)/2
                z=-sgn*(alpha+1d0)*exp(f*ell)
                q(i,k)=sgn*(-1d0-(alpha+1d0)*expm1(f*ell))
                w(i,k)=wg_axis(i-first+1)*(-abs(z)*ell/2)
                d(i,k)=(center(k)-t(k))+z
              enddo
            endif
          enddo
          np=n_axis(k)
          call legendrevec(15_8,(q(:np,k)-c(k))/scale(k),P(:np,:,k))
          if (.not.naive) then
            ! Change .false. to .true. (here and below) to force the original evaluation for comparison.
            if (.false.) then
            do j_sigma=1,16
              line((j_sigma-1)*np+1:j_sigma*np)=(t(k)+d(:np,k)*sigma(j_sigma)-c(k))/scale(k)
            enddo
            call legendrevec(15_8,line(:np*16),B_line(:np*16,:))
            do l=1,16
              T1(:np,l,k)=matmul(reshape(B_line(:np*16,l),[np,16_8]),w_sigma*sigma)
              if (k==1) T0L(:,l)=matmul(reshape(B_line(:,l),[nlong,16_8]),w_sigma)
            enddo
            else
              ! Same GL16 sums over sigma, accumulated inside the Legendre recurrence:
              ! T1(i,l) = sum_j w_sigma(j)*sigma(j)*P_{l-1}(z_ij), T0(i,l) = sum_j w_sigma(j)*P_{l-1}(z_ij).
              T1(:np,:,k)=0d0; T0_line(:np,:)=0d0
              do j_sigma=1,16
                w_T1=w_sigma(j_sigma)*sigma(j_sigma); w_T0=w_sigma(j_sigma)
                do i=1,np
                  z_line(i)=(t(k)+d(i,k)*sigma(j_sigma)-c(k))/scale(k)
                  P_prev(i)=1d0; P_cur(i)=z_line(i)
                  T1(i,1,k)=T1(i,1,k)+w_T1; T1(i,2,k)=T1(i,2,k)+P_cur(i)*w_T1
                  T0_line(i,1)=T0_line(i,1)+w_T0; T0_line(i,2)=T0_line(i,2)+P_cur(i)*w_T0
                enddo
                do l=1,14
                  rec1=real(2*l+1,8)/real(l+1,8); rec2=real(l,8)/real(l+1,8)
                  !GCC$ ivdep
                  do i=1,np
                    P_next(i)=rec1*z_line(i)*P_cur(i)-rec2*P_prev(i)
                    T1(i,l+2,k)=T1(i,l+2,k)+P_next(i)*w_T1
                    T0_line(i,l+2)=T0_line(i,l+2)+P_next(i)*w_T0
                    P_prev(i)=P_cur(i); P_cur(i)=P_next(i)
                  enddo
                enddo
              enddo
              if (k==1) T0L=T0_line
            endif
          endif
        enddo

        ! Naive targets keep the original evaluation; .false. -> .true. forces it for every target.
        if (.false. .or. naive) then
        do k=1,3
          XL(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,1)),transpose(P(:nshort,:,2)))
          XS(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,2)),transpose(P(:nshort,:,2)))
        enddo
        if (naive) then
          do k=1,3
            chord(:,:,k)=matmul(matmul(P(:,:,1),C_pos(:,:,k)),transpose(P(:nshort,:,2)))-tx(k,itarg)
          enddo
          R=chord
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          dL=spread(d(:,1),2,nshort); dS=spread(d(:nshort,2),1,nlong)
          do k=1,3
            HSS=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
            HLS=matmul(matmul(T0L,Xij(:,:,k,2)),transpose(P(:nshort,:,2)))
            HLL=matmul(matmul(T1(:,:,1),Xij(:,:,k,1)),transpose(P(:nshort,:,2)))
            R(:,:,k)=offset(k)-dS**2*spread(HSS,1,nlong)-dL*dS*HLS-dL**2*HLL
            chord(:,:,k)=dL*XL(:,:,k)+dS*XS(:,:,k)+R(:,:,k)
          enddo
        endif

        Nvec(:,:,1)=orientation*(XS(:,:,2)*XL(:,:,3)-XS(:,:,3)*XL(:,:,2))
        Nvec(:,:,2)=orientation*(XS(:,:,3)*XL(:,:,1)-XS(:,:,1)*XL(:,:,3))
        Nvec(:,:,3)=orientation*(XS(:,:,1)*XL(:,:,2)-XS(:,:,2)*XL(:,:,1))
        distance2=sum(chord**2,dim=3)
        rinv=1d0/sqrt(distance2)
        omega=inv2pi*spread(w(:,1),2,nshort)*spread(w(:nshort,2),1,nlong)
        kernel_S=omega*sqrt(sum(Nvec**2,dim=3))*rinv
        kernel_D=omega*sum(Nvec*R,dim=3)*rinv/distance2
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          ! Coefficients stacked as (long mode, field, short mode): XL (fields 1:3) and XS (4:6); HLS; HLL.
          do l=1,16
            do k=1,3
              C_XLS(:,k,l)=Xi(:,l,k,1); C_XLS(:,3+k,l)=Xi(:,l,k,2)
              C_HLS(:,k,l)=Xij(:,l,k,2); C_HLL(:,k,l)=Xij(:,l,k,1)
            enddo
          enddo
          ! X = P_long*(C*P_short^T): the short direction first, then one GEMM per group along the long direction.
          call dgemm('N','T',96_8,nshort,16_8,1d0,C_XLS,96_8,P(1,1,2),nlong,0d0,S_XLS,96_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLS,48_8,P(1,1,2),nlong,0d0,S_HLS,48_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLL,48_8,P(1,1,2),nlong,0d0,S_HLL,48_8)
          call dgemm('N','N',nlong,6*nshort,16_8,1d0,P(1,1,1),nlong,S_XLS,16_8,0d0,XLS,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T0L,nlong,S_HLS,16_8,0d0,HLS3,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T1(1,1,1),nlong,S_HLL,16_8,0d0,HLL3,nlong)
          do k=1,3
            HSS3(:,k)=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
          enddo
          ! R, chord, normal and both kernels at each node, same expressions as the reference.
          do js=1,nshort
            d_S=d(js,2)
            hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
            !GCC$ ivdep
            do il=1,nlong
              d_L=d(il,1)
              r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
              r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
              r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
              c1=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
              c2=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
              c3=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
              Nv1=orientation*(XLS(il,5,js)*XLS(il,3,js)-XLS(il,6,js)*XLS(il,2,js))
              Nv2=orientation*(XLS(il,6,js)*XLS(il,1,js)-XLS(il,4,js)*XLS(il,3,js))
              Nv3=orientation*(XLS(il,4,js)*XLS(il,2,js)-XLS(il,5,js)*XLS(il,1,js))
              dist2=c1**2+c2**2+c3**2
              r_inv=1d0/sqrt(dist2)
              om=inv2pi*w(il,1)*w(js,2)
              kernel_S(il,js)=om*sqrt(Nv1**2+Nv2**2+Nv3**2)*r_inv
              kernel_D(il,js)=om*(Nv1*r1+Nv2*r2+Nv3*r3)*r_inv/dist2
            enddo
          enddo
          ! chord and distance2 at every node: the Helmholtz factors need the distance, SLPn and DLPn the chord.
            do js=1,nshort
              d_S=d(js,2)
              hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
              do il=1,nlong
                d_L=d(il,1)
                r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
                r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
                r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
                chord(il,js,1)=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
                chord(il,js,2)=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
                chord(il,js,3)=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
                distance2(il,js)=chord(il,js,1)**2+chord(il,js,2)**2+chord(il,js,3)**2
              enddo
            enddo
        endif

        ! Helmholtz kernels at the nodes from the Laplace SLP, DLP values and the chord.
        rinv=1d0/sqrt(distance2); r_dist=distance2*rinv
        call compute_helm_kernels(nlong*nshort,kernel_S,kernel_D,chord,rinv,r_dist,K7)

        call legendrevec(15_8,q(:,1),P(:,:,1))
        call legendrevec(15_8,q(:nshort,2),P(:nshort,:,2))
        LL=matmul(P(:,:,1),transpose(density_map)); LS=matmul(P(:nshort,:,2),transpose(density_map))
        do k=1,n_kernel   ! real and imaginary parts contracted separately
          do part=1,2
            call dgemm('T','N',16_8,nshort,nlong,1d0,LL,nlong,K7(1,k+n_kernel*(part-1)),nlong,0d0,Y_leg,16_8)
            call dgemm('N','N',16_8,16_8,nshort,1d0,Y_leg,16_8,LS,nshort,0d0,W_part7(1,1,part),16_8)
            if (axis(1)==2) W_part7(:,:,part)=transpose(W_part7(:,:,part))
          enddo
          W_c7=reshape(cmplx(W_part7(:,:,1),W_part7(:,:,2),kind=8),[256])
          if (k==1) then
            WS(itarg,:)=W_c7
          else if (k==2) then
            WD(itarg,:)=W_c7
          else if (k<=k_Dn) then
            WSn(itarg,:,k-k_Sn)=W_c7
          else
            WDn(itarg,:,i_Dn(k-k_Dn),j_Dn(k-k_Dn))=W_c7; WDn(itarg,:,j_Dn(k-k_Dn),i_Dn(k-k_Dn))=W_c7
          endif
        enddo
      enddo
      deallocate(rows)
      deallocate(K7,r_dist)
    end block
  endif
  remaining=remaining .and. .not.use_sinh
  where (remaining) ier=1_8


contains

  ! Helmholtz kernels at nb nodes from the Laplace SLP value v_in (weight/(2*pi*r)), DLP value dlp_in, the chord
  ! d_in = y-x, 1/r and r: K_out(:,k) is the real part of kernel k, K_out(:,n_kernel+k) its imaginary part.
  ! exp(1i*zk*r) = exp(-Im(zk)*r)*(cos+1i*sin)(Re(zk)*r), evaluated on the whole batch.
  subroutine compute_helm_kernels(nb,v_in,dlp_in,d_in,rinv_in,r_in,K_out)
    integer(8), intent(in) :: nb
    real(8), intent(in) :: v_in(nb),dlp_in(nb),d_in(nb,3),rinv_in(nb),r_in(nb)
    real(8), intent(out) :: K_out(nb,n_row)
    real(8), allocatable :: exp_decay(:),e_re(:),e_im(:),c_re(:),c_im(:)
    real(8), allocatable :: f_in(:),P_in(:),dd(:)
    real(8) :: kr,ki,k2r,k2i,delta
    integer(8) :: kn_h,a,b
    kr=real(zk,8); ki=aimag(zk); k2r=kr*kr-ki*ki; k2i=2d0*kr*ki
    allocate(exp_decay(nb),e_re(nb),e_im(nb),c_re(nb),c_im(nb))
    ! e = exp(1i*z), c = e*(1-1i*z) = e*((1+Im(zk)*r)-1i*Re(zk)*r)
    ! Real wavenumbers have no exponential decay; avoid a scalar exp(0)
    ! evaluation at every quadrature node in this hot path.
    if (ki==0d0) then
      exp_decay=1d0
    else
#ifdef SURFQUADR_USE_VFORCE
      c_re=-ki*r_in
      call vvexp(exp_decay,c_re,int(nb,c_int))
#else
      exp_decay=exp(-ki*r_in)
#endif
    endif
#ifdef SURFQUADR_USE_VFORCE
    ! Reuse the factor buffer for the phase; vForce expects a C int by reference.
    c_re=kr*r_in
    call vvsincos(e_im,e_re,c_re,int(nb,c_int))
    e_re=exp_decay*e_re; e_im=exp_decay*e_im
#else
    e_re=exp_decay*cos(kr*r_in); e_im=exp_decay*sin(kr*r_in)
#endif
    c_re=e_re*(1d0+ki*r_in)+e_im*(kr*r_in); c_im=e_im*(1d0+ki*r_in)-e_re*(kr*r_in)
    K_out(:,1)=v_in*e_re; K_out(:,n_kernel+1)=v_in*e_im
    K_out(:,2)=dlp_in*c_re; K_out(:,n_kernel+2)=dlp_in*c_im
    if (n_kernel>2_8) then
      allocate(f_in(nb),P_in(nb),dd(nb))
      f_in=v_in*rinv_in*rinv_in
      do kn_h=1,3*ifSn
        K_out(:,k_Sn+kn_h)=f_in*d_in(:,kn_h)*c_re; K_out(:,n_kernel+k_Sn+kn_h)=f_in*d_in(:,kn_h)*c_im
      enddo
      do kn_h=1,6*ifDn
        a=i_Dn(kn_h); b=j_Dn(kn_h); delta=merge(1d0,0d0,a==b)
        dd=d_in(:,a)*d_in(:,b); P_in=3d0*dd*rinv_in*rinv_in-delta
        K_out(:,k_Dn+kn_h)=f_in*(P_in*c_re-dd*(k2r*e_re-k2i*e_im))
        K_out(:,n_kernel+k_Dn+kn_h)=f_in*(P_in*c_im-dd*(k2r*e_im+k2i*e_re))
      enddo
    endif
  end subroutine compute_helm_kernels

end subroutine quadSD_helm_mid_w16


! Stokes version of quadSD_mid_w16 (viscosity 1, the kernels of quadSD_stok_w16): SLP, DLP and SLPn (the traction of the SLP at the
! target normal tnx) evaluated at the integration nodes of every region, then contracted as in quadSD_mid_w16.
! At a node with d = y-x, r = |d|, Laplace SLP value v = w/(2*pi*r) and DLP value dlp = v*(d.n_s)/r^2 (n_s the normal of
! WD, snx), dd = d_a*d_b/r^2, dtn = d.tnx:
!   SLP_ab = (1/4)*v*(delta_ab+dd),  DLP_ab = (3/2)*dlp*dd,  SLPn_ab = (3/2)*v*dtn*dd/r^2   (compute_stok_kernels),
! i.e. (I/r+r*r'/r^3)/(8*pi), -3*((x-y).n_s)*r*r'/(4*pi*r^5), -3*((x-y).tnx)*r*r'/(4*pi*r^5) times the weight.
! Six components per kernel, (a,b) = (1,1),(1,2),(1,3),(2,2),(2,3),(3,3), stacked as kernel rows k_S+c, k_D+c, k_Sn+c;
! the outputs WS, WD, WSn are (m,n,3,3), symmetric. Region 7 uses its stable kernel_D and chord.
! DLPn (the traction of the DLP at tnx, the negative of Sto3dDLPmat's second output with source normal n_s, 1/(4*pi),
! as quadDn_stok_{on,off}_w16): with tns = tnx.n_s, dns = d.n_s (v*dns = dlp*r^2),
!   DLPn_ij = -[v/r^2*(tnx_i*n_s_j+(3/2)*tns*dd_ij+(3/2)*dtn*n_s_i*d_j/r^2)
!               +dlp/r^2*((3/2)*d_i*tnx_j+dtn*((3/2)*delta_ij-15*dd_ij))],
! nine components (i,j) at i+3*(j-1), kernel rows k_Dn+c; WDn is (m,n,3,3), not symmetric.
subroutine quadSD_stok_mid_w16(m,tx,tnx,n,sx,snx,sw,WS,WD,WSn,WDn,ifS,ifD,ifSn,ifDn, &
    t0,h,C_X,C_Xi,C_Xij,geometry_center,geometry_scale,classification, &
    C_Xi_typeb,C_Xij_typeb,geometry_center_typeb,geometry_scale_typeb,ier,root_ok)
  use surfquadr_mid, only: legendrevec, gauss => compute_gaussleg
  use, intrinsic :: iso_c_binding, only: c_double,c_bool
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  implicit none
  integer(8), intent(in) :: m,n
  real(8), intent(in) :: tx(3,m),tnx(3,m),sx(3,n),snx(3,n),sw(n),t0(2,m),h(m)
  real(8), intent(in) :: C_X(16,16,3),C_Xi(16,16,3,2),C_Xij(16,16,3,3)
  real(8), intent(in) :: geometry_center(2),geometry_scale(2)
  ! *_typeb: geometry of the type-b chart the source sits in; geometry_scale_typeb = 0: no type-b chart.
  real(8), intent(in) :: C_Xi_typeb(16,16,3,2),C_Xij_typeb(16,16,3,3)
  real(8), intent(in) :: geometry_center_typeb(2),geometry_scale_typeb(2)
  integer(8), intent(in) :: ifS,ifD,ifSn,ifDn
  real(8), intent(inout) :: WS(m,n*ifS,3,3),WD(m,n*ifD,3,3),WSn(m,n*ifSn,3,3),WDn(m,n*ifDn,3,3)
  integer(8), intent(inout) :: classification(m),ier(m)
  logical(c_bool), intent(in) :: root_ok(m)
  real(8) :: t0_local(2,m),h_local(m)
  ! The 256 source nodes are ordered with u varying fastest.
  ! The caller excludes near/self and supplies the Region 7 projection and signed height.
  ! Regions 1:3 use 1D parity factors for every kernel.
  ! Kernel rows k_S+c (SLP), k_D+c (DLP), k_Sn+c (SLPn), component c = (a_comp(c),b_comp(c)); k_Dn+c (DLPn), c = 1..9.
  integer(8), parameter :: a_comp(6)=[1_8,1_8,1_8,2_8,2_8,3_8], b_comp(6)=[1_8,2_8,3_8,2_8,3_8,3_8]
  integer(8) :: n_kernel,k_S,k_D,k_Sn,k_Dn,kn,q_node
  integer(8) :: i,j,k,m16,mr,nq,nh,ld,layer,info,piv(16)
  ! integer(8), parameter :: order_p(3)=[32,64,128]
  integer(8), parameter :: order_p(3)=[32,48,96]
  integer(8), allocatable :: idx16(:),rows(:)
  real(8) :: g16(16),wg16(16),V16(16,16),Vsolve(16,16),density_map(16,16),source_t(2,256)
  real(8) :: xc(3),rt0(3),rp0(3),e1(3),e2(3),e3(3),frame(3,3),r(3),xyz(3,m)
  real(8) :: ratio
  real(8) :: lambda1,d16,d32,d64,d128,abc32(3),abc64(3),abc128(3)
  real(8), parameter :: pi=acos(-1.0d0),inv2pi=0.5d0/pi
  logical(c_bool) :: outside16(m),inside16(m)
  logical :: inside32(m),inside64(m),inside128(m),shell32(m),shell64(m),shell128(m)
  logical :: target_groups(m,3),remaining(m),use_sinh(m),off_column4(m)
  real(8), allocatable :: Pu(:,:),Pv(:,:),X(:,:),rt(:,:,:),rp(:,:,:),Xtmp(:,:)
  real(8), allocatable :: g(:),wg(:),L(:,:,:),sq_x(:,:),sq_nx(:,:),sq_w(:),target(:,:),Y(:,:)
  real(8) :: P_half(maxval(order_p)/2,16),L_half(maxval(order_p)/2,16)
  real(8) :: dx,dy,dz,distance2_inv,distance_inv,kernel_value
  integer(8) :: iu,iv,mode,iq,is,node(4,4096),column(4,64)
  real(8) :: v1,v2,v3,v4,p1,p2,p3,p4,d1,d2,d3,d4
  real(8), allocatable :: parity_K(:,:,:),parity_W(:,:,:)
  external :: dgemm,dgesv
  ! Laplace values, chord, 1/r, chord.tnx and the Stokes kernels at the nodes of one evaluation step
  ! (compute_stok_kernels).
  integer(8) :: n_row,row0,kk,jn,ie,ln,off_node,base_node
  real(8) :: w_node
  real(8), allocatable :: v_node(:),dlp_node(:),d_node(:,:),rinv_node(:),dtn_node(:),K_val(:,:),W_part(:,:)
  ! DLPn only: the source normal n_s and the target normal tnx at the same entries
  real(8), allocatable :: snx_node(:,:),tnx_node(:,:)
  real(8), allocatable :: target_n(:,:)

  ier=0_8
  n_kernel=6_8*(ifS+ifD+ifSn)+9_8*ifDn; k_S=0_8; k_D=6_8*ifS; k_Sn=6_8*(ifS+ifD); k_Dn=6_8*(ifS+ifD+ifSn)
  n_row=n_kernel
  if (n_kernel==0_8) return
  call prepare_quadSD_targets_w16(m,tx,n,sx,C_X,C_Xi,geometry_center,geometry_scale, &
      inside16,outside16,xc,rt0,rp0,e1,e2,e3,r,xyz,d16,d128,lambda1)
  frame(:,1)=e1; frame(:,2)=e2; frame(:,3)=e3

  d32=0.50d0*d16; d64=0.25d0*d16
  abc32=lambda1*sqrt(r*r+d32*d32)
  abc64=lambda1*sqrt(r*r+d64*d64); abc128=lambda1*sqrt(r*r+d128*d128)

  ! 1st layer of guard.
  where (outside16) classification=0_8
  m16=count(outside16); idx16=pack([(k,k=1,m)],outside16)
  if (m16>0) then
    target=tx(:,idx16)
    allocate(v_node(m16),dlp_node(m16),d_node(m16,3),rinv_node(m16),dtn_node(m16),K_val(m16,n_row))
    allocate(snx_node(m16,3),tnx_node(m16,3))
    if (ifDn/=0_8) tnx_node=transpose(tnx(:,idx16))
    do j=1,n
      if (ifDn/=0_8) snx_node=spread(snx(:,j),1,m16)
      !GCC$ ivdep
      do i=1,m16
        dx=sx(1,j)-target(1,i); dy=sx(2,j)-target(2,i); dz=sx(3,j)-target(3,i)
        d_node(i,1)=dx; d_node(i,2)=dy; d_node(i,3)=dz
        distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz); rinv_node(i)=sqrt(distance2_inv)
        dtn_node(i)=(dx*tnx(1,idx16(i))+dy*tnx(2,idx16(i)))+dz*tnx(3,idx16(i))
        v_node(i)=(sw(j)*inv2pi)*rinv_node(i)
        dlp_node(i)=v_node(i)*(dx*snx(1,j)+dy*snx(2,j)+dz*snx(3,j))*distance2_inv
      enddo
      call compute_stok_kernels(m16,v_node,dlp_node,d_node,rinv_node,dtn_node,snx_node,tnx_node,K_val)
      do k=1,n_kernel
        call store_stok_weights(idx16,j,j,k,K_val(:,k))
      enddo
    enddo
    deallocate(target,v_node,dlp_node,d_node,rinv_node,dtn_node,K_val,snx_node,tnx_node)
  endif
  deallocate(idx16)

  inside32=inside16 .and. sum((xyz/spread(abc32,2,m))**2,dim=1)<=1.0d0
  shell32=inside16 .and. .not.inside32
  inside64=inside32 .and. sum((xyz/spread(abc64,2,m))**2,dim=1)<=1.0d0
  shell64=inside32 .and. .not.inside64
  inside128=inside64 .and. sum((xyz/spread(abc128,2,m))**2,dim=1)<=1.0d0
  shell128=inside64 .and. .not.inside128
  target_groups(:,1)=shell32; target_groups(:,2)=shell64; target_groups(:,3)=shell128

  if (any(inside16)) then
    call gauss(16_8,g16,wg16)
    call legendrevec(15_8,g16,V16)
    Vsolve=transpose(V16); density_map=0.0d0
    do j=1,16
      density_map(j,j)=1.0d0
    enddo
    call dgesv(16_8,16_8,Vsolve,16_8,piv,density_map,16_8,info)
    if (info/=0) return
  endif

  ! Common p refinement: contract v, then u for order_p.
  do layer=1,3
    where (target_groups(:,layer)) classification=layer
    mr=count(target_groups(:,layer))
    if (mr==0) cycle
    rows=pack([(k,k=1,m)],target_groups(:,layer))
    ! nq=16*2**layer; nh=nq/2
    nq=order_p(layer); nh=nq/2
    allocate(g(nq),wg(nq),L(nh,8,2))
    call gauss(nq,g,wg)
    call legendrevec(15_8,g(:nh),P_half(:nh,:))
    call dgemm('N','T',nh,16_8,16_8,1d0,P_half,size(P_half,1,kind=8), &
      density_map,16_8,0d0,L_half,size(L_half,1,kind=8))
    L(:,:,1)=0.5d0*(L_half(:nh,1:8)+L_half(:nh,16:9:-1))
    L(:,:,2)=0.5d0*(L_half(:nh,1:8)-L_half(:nh,16:9:-1))
    allocate(Pu(nq,16),Pv(nq,16),X(nq,nq),rt(nq,nq,3),rp(nq,nq,3),Xtmp(nq,16))
    allocate(sq_x(3,nq*nq),sq_nx(3,nq*nq),sq_w(nq*nq))
    call legendrevec(15_8,(g-geometry_center(1))/geometry_scale(1),Pu)
    call legendrevec(15_8,(g-geometry_center(2))/geometry_scale(2),Pv)
    do k=1,3
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_X(:,:,k),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,X,nq)
      sq_x(k,:)=reshape(X,[nq*nq])
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,1),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rt(:,:,k),nq)
      call dgemm('N','N',nq,16_8,16_8,1.0d0,Pu,nq,C_Xi(:,:,k,2),16_8,0.0d0,Xtmp,nq)
      call dgemm('N','T',nq,nq,16_8,1.0d0,Xtmp,nq,Pv,nq,0.0d0,rp(:,:,k),nq)
    enddo
    sq_nx(1,:)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[nq*nq])
    sq_nx(2,:)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[nq*nq])
    sq_nx(3,:)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[nq*nq])
    sq_w=sqrt(sum(sq_nx*sq_nx,dim=1)); sq_nx=sq_nx/spread(sq_w,1,3)
    sq_w=sq_w*reshape(spread(wg,2,nq)*spread(wg,1,nq),[nq*nq])
    target=transpose(tx(:,rows)); target_n=transpose(tnx(:,rows))
    ! Keep targets contiguous; contract v, then u without transposing kernels.
    ! Row i+mr*(k-1): target i, kernel row k.
    ld=8*((n_row*mr+7)/8)
    allocate(parity_K(ld,nh*nh,4),parity_W(ld,64,4),Y(ld*nh,8))
    allocate(v_node(4*mr),dlp_node(4*mr),d_node(4*mr,3),rinv_node(4*mr),dtn_node(4*mr),K_val(4*mr,n_row))
    allocate(snx_node(4*mr,3),tnx_node(4*mr,3))
    if (ifDn/=0_8) then
      do q_node=1,4
        tnx_node(mr*(q_node-1)+1:mr*q_node,:)=target_n(1:mr,:)
      enddo
    endif
    parity_K(n_row*mr+1:ld,:,:)=0.0d0
    do iv=1,nh
      do iu=1,nh
        iq=iu+nh*(iv-1)
        node(:,iq)=[iu+nq*(iv-1),nq+1-iu+nq*(iv-1),iu+nq*(nq-iv),nq+1-iu+nq*(nq-iv)]
      enddo
    enddo
    do iv=1,8
      do iu=1,8
        is=iu+8*(iv-1)
        column(:,is)=[iu+16*(iv-1),17-iu+16*(iv-1),iu+16*(16-iv),17-iu+16*(16-iv)]
      enddo
    enddo
    do j=1,nh*nh
      ! The four mirrored nodes (entries i+mr*(q_node-1)), their kernel values, then the parity sums.
      do q_node=1,4
        jn=node(q_node,j); w_node=sq_w(jn)*inv2pi
        !GCC$ ivdep
        do i=1,mr
          dx=sq_x(1,jn)-target(i,1); dy=sq_x(2,jn)-target(i,2); dz=sq_x(3,jn)-target(i,3)
          ie=i+mr*(q_node-1)
          d_node(ie,1)=dx; d_node(ie,2)=dy; d_node(ie,3)=dz
          distance2_inv=1.0d0/(dx*dx+dy*dy+dz*dz); rinv_node(ie)=sqrt(distance2_inv)
          dtn_node(ie)=(dx*target_n(i,1)+dy*target_n(i,2))+dz*target_n(i,3)
          v_node(ie)=w_node*rinv_node(ie)
          dlp_node(ie)=v_node(ie)*(dx*sq_nx(1,jn)+dy*sq_nx(2,jn)+dz*sq_nx(3,jn))*distance2_inv
        enddo
        if (ifDn/=0_8) snx_node(mr*(q_node-1)+1:mr*q_node,:)=spread(sq_nx(:,jn),1,mr)
      enddo
      call compute_stok_kernels(4*mr,v_node,dlp_node,d_node,rinv_node,dtn_node,snx_node,tnx_node,K_val)
      do kk=1,n_row
        row0=mr*(kk-1)
        !GCC$ ivdep
        do i=1,mr
          p1=K_val(i,kk)+K_val(mr+i,kk); p2=K_val(i,kk)-K_val(mr+i,kk)
          p3=K_val(2*mr+i,kk)+K_val(3*mr+i,kk); p4=K_val(2*mr+i,kk)-K_val(3*mr+i,kk)
          parity_K(row0+i,j,1)=p1+p3; parity_K(row0+i,j,2)=p2+p4
          parity_K(row0+i,j,3)=p1-p3; parity_K(row0+i,j,4)=p2-p4
        enddo
      enddo
    enddo
    do mode=1,4
      call dgemm('N','N',ld*nh,8_8,nh,1.0d0,parity_K(1,1,mode),ld*nh, &
        L(1,1,1+(mode-1)/2),nh,0.0d0,Y,ld*nh)
      do iv=1,8
        call dgemm('N','N',n_row*mr,8_8,nh,1.0d0,Y(1,iv),ld,L(1,1,1+mod(mode-1,2)),nh, &
          0.0d0,parity_W(1,1+8*(iv-1),mode),ld)
      enddo
    enddo
    allocate(W_part(mr,n))
    do k=1,n_kernel
      row0=mr*(k-1)
      do j=1,64
        !GCC$ ivdep
        do i=1,mr
          v1=parity_W(row0+i,j,1); v2=parity_W(row0+i,j,2)
          v3=parity_W(row0+i,j,3); v4=parity_W(row0+i,j,4)
          p1=v1+v2; p2=v3+v4; p3=v1-v2; p4=v3-v4
          W_part(i,column(1,j))=p1+p2; W_part(i,column(2,j))=p3+p4
          W_part(i,column(3,j))=p1-p2; W_part(i,column(4,j))=p3-p4
        enddo
      enddo
      do j=1,n
        call store_stok_weights(rows,j,j,k,W_part(:,j))
      enddo
    enddo
    deallocate(W_part,v_node,dlp_node,d_node,rinv_node,dtn_node,K_val,snx_node,tnx_node)
    deallocate(parity_K,parity_W,Y)
    deallocate(Pu,Pv,X,Xtmp,rt,rp,sq_x,sq_nx,sq_w,target,target_n)
    deallocate(g,wg,L,rows)
  enddo

  remaining=inside128
  off_column4=.false.
  t0_local=t0; h_local=h
  if (any(inside128)) then
    do iv=1,16
      do iu=1,16
        source_t(:,iu+16*(iv-1))=[g16(iu),g16(iv)]
      enddo
    enddo
  endif

  ! Regions 4:6: panel refinement for patches up to 10:1.
  ratio = norm2(rp0)/norm2(rt0)
  ratio = max(ratio,1d0/ratio)
  ! if (any(remaining) .and. ratio <= 10d0) then
  if (any(remaining) .and. ratio <= 8d0) then ! heuristic
    block
      integer(8), parameter :: order4(2)=[64,16] ! Common orders for regions 4:6.
      integer(8), parameter :: nv=order4(2),npq=product(order4)
      integer(8), parameter :: preset=2 ! 1: five rows; 2: three rows.
      integer(8) :: p,ip,ip_right,pr,direction,n_panel,np,cols,parameter_order(2),density_column(16,16)
      integer(8) :: panel_count(3),row,nparent,nchild,nsplit,nlayout,first_panel,last_panel
      integer(8), allocatable :: ncols(:),panel_id(:,:),panel_row(:),needed(:),panel_mirror(:),panel_pair(:,:,:)
      integer(8) :: node_order(npq),reverse_nodes(npq)
      real(8) :: box_center(3),box_frame(3,3),box_halfsize4(3),box_halfsize5(3),box_halfsize6(3)
      real(8) :: scale_box4,scale_box5,scale_box6,d_box4,d_box5,d_box6
      real(8), allocatable :: row_ratio(:),split_ratio(:),v_edges(:),panels_all(:,:)
      real(8) :: width,top_height
      real(8) :: b(4),c(2),half(2),rotation(2,2,4)
      real(8) :: tu(order4(1)),tv(nv),g_h(order4(1)),wg_h(order4(1)),g_v(nv),wg_v(nv),wq(npq)
      real(8) :: C_Xp(16,16,3),C_Xip(16,16,3,2)
      real(8) :: normal_h(npq,3),weight_h(npq)
      real(8), allocatable :: Lu(:,:,:),Lv(:,:,:),W_mid(:,:)
      logical :: inside_box4(m),inside_box5(m),inside_box6(m),region(m,3),direction_groups(m,4),active(3)
      logical, allocatable :: use_panel(:)
      ! Cylindrical bounds around the self patch for regions 4:6.
      integer(8), parameter :: p_R=3,n_R=(p_R+1)*(p_R+2)/2,n_b=68
      integer(8) :: i_R(n_R),j_R(n_R),iR,jR,kR,itarg
      real(8) :: m_box(3),H_box,gram(3),t_b(n_b,2),P_bu(n_b,16),P_bv(n_b,16),x_b(3,n_b),xyz_s(3,n+n_b),s_R(2)
      real(8) :: B_s(n+n_b,n_R),B_fit(n,n_R),c_R(n,2),r_s(2),zeta_lo,zeta_hi,ell_R(2),b_box(3,2),lapack_work(4096)
      real(8), allocatable :: P_a(:,:),P_b(:,:),B_t(:,:),R_t(:,:)
      real(8) :: P_est(2,16),X_est(3),J_est(3,2),n_est(3),t_direction(2,m)
      logical :: in_slab(m),beside(m)
      external :: dgels

      box_center=xc; box_frame=frame
      t_direction=t0_local
      if (.false.) then ! old rule, turns out to be not valid for sphere, near triple junction
        ! scale_box4=0.95d0; d_box4=0.15d0*maxval(r(1:2))
        scale_box4=0.95d0; d_box4=0.0425d0*maxval(r(1:2))
        ! scale_box5=0.95d0; d_box5=0.080d0*maxval(r(1:2))
        scale_box5=0.65d0; d_box5=0.0d0*maxval(r(1:2))
        ! scale_box6=0.95d0; d_box6=0.0d0*maxval(r(1:2))
        ! scale_box6=0.90d0; d_box6=0.0d0*maxval(r(1:2))
        scale_box6=0.60d0; d_box6=0.0d0*maxval(r(1:2))
        box_halfsize4=scale_box4*r+d_box4
        box_halfsize5=scale_box5*r+d_box5
        box_halfsize6=scale_box6*r+d_box6
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        inside_box4=all(abs(xyz)<=spread(box_halfsize4,2,m),dim=1)
        inside_box5=all(abs(xyz)<=spread(box_halfsize5,2,m),dim=1)
        inside_box6=all(abs(xyz)<=spread(box_halfsize6,2,m),dim=1)
      else
        ! Boxes 4:6 as cylinders around the self patch.
        ! (xi,eta,zeta) = box_frame^T (x-box_center); (u,v) = R(xi,eta): least-squares inverse fit of degree p_R
        ! on the source nodes. box_k = {|R_i| <= r_s,i + m_k*l_other/l_i, zeta_lo-H_box <= zeta <= zeta_hi+H_box};
        ! r_s = extent of R over self, so self is inside every box.
        ! Region 7 = box6; a target with trusted (t0,h) must also satisfy |t0_i| <= b6_i
        ! (keeps targets beside a curved edge out of Region 7).
        !
        ! things to tweak inside, heuristic: p_R (parameter above), m_box, H_box
        m_box=[0.10d0,0.06d0,0.02d0]
        gram=[dot_product(rt0,rt0),dot_product(rt0,rp0),dot_product(rp0,rp0)]
        H_box=0.6d0*sqrt(0.5d0*(gram(1)+gram(3))+sqrt(0.25d0*(gram(1)-gram(3))**2+gram(2)**2)) ! 0.6*||[rt0,rp0]||_2
        !
        kR=0 ! total degree <= p_R
        do jR=1,p_R+1
          do iR=1,p_R+2-jR
            kR=kR+1; i_R(kR)=iR; j_R(kR)=jR
          enddo
        enddo
        t_b(1:16,1)=g16; t_b(1:16,2)=-1d0 ! self boundary: GL16 on each edge + corners
        t_b(17:32,1)=1d0; t_b(17:32,2)=g16
        t_b(33:48,1)=g16; t_b(33:48,2)=1d0
        t_b(49:64,1)=-1d0; t_b(49:64,2)=g16
        t_b(65:68,1)=[-1d0,1d0,1d0,-1d0]; t_b(65:68,2)=[-1d0,-1d0,1d0,1d0]
        call legendrevec(15_8,(t_b(:,1)-geometry_center(1))/geometry_scale(1),P_bu)
        call legendrevec(15_8,(t_b(:,2)-geometry_center(2))/geometry_scale(2),P_bv)
        do k=1,3
          x_b(k,:)=sum(matmul(P_bu,C_X(:,:,k))*P_bv,dim=2)
        enddo
        xyz_s(:,1:n)=matmul(transpose(box_frame),sx-spread(box_center,2,n)) ! self samples: n nodes, then boundary
        xyz_s(:,n+1:)=matmul(transpose(box_frame),x_b-spread(box_center,2,n_b))
        s_R=maxval(abs(xyz_s(1:2,1:n)),dim=2)
        allocate(P_a(n+n_b,p_R+1),P_b(n+n_b,p_R+1))
        call legendrevec(p_R,xyz_s(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz_s(2,:)/s_R(2),P_b)
        B_s=P_a(:,i_R)*P_b(:,j_R)
        deallocate(P_a,P_b)
        B_fit=B_s(1:n,:); c_R=transpose(source_t) ! fit on the n source nodes
        call dgels('N',n,n_R,2_8,B_fit,n,c_R,n,lapack_work,4096_8,info)
        r_s=maxval(abs(matmul(B_s,c_R(1:n_R,:))),dim=1) ! self extent of R (nodes and boundary)
        zeta_lo=minval(xyz_s(3,:)); zeta_hi=maxval(xyz_s(3,:))
        xyz=matmul(transpose(box_frame),tx-spread(box_center,2,m))
        allocate(P_a(m,p_R+1),P_b(m,p_R+1),B_t(m,n_R),R_t(m,2))
        call legendrevec(p_R,xyz(1,:)/s_R(1),P_a)
        call legendrevec(p_R,xyz(2,:)/s_R(2),P_b)
        B_t=P_a(:,i_R)*P_b(:,j_R)
        R_t=matmul(B_t,c_R(1:n_R,:))
        do itarg=1,m
          if (.not.root_ok(itarg) .or. .not.all(ieee_is_finite(t0_local(:,itarg)))) t_direction(:,itarg)=R_t(itarg,:)
        enddo
        in_slab=xyz(3,:)>=zeta_lo-H_box .and. xyz(3,:)<=zeta_hi+H_box
        ! facing-edge guard: the margin across an edge scales with that edge's length (anisotropic patches)
        ell_R=[norm2(rt0),norm2(rp0)] ! physical length per unit t1, t2
        b_box(:,1)=r_s(1)+m_box*ell_R(2)/ell_R(1) ! rows box4/5/6, columns t1/t2
        b_box(:,2)=r_s(2)+m_box*ell_R(1)/ell_R(2)
        ! Off-surface targets in column4 bypass the slab and regions 4:6.
        off_column4=remaining .and. abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. h_local/=0d0
        do itarg=1,m
          if (.not.off_column4(itarg)) cycle
          if (all(ieee_is_finite(t0_local(:,itarg))) .and. ieee_is_finite(h_local(itarg))) cycle
          t0_local(:,itarg)=R_t(itarg,:)
          call legendrevec(15_8,(t0_local(:,itarg)-geometry_center)/geometry_scale,P_est)
          do k=1,3
            X_est(k)=dot_product(matmul(P_est(1,:),C_X(:,:,k)),P_est(2,:))
            J_est(k,1)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,1)),P_est(2,:))
            J_est(k,2)=dot_product(matmul(P_est(1,:),C_Xi(:,:,k,2)),P_est(2,:))
          enddo
          n_est=[J_est(2,2)*J_est(3,1)-J_est(3,2)*J_est(2,1), &
                 J_est(3,2)*J_est(1,1)-J_est(1,2)*J_est(3,1), &
                 J_est(1,2)*J_est(2,1)-J_est(2,2)*J_est(1,1)]
          n_est=n_est/norm2(n_est)
          h_local(itarg)=dot_product(tx(:,itarg)-X_est,n_est)
        enddo
        off_column4=off_column4 .and. ieee_is_finite(h_local) .and. all(ieee_is_finite(t0_local),dim=1)
        remaining=remaining .and. .not.off_column4
        inside_box4=abs(R_t(:,1))<=b_box(1,1) .and. abs(R_t(:,2))<=b_box(1,2) .and. in_slab
        inside_box5=abs(R_t(:,1))<=b_box(2,1) .and. abs(R_t(:,2))<=b_box(2,2) .and. in_slab
        inside_box6=abs(R_t(:,1))<=b_box(3,1) .and. abs(R_t(:,2))<=b_box(3,2) .and. in_slab
        ! Region 7 only above the patch: a target with trusted t0,h must also have its root inside box6
        beside=root_ok .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local) .and. &
               (abs(t0_local(1,:))>b_box(3,1) .or. abs(t0_local(2,:))>b_box(3,2))
        inside_box6=inside_box6 .and. .not.beside
        deallocate(P_a,P_b,B_t,R_t)
      endif
      region(:,1)=remaining .and. .not.inside_box4; remaining=remaining .and. inside_box4
      region(:,2)=remaining .and. .not.inside_box5; remaining=remaining .and. inside_box5
      region(:,3)=remaining .and. .not.inside_box6; remaining=remaining .and. inside_box6
      where (region(:,1)) classification=4_8
      where (region(:,2)) classification=5_8
      where (region(:,3)) classification=6_8

      ! List row heights and even column counts from the target-facing edge down.
      select case (preset)
      case (1)
        row_ratio=[1d0,1d0,1d0,3d0,3d0]
        ncols=[4,4,4,2,2]
        split_ratio=[1d0,1d0,1d0]
      case (2)
        row_ratio=[1d0,2d0,6d0]
        ncols=[4,4,2]
        split_ratio=[1d0,2d0]
      end select
      nsplit=size(split_ratio)
      panel_count(1)=sum(ncols); nlayout=panel_count(1); nparent=ncols(1)
      do layer=2,3
        nchild=2*nparent*nsplit
        panel_count(layer)=panel_count(layer-1)-nparent+nchild
        nlayout=nlayout+nchild; nparent=2*nparent
      enddo
      allocate(panels_all(4,nlayout),panel_id(maxval(panel_count),3),panel_mirror(nlayout))
      allocate(panel_pair(2,maxval(panel_count)/2,3),panel_row(nlayout),needed(nlayout),use_panel(nlayout))
      v_edges=[(1d0-2d0*sum(row_ratio(:row))/sum(row_ratio),row=0,size(row_ratio))]
      last_panel=0
      do layer=1,3
        if (layer>1) then
          nparent=ncols(1); top_height=1d0-v_edges(2)
          v_edges=[(1d0-top_height*sum(split_ratio(:row))/sum(split_ratio),row=0,nsplit)]
          ncols=spread(2*nparent,1,nsplit)
        endif
        first_panel=last_panel+1
        do row=1,size(ncols)
          width=2d0/ncols(row)
          do j=1,ncols(row)
            p=last_panel+j
            panels_all(:,p)=[-1d0+(j-1)*width,-1d0+j*width,v_edges(row+1),v_edges(row)]
            panel_mirror(p)=last_panel+ncols(row)+1-j
          enddo
          last_panel=last_panel+ncols(row)
        enddo
        nchild=last_panel-first_panel+1
        panel_id(1:nchild,layer)=[(p,p=first_panel,last_panel)]
        ! Lower rows retain their IDs from the previous layout.
        if (layer>1) panel_id(nchild+1:panel_count(layer),layer)=panel_id(nparent+1:panel_count(layer-1),layer-1)
      enddo
      ! Pair canonical left/right rectangles; geometry and targets need no symmetry.
      do layer=1,3
        is=0
        do p=1,panel_count(layer)
          ip=panel_id(p,layer)
          if (panels_all(1,ip)+panels_all(2,ip)>=0d0) cycle
          is=is+1
          panel_pair(:,is,layer)=[p,findloc(panel_id(1:panel_count(layer),layer),panel_mirror(ip),dim=1)]
        enddo
      enddo
      direction_groups(:,1)=t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,2)=t_direction(1,:)>abs(t_direction(2,:))
      direction_groups(:,3)=t_direction(2,:)<0 .and. -t_direction(2,:)>=abs(t_direction(1,:))
      direction_groups(:,4)=-t_direction(1,:)>abs(t_direction(2,:))
      rotation(:,:,1)=reshape([1d0,0d0,0d0,1d0],[2,2])
      rotation(:,:,2)=reshape([0d0,-1d0,1d0,0d0],[2,2])
      rotation(:,:,3)=-rotation(:,:,1); rotation(:,:,4)=-rotation(:,:,2)
      nq=order4(1)
      call gauss(nq,g_h,wg_h)
      call gauss(nv,g_v,wg_v)
      wq=reshape(spread(wg_h,2,nv)*spread(wg_v,1,nq),[npq])

      reverse_nodes=[(nq-mod(j-1,nq)+nq*((j-1)/nq),j=1,npq)]

      do direction=1,4
        use_panel=.false.
        do layer=1,3
          target_groups(:,layer)=region(:,layer) .and. direction_groups(:,direction)
          active(layer)=any(target_groups(:,layer))
          if (active(layer)) use_panel(panel_id(1:panel_count(layer),layer))=.true.
        enddo
        n_panel=count(use_panel)
        if (n_panel==0) cycle
        needed(1:n_panel)=pack([(p,p=1,nlayout)],use_panel)
        panel_row=0; panel_row(needed(1:n_panel))=[(p,p=1,n_panel)]
        parameter_order=maxloc(abs(rotation(:,:,direction)),dim=1)
        C_Xp=C_X; C_Xip=C_Xi
        density_column=reshape([(j,j=1,n)],[16,16])
        if (parameter_order(1)==2) then
          do k=1,3
            C_Xp(:,:,k)=transpose(C_X(:,:,k))
            do j=1,2
              C_Xip(:,:,k,j)=transpose(C_Xi(:,:,k,j))
            enddo
          enddo
          density_column=transpose(density_column)
        endif
        allocate(sq_x(npq*n_panel,3),sq_nx(npq*n_panel,3),sq_w(npq*n_panel))
        allocate(Lu(nq,16,n_panel),Lv(nv,16,n_panel),Pu(nq,16),Pv(nv,16),X(nq,nv),Xtmp(nq,16))
        allocate(rt(nq,nv,3),rp(nq,nv,3))
        do p=1,n_panel
          b=panels_all(:,needed(p)); c=(b([1,3])+b([2,4]))/2; half=(b([2,4])-b([1,3]))/2
          cols=(p-1)*npq
          node_order=[(j,j=1,npq)]
          if (b(1)+b(2)>0d0) node_order=reverse_nodes
          tu=rotation(parameter_order(1),1,direction)*(c(1)+half(1)*g_h)
          tv=rotation(parameter_order(2),2,direction)*(c(2)+half(2)*g_v)
          call legendrevec(15_8,(tu-geometry_center(parameter_order(1))) &
            /geometry_scale(parameter_order(1)),Pu)
          call legendrevec(15_8,(tv-geometry_center(parameter_order(2))) &
            /geometry_scale(parameter_order(2)),Pv)
          do k=1,3
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xp(:,:,k),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,X,nq)
            sq_x(cols+node_order,k)=reshape(X,[npq])
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,1),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rt(:,:,k),nq)
            call dgemm('N','N',nq,16_8,16_8,1d0,Pu,nq,C_Xip(:,:,k,2),16_8,0d0,Xtmp,nq)
            call dgemm('N','T',nq,nv,16_8,1d0,Xtmp,nq,Pv,nv,0d0,rp(:,:,k),nq)
          enddo
          normal_h(:,1)=reshape(rp(:,:,2)*rt(:,:,3)-rp(:,:,3)*rt(:,:,2),[npq])
          normal_h(:,2)=reshape(rp(:,:,3)*rt(:,:,1)-rp(:,:,1)*rt(:,:,3),[npq])
          normal_h(:,3)=reshape(rp(:,:,1)*rt(:,:,2)-rp(:,:,2)*rt(:,:,1),[npq])
          weight_h=sqrt(sum(normal_h**2,dim=2))
          do k=1,3
            sq_nx(cols+node_order,k)=normal_h(:,k)/weight_h
          enddo
          sq_w(cols+node_order)=product(half)*weight_h*wq
          ! Only the canonical-left panel supplies the paired density operators.
          if (b(1)+b(2)<0d0) then
            call legendrevec(15_8,tu,Pu)
            call legendrevec(15_8,tv,Pv)
            call dgemm('N','T',nq,16_8,16_8,1d0,Pu,nq,density_map,16_8,0d0,Lu(:,:,p),nq)
            Xtmp=Lu(:,:,p)
            Lu(:,1:8,p)=(Xtmp(:,1:8)+Xtmp(:,16:9:-1))/2d0
            Lu(:,9:16,p)=(Xtmp(:,1:8)-Xtmp(:,16:9:-1))/2d0
            call dgemm('N','T',nv,16_8,16_8,1d0,Pv,nv,density_map,16_8,0d0,Lv(:,:,p),nv)
          endif
        enddo
        deallocate(Pu,Pv,X,Xtmp,rt,rp)

        do layer=1,3
          if (.not.active(layer)) cycle
          mr=count(target_groups(:,layer)); rows=pack([(i,i=1,m)],target_groups(:,layer))
          target=transpose(tx(:,rows)); np=panel_count(layer)
          ! Column i+mr*(k-1)+n_row*mr*(iv-1): target i, kernel row k, short node iv.
          allocate(W_mid(n_row*mr,n),Y(8,n_row*mr*nv))
          allocate(parity_K(nq,n_row*mr*nv,2),parity_W(8,n_row*mr*16,2))
          allocate(v_node(2*npq),dlp_node(2*npq),d_node(2*npq,3),rinv_node(2*npq),dtn_node(2*npq),K_val(2*npq,n_row))
          allocate(snx_node(2*npq,3),tnx_node(2*npq,3))
          parity_W=0d0
          do p=1,np/2
            ip=panel_row(panel_id(panel_pair(1,p,layer),layer))
            ip_right=panel_row(panel_id(panel_pair(2,p,layer),layer))
            do i=1,mr
              ! The mirrored pair (entries ln and npq+ln, ln = iu+nq*(iv-1)), its kernel values, then the pair sums.
              do q_node=1,2
                off_node=(merge(ip,ip_right,q_node==1)-1)*npq; base_node=npq*(q_node-1)
                !GCC$ ivdep
                do ln=1,npq
                  dx=sq_x(off_node+ln,1)-target(i,1); dy=sq_x(off_node+ln,2)-target(i,2); dz=sq_x(off_node+ln,3)-target(i,3)
                  ie=base_node+ln
                  d_node(ie,1)=dx; d_node(ie,2)=dy; d_node(ie,3)=dz
                  distance2_inv=1d0/(dx*dx+dy*dy+dz*dz); rinv_node(ie)=sqrt(distance2_inv)
                  dtn_node(ie)=(dx*tnx(1,rows(i))+dy*tnx(2,rows(i)))+dz*tnx(3,rows(i))
                  v_node(ie)=(sq_w(off_node+ln)*inv2pi)*rinv_node(ie)
                  dlp_node(ie)=v_node(ie)*(dx*sq_nx(off_node+ln,1)+dy*sq_nx(off_node+ln,2)+dz*sq_nx(off_node+ln,3))*distance2_inv
                enddo
                if (ifDn/=0_8) snx_node(base_node+1:base_node+npq,:)=sq_nx(off_node+1:off_node+npq,:)
              enddo
              if (ifDn/=0_8) tnx_node=spread(tnx(:,rows(i)),1,2*npq)
              call compute_stok_kernels(2*npq,v_node,dlp_node,d_node,rinv_node,dtn_node,snx_node,tnx_node,K_val)
              do kk=1,n_row
                do iv=1,nv
                  cols=i+mr*(kk-1)+n_row*mr*(iv-1)
                  !GCC$ ivdep
                  do iu=1,nq
                    ln=iu+nq*(iv-1)
                    parity_K(iu,cols,1)=K_val(ln,kk)+K_val(npq+ln,kk)
                    parity_K(iu,cols,2)=K_val(ln,kk)-K_val(npq+ln,kk)
                  enddo
                enddo
              enddo
            enddo
            ! Reduce the paired coordinate to eight density nodes first.
            do mode=1,2
              call dgemm('T','N',8_8,n_row*mr*nv,nq,1d0,Lu(1,1+8*(mode-1),ip),nq, &
                parity_K(1,1,mode),nq,0d0,Y,8_8)
              call dgemm('N','N',8*n_row*mr,16_8,nv,1d0,Y,8*n_row*mr,Lv(:,:,ip),nv, &
                1d0,parity_W(1,1,mode),8*n_row*mr)
            enddo
          enddo
          do iv=1,16
            do iu=1,8
              W_mid(:,density_column(iu,iv))=parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,1) &
                +parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,2)
              W_mid(:,density_column(17-iu,iv))=parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,1) &
                -parity_W(iu,1+n_row*mr*(iv-1):n_row*mr*iv,2)
            enddo
          enddo
          ! Kernel row k in rows mr*(k-1)+1:mr*k.
          do k=1,n_kernel
            row0=mr*(k-1)
            do j=1,n
              call store_stok_weights(rows,j,j,k,W_mid(row0+1:row0+mr,j))
            enddo
          enddo
          deallocate(v_node,dlp_node,d_node,rinv_node,dtn_node,K_val,snx_node,tnx_node)
          deallocate(W_mid,Y,parity_K,parity_W,target,rows)
        enddo
        deallocate(sq_x,sq_nx,sq_w,Lu,Lv)
      enddo
    end block
  endif

  ! Region 7: tensor sinh for targets with finite source-chart coordinates.
  remaining=remaining .or. off_column4
  where (remaining) classification=7_8
  use_sinh=remaining .and. all(ieee_is_finite(t0_local),dim=1) .and. ieee_is_finite(h_local)
  if (any(use_sinh)) then
    block
      integer(8), parameter :: nlong=192,nshort=32
      integer(8) :: axis(2),n_axis(2),n_panels,it,itarg,np,first,last,l,j_sigma
      real(8) :: c(2),scale(2),J_center(3,2),a2,b2,shear,area,t(2),dc(2),center(2),beta(2)
      real(8) :: C_pos(16,16,3),Xi(16,16,3,2),Xij(16,16,3,3),orientation
      real(8) :: c_s(2),scale_s(2),Xi_s(16,16,3,2),Xij_s(16,16,3,3)
      real(8) :: c_b(2),scale_b(2),Xi_b(16,16,3,2),Xij_b(16,16,3,3)
      logical :: use_typeb,inside
      real(8) :: g32(32),wg32(32),g96(96),wg96(96),g192(192),wg192(192),sigma(16),w_sigma(16)
      real(8) :: edges(3),g_axis(nlong),wg_axis(nlong),a,b,eta0,span,step,eta,z,A0,r0,p0_axis,r0_axis
      real(8) :: alpha,ell,sgn,f,q(nlong,2),d(nlong,2),w(nlong,2),P(nlong,16,2)
      real(8) :: line(nlong*16),B_line(nlong*16,16),T1(nlong,16,2),T0L(nlong,16)
      real(8) :: P0(2,16),J0(3,2),n0(3),offset(3)
      real(8) :: dL(nlong,nshort),dS(nlong,nshort),HSS(nshort),HLS(nlong,nshort),HLL(nlong,nshort)
      real(8) :: XL(nlong,nshort,3),XS(nlong,nshort,3),chord(nlong,nshort,3),R(nlong,nshort,3)
      real(8) :: Nvec(nlong,nshort,3),distance2(nlong,nshort),rinv(nlong,nshort),omega(nlong,nshort)
      real(8) :: kernel_S(nlong,nshort),kernel_D(nlong,nshort),LL(nlong,16),LS(nshort,16)
      real(8) :: Y_leg(16,nshort),W_part7(16,16,2)
      ! K7(:,k): Stokes kernel row k at the nodes; dtn7 = chord.tnx; tnx7: tnx at every node (DLPn)
      real(8), allocatable :: K7(:,:),dtn7(:,:),tnx7(:,:)
      real(8) :: W_r7(256)
      logical :: naive
      ! Restructured evaluation (same nodes and formulas as the reference branches below).
      integer(8) :: il,js
      real(8) :: z_line(nlong),P_prev(nlong),P_cur(nlong),P_next(nlong),T0_line(nlong,16),w_T0,w_T1,rec1,rec2
      real(8) :: C_XLS(16,6,16),C_HLS(16,3,16),C_HLL(16,3,16)
      real(8) :: S_XLS(16,6,nshort),S_HLS(16,3,nshort),S_HLL(16,3,nshort)
      real(8) :: XLS(nlong,6,nshort),HLS3(nlong,3,nshort),HLL3(nlong,3,nshort),HSS3(nshort,3)
      real(8) :: d_L,d_S,hs1,hs2,hs3,r1,r2,r3,c1,c2,c3,Nv1,Nv2,Nv3,dist2,r_inv,om
      interface
        pure function expm1(x) bind(C,name='expm1') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function expm1
        pure function log1p(x) bind(C,name='log1p') result(y)
          import c_double
          real(c_double), value, intent(in) :: x
          real(c_double) :: y
        end function log1p
      end interface

      axis=[1_8,2_8]
      if (dot_product(rp0,rp0)>dot_product(rt0,rt0)) axis=[2_8,1_8]
      orientation=real(3-2*axis(1),8)
      c=geometry_center(axis); scale=geometry_scale(axis)
      J_center(:,1)=rt0; J_center(:,2)=rp0; J_center=J_center(:,axis)
      a2=dot_product(J_center(:,1),J_center(:,1)); b2=dot_product(J_center(:,2),J_center(:,2))
      shear=dot_product(J_center(:,1),J_center(:,2))
      area=norm2([rt0(2)*rp0(3)-rt0(3)*rp0(2),rt0(3)*rp0(1)-rt0(1)*rp0(3),rt0(1)*rp0(2)-rt0(2)*rp0(1)])
      C_pos=C_X
      Xi=C_Xi(:,:,:,axis); Xij=C_Xij(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
      if (axis(1)==2) then
        do k=1,3
          C_pos(:,:,k)=transpose(C_pos(:,:,k))
          do j=1,2
            Xi(:,:,k,j)=transpose(Xi(:,:,k,j))
          enddo
          do j=1,3
            Xij(:,:,k,j)=transpose(Xij(:,:,k,j))
          enddo
        enddo
      endif
      c_s=c; scale_s=scale; Xi_s=Xi; Xij_s=Xij ! source chart
      use_typeb=all(geometry_scale_typeb>0d0)
      if (use_typeb) then ! type-b chart: the source sits inside it; used when t0 falls inside
        c_b=geometry_center_typeb(axis); scale_b=geometry_scale_typeb(axis)
        Xi_b=C_Xi_typeb(:,:,:,axis); Xij_b=C_Xij_typeb(:,:,:,[2*axis(1)-1,2_8,2*axis(2)-1])
        if (axis(1)==2) then
          do k=1,3
            do j=1,2
              Xi_b(:,:,k,j)=transpose(Xi_b(:,:,k,j))
            enddo
            do j=1,3
              Xij_b(:,:,k,j)=transpose(Xij_b(:,:,k,j))
            enddo
          enddo
        endif
      endif
      call gauss(32_8,g32,wg32)
      call gauss(96_8,g96,wg96)
      call gauss(192_8,g192,wg192)
      sigma=(g16+1d0)/2; w_sigma=wg16/2
      rows=pack([(j,j=1,m)],use_sinh)
      n_axis=[nlong,nshort]
      allocate(K7(nlong*nshort,n_row),dtn7(nlong,nshort),tnx7(nlong*nshort,3))

      do it=1,size(rows,kind=8)
        itarg=rows(it)
        t=t0_local(axis,itarg)
        ! naive=abs(t(2))>2.5d0
        ! naive=off_column4(itarg) .or. abs(t(2))>1.5d0 ! heuristic
        inside=.false.
        if (use_typeb) inside=all(abs(t-c_b)<=scale_b*(1d0+1d-12)) ! t0 inside the type-b chart: geometry is valid at t0
        if (inside .and. .not.off_column4(itarg) .and. root_ok(itarg)) then
          c=c_b; scale=scale_b; Xi=Xi_b; Xij=Xij_b; naive=.false.
        else
          c=c_s; scale=scale_s; Xi=Xi_s; Xij=Xij_s; naive=.not.root_ok(itarg) .or. off_column4(itarg) .or. abs(t(2))>1.5d0
        endif
        dc=min(1d0,max(-1d0,t))-t
        center=t-[shear/a2*dc(2),shear/b2*dc(1)]
        beta=[hypot(h_local(itarg),area/sqrt(a2)*dc(2))/sqrt(a2),hypot(h_local(itarg),area/sqrt(b2)*dc(1))/sqrt(b2)]
        if (.not.root_ok(itarg) .and. (any(.not.ieee_is_finite(center)) .or. any(.not.ieee_is_finite(beta)) .or. &
            any(beta==0d0 .and. abs(center)<=1d0))) then
          use_sinh(itarg)=.false.
          cycle
        endif

        do k=1,2
          n_panels=1; edges(1:2)=[-1d0,1d0]
          if (k==2) then
            np=nshort; g_axis(:np)=g32; wg_axis(:np)=wg32
          else if (abs(center(k))<1d0) then
            np=nlong/2; g_axis(:np)=g96; wg_axis(:np)=wg96
            n_panels=2; edges=[-1d0,center(k),1d0]
          else
            np=nlong; g_axis=g192; wg_axis=wg192
          endif
          do j=1,n_panels
            a=edges(j); b=edges(j+1)
            first=(j-1)*np+1; last=j*np
            if (beta(k)>0d0) then
              eta0=asinh((a-center(k))/beta(k))
              span=asinh((b-center(k))/beta(k))-eta0
              if (abs(center(k))>=2d0) then
                span=2*atanh((b-a)/(hypot(a-center(k),beta(k))+hypot(b-center(k),beta(k))))
              endif
              A0=a-center(k); r0=hypot(A0,beta(k))
              if (A0>=0d0) then
                p0_axis=r0+A0; r0_axis=beta(k)*(beta(k)/p0_axis)
              else
                r0_axis=r0-A0; p0_axis=beta(k)*(beta(k)/r0_axis)
              endif
              do i=first,last
                step=(g_axis(i-first+1)+1d0)*(span/2)
                eta=eta0+step; z=beta(k)*sinh(eta)
                w(i,k)=wg_axis(i-first+1)*((span/2)*beta(k)*cosh(eta))
                if (abs(center(k))<=1d0) then
                  q(i,k)=center(k)+z
                else
                  q(i,k)=a+(p0_axis*expm1(step)-r0_axis*expm1(-step))/2
                endif
                d(i,k)=(center(k)-t(k))+z
              enddo
            else
              alpha=abs(center(k)); sgn=sign(1d0,center(k))
              ell=-log1p(2/(alpha-1d0))
              do i=first,last
                f=(sgn*g_axis(i-first+1)+1d0)/2
                z=-sgn*(alpha+1d0)*exp(f*ell)
                q(i,k)=sgn*(-1d0-(alpha+1d0)*expm1(f*ell))
                w(i,k)=wg_axis(i-first+1)*(-abs(z)*ell/2)
                d(i,k)=(center(k)-t(k))+z
              enddo
            endif
          enddo
          np=n_axis(k)
          call legendrevec(15_8,(q(:np,k)-c(k))/scale(k),P(:np,:,k))
          if (.not.naive) then
            ! Change .false. to .true. (here and below) to force the original evaluation for comparison.
            if (.false.) then
            do j_sigma=1,16
              line((j_sigma-1)*np+1:j_sigma*np)=(t(k)+d(:np,k)*sigma(j_sigma)-c(k))/scale(k)
            enddo
            call legendrevec(15_8,line(:np*16),B_line(:np*16,:))
            do l=1,16
              T1(:np,l,k)=matmul(reshape(B_line(:np*16,l),[np,16_8]),w_sigma*sigma)
              if (k==1) T0L(:,l)=matmul(reshape(B_line(:,l),[nlong,16_8]),w_sigma)
            enddo
            else
              ! Same GL16 sums over sigma, accumulated inside the Legendre recurrence:
              ! T1(i,l) = sum_j w_sigma(j)*sigma(j)*P_{l-1}(z_ij), T0(i,l) = sum_j w_sigma(j)*P_{l-1}(z_ij).
              T1(:np,:,k)=0d0; T0_line(:np,:)=0d0
              do j_sigma=1,16
                w_T1=w_sigma(j_sigma)*sigma(j_sigma); w_T0=w_sigma(j_sigma)
                do i=1,np
                  z_line(i)=(t(k)+d(i,k)*sigma(j_sigma)-c(k))/scale(k)
                  P_prev(i)=1d0; P_cur(i)=z_line(i)
                  T1(i,1,k)=T1(i,1,k)+w_T1; T1(i,2,k)=T1(i,2,k)+P_cur(i)*w_T1
                  T0_line(i,1)=T0_line(i,1)+w_T0; T0_line(i,2)=T0_line(i,2)+P_cur(i)*w_T0
                enddo
                do l=1,14
                  rec1=real(2*l+1,8)/real(l+1,8); rec2=real(l,8)/real(l+1,8)
                  !GCC$ ivdep
                  do i=1,np
                    P_next(i)=rec1*z_line(i)*P_cur(i)-rec2*P_prev(i)
                    T1(i,l+2,k)=T1(i,l+2,k)+P_next(i)*w_T1
                    T0_line(i,l+2)=T0_line(i,l+2)+P_next(i)*w_T0
                    P_prev(i)=P_cur(i); P_cur(i)=P_next(i)
                  enddo
                enddo
              enddo
              if (k==1) T0L=T0_line
            endif
          endif
        enddo

        ! Naive targets keep the original evaluation; .false. -> .true. forces it for every target.
        if (.false. .or. naive) then
        do k=1,3
          XL(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,1)),transpose(P(:nshort,:,2)))
          XS(:,:,k)=matmul(matmul(P(:,:,1),Xi(:,:,k,2)),transpose(P(:nshort,:,2)))
        enddo
        if (naive) then
          do k=1,3
            chord(:,:,k)=matmul(matmul(P(:,:,1),C_pos(:,:,k)),transpose(P(:nshort,:,2)))-tx(k,itarg)
          enddo
          R=chord
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          dL=spread(d(:,1),2,nshort); dS=spread(d(:nshort,2),1,nlong)
          do k=1,3
            HSS=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
            HLS=matmul(matmul(T0L,Xij(:,:,k,2)),transpose(P(:nshort,:,2)))
            HLL=matmul(matmul(T1(:,:,1),Xij(:,:,k,1)),transpose(P(:nshort,:,2)))
            R(:,:,k)=offset(k)-dS**2*spread(HSS,1,nlong)-dL*dS*HLS-dL**2*HLL
            chord(:,:,k)=dL*XL(:,:,k)+dS*XS(:,:,k)+R(:,:,k)
          enddo
        endif

        Nvec(:,:,1)=orientation*(XS(:,:,2)*XL(:,:,3)-XS(:,:,3)*XL(:,:,2))
        Nvec(:,:,2)=orientation*(XS(:,:,3)*XL(:,:,1)-XS(:,:,1)*XL(:,:,3))
        Nvec(:,:,3)=orientation*(XS(:,:,1)*XL(:,:,2)-XS(:,:,2)*XL(:,:,1))
        distance2=sum(chord**2,dim=3)
        rinv=1d0/sqrt(distance2)
        omega=inv2pi*spread(w(:,1),2,nshort)*spread(w(:nshort,2),1,nlong)
        kernel_S=omega*sqrt(sum(Nvec**2,dim=3))*rinv
        kernel_D=omega*sum(Nvec*R,dim=3)*rinv/distance2
        else
          call legendrevec(15_8,(t-c)/scale,P0)
          do k=1,3
            J0(k,1)=dot_product(matmul(P0(1,:),Xi(:,:,k,1)),P0(2,:))
            J0(k,2)=dot_product(matmul(P0(1,:),Xi(:,:,k,2)),P0(2,:))
          enddo
          n0=orientation*[J0(2,2)*J0(3,1)-J0(3,2)*J0(2,1),J0(3,2)*J0(1,1)-J0(1,2)*J0(3,1), &
                         J0(1,2)*J0(2,1)-J0(2,2)*J0(1,1)]
          offset=-h_local(itarg)*n0/norm2(n0)
          ! Coefficients stacked as (long mode, field, short mode): XL (fields 1:3) and XS (4:6); HLS; HLL.
          do l=1,16
            do k=1,3
              C_XLS(:,k,l)=Xi(:,l,k,1); C_XLS(:,3+k,l)=Xi(:,l,k,2)
              C_HLS(:,k,l)=Xij(:,l,k,2); C_HLL(:,k,l)=Xij(:,l,k,1)
            enddo
          enddo
          ! X = P_long*(C*P_short^T): the short direction first, then one GEMM per group along the long direction.
          call dgemm('N','T',96_8,nshort,16_8,1d0,C_XLS,96_8,P(1,1,2),nlong,0d0,S_XLS,96_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLS,48_8,P(1,1,2),nlong,0d0,S_HLS,48_8)
          call dgemm('N','T',48_8,nshort,16_8,1d0,C_HLL,48_8,P(1,1,2),nlong,0d0,S_HLL,48_8)
          call dgemm('N','N',nlong,6*nshort,16_8,1d0,P(1,1,1),nlong,S_XLS,16_8,0d0,XLS,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T0L,nlong,S_HLS,16_8,0d0,HLS3,nlong)
          call dgemm('N','N',nlong,3*nshort,16_8,1d0,T1(1,1,1),nlong,S_HLL,16_8,0d0,HLL3,nlong)
          do k=1,3
            HSS3(:,k)=matmul(matmul(P0(1,:),Xij(:,:,k,3)),transpose(T1(:nshort,:,2)))
          enddo
          ! R, chord, normal and both kernels at each node, same expressions as the reference.
          do js=1,nshort
            d_S=d(js,2)
            hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
            !GCC$ ivdep
            do il=1,nlong
              d_L=d(il,1)
              r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
              r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
              r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
              c1=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
              c2=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
              c3=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
              Nv1=orientation*(XLS(il,5,js)*XLS(il,3,js)-XLS(il,6,js)*XLS(il,2,js))
              Nv2=orientation*(XLS(il,6,js)*XLS(il,1,js)-XLS(il,4,js)*XLS(il,3,js))
              Nv3=orientation*(XLS(il,4,js)*XLS(il,2,js)-XLS(il,5,js)*XLS(il,1,js))
              dist2=c1**2+c2**2+c3**2
              r_inv=1d0/sqrt(dist2)
              om=inv2pi*w(il,1)*w(js,2)
              kernel_S(il,js)=om*sqrt(Nv1**2+Nv2**2+Nv3**2)*r_inv
              kernel_D(il,js)=om*(Nv1*r1+Nv2*r2+Nv3*r3)*r_inv/dist2
            enddo
          enddo
          ! chord and distance2 at every node: the Stokes kernels need the chord.
            do js=1,nshort
              d_S=d(js,2)
              hs1=d_S**2*HSS3(js,1); hs2=d_S**2*HSS3(js,2); hs3=d_S**2*HSS3(js,3)
              do il=1,nlong
                d_L=d(il,1)
                r1=offset(1)-hs1-d_L*d_S*HLS3(il,1,js)-d_L**2*HLL3(il,1,js)
                r2=offset(2)-hs2-d_L*d_S*HLS3(il,2,js)-d_L**2*HLL3(il,2,js)
                r3=offset(3)-hs3-d_L*d_S*HLS3(il,3,js)-d_L**2*HLL3(il,3,js)
                chord(il,js,1)=d_L*XLS(il,1,js)+d_S*XLS(il,4,js)+r1
                chord(il,js,2)=d_L*XLS(il,2,js)+d_S*XLS(il,5,js)+r2
                chord(il,js,3)=d_L*XLS(il,3,js)+d_S*XLS(il,6,js)+r3
                distance2(il,js)=chord(il,js,1)**2+chord(il,js,2)**2+chord(il,js,3)**2
              enddo
            enddo
            if (ifDn/=0_8) then   ! the normal at every node: DLPn needs the vector
              do js=1,nshort
                do il=1,nlong
                  Nvec(il,js,1)=orientation*(XLS(il,5,js)*XLS(il,3,js)-XLS(il,6,js)*XLS(il,2,js))
                  Nvec(il,js,2)=orientation*(XLS(il,6,js)*XLS(il,1,js)-XLS(il,4,js)*XLS(il,3,js))
                  Nvec(il,js,3)=orientation*(XLS(il,4,js)*XLS(il,2,js)-XLS(il,5,js)*XLS(il,1,js))
                enddo
              enddo
            endif
        endif

        ! Stokes kernels at the nodes from the Laplace SLP, DLP values (kernel_D: the stable normal remainder) and the chord.
        rinv=1d0/sqrt(distance2)
        dtn7=(chord(:,:,1)*tnx(1,itarg)+chord(:,:,2)*tnx(2,itarg))+chord(:,:,3)*tnx(3,itarg)
        if (ifDn/=0_8) then   ! unit normal, as sq_nx at the other regions
          Nvec=Nvec/spread(sqrt(sum(Nvec**2,dim=3)),3,3)
          tnx7=spread(tnx(:,itarg),1,nlong*nshort)
        endif
        call compute_stok_kernels(nlong*nshort,kernel_S,kernel_D,chord,rinv,dtn7,Nvec,tnx7,K7)

        call legendrevec(15_8,q(:,1),P(:,:,1))
        call legendrevec(15_8,q(:nshort,2),P(:nshort,:,2))
        LL=matmul(P(:,:,1),transpose(density_map)); LS=matmul(P(:nshort,:,2),transpose(density_map))
        do k=1,n_kernel
          call dgemm('T','N',16_8,nshort,nlong,1d0,LL,nlong,K7(1,k),nlong,0d0,Y_leg,16_8)
          call dgemm('N','N',16_8,16_8,nshort,1d0,Y_leg,16_8,LS,nshort,0d0,W_part7(1,1,1),16_8)
          if (axis(1)==2) W_part7(:,:,1)=transpose(W_part7(:,:,1))
          W_r7=reshape(W_part7(:,:,1),[256])
          call store_stok_weights([itarg],1_8,n,k,W_r7)
        enddo
      enddo
      deallocate(rows)
      deallocate(K7,dtn7,tnx7)
    end block
  endif
  remaining=remaining .and. .not.use_sinh
  where (remaining) ier=1_8


contains

  ! Stokes kernel rows at nb nodes from the Laplace SLP value v_in (weight/(2*pi*r)), DLP value dlp_in, the chord
  ! d_in = y-x, 1/r and dtn_in = d_in.tnx: K_out(:,k_S+c) SLP, K_out(:,k_D+c) DLP, K_out(:,k_Sn+c) SLPn, c = 1..6;
  ! K_out(:,k_Dn+i+3*(j-1)) DLPn from the unit source normal snx_in and the target normal tnx_in (read for DLPn only).
  subroutine compute_stok_kernels(nb,v_in,dlp_in,d_in,rinv_in,dtn_in,snx_in,tnx_in,K_out)
    integer(8), intent(in) :: nb
    real(8), intent(in) :: v_in(nb),dlp_in(nb),d_in(nb,3),rinv_in(nb),dtn_in(nb),snx_in(nb,3),tnx_in(nb,3)
    real(8), intent(out) :: K_out(nb,n_row)
    integer(8) :: c_s,i_s,ic_s,jc_s
    real(8) :: delta,dd,rinv2,rinv2_s(nb),v_r2(nb),dlp_r2(nb),dtn_r2(nb),tns(nb)
    do c_s=1,6
      delta=merge(1d0,0d0,a_comp(c_s)==b_comp(c_s))
      !GCC$ ivdep
      do i_s=1,nb
        rinv2=rinv_in(i_s)*rinv_in(i_s)
        dd=d_in(i_s,a_comp(c_s))*d_in(i_s,b_comp(c_s))*rinv2
        if (ifS/=0_8) K_out(i_s,k_S+c_s)=0.25d0*v_in(i_s)*(delta+dd)
        if (ifD/=0_8) K_out(i_s,k_D+c_s)=1.5d0*dlp_in(i_s)*dd
        if (ifSn/=0_8) K_out(i_s,k_Sn+c_s)=1.5d0*v_in(i_s)*rinv2*dtn_in(i_s)*dd
      enddo
    enddo
    if (ifDn==0_8) return
    ! per node, then per component: v/r^2, dlp/r^2, dtn/r^2, tnx.n_s
    !GCC$ ivdep
    do i_s=1,nb
      rinv2_s(i_s)=rinv_in(i_s)*rinv_in(i_s)
      v_r2(i_s)=v_in(i_s)*rinv2_s(i_s)
      dlp_r2(i_s)=dlp_in(i_s)*rinv2_s(i_s)
      dtn_r2(i_s)=dtn_in(i_s)*rinv2_s(i_s)
      tns(i_s)=(tnx_in(i_s,1)*snx_in(i_s,1)+tnx_in(i_s,2)*snx_in(i_s,2))+tnx_in(i_s,3)*snx_in(i_s,3)
    enddo
    do jc_s=1,3
      do ic_s=1,3
        c_s=ic_s+3*(jc_s-1)
        delta=merge(1d0,0d0,ic_s==jc_s)
        !GCC$ ivdep
        do i_s=1,nb
          dd=d_in(i_s,ic_s)*d_in(i_s,jc_s)*rinv2_s(i_s)
          K_out(i_s,k_Dn+c_s)=-(v_r2(i_s)*((tnx_in(i_s,ic_s)*snx_in(i_s,jc_s)+1.5d0*tns(i_s)*dd) &
            +1.5d0*dtn_r2(i_s)*snx_in(i_s,ic_s)*d_in(i_s,jc_s)) &
            +dlp_r2(i_s)*(1.5d0*d_in(i_s,ic_s)*tnx_in(i_s,jc_s)+dtn_in(i_s)*(1.5d0*delta-15.0d0*dd)))
        enddo
      enddo
    enddo
  end subroutine compute_stok_kernels

  ! Weights of kernel row k for the targets target_rows and source nodes j_first:j_last into WS, WD or WSn
  ! (component (a_comp,b_comp) and its symmetric copy) or WDn (component (i,j), k = k_Dn+i+3*(j-1));
  ! W_k(size(target_rows)*(j_last-j_first+1)), targets fastest.
  subroutine store_stok_weights(target_rows,j_first,j_last,k,W_k)
    integer(8), intent(in) :: target_rows(:),j_first,j_last,k
    real(8), intent(in) :: W_k(size(target_rows),j_first:j_last)
    integer(8) :: c_s
    c_s=mod(k-1_8,6_8)+1_8
    if (k<=k_D) then
      WS(target_rows,j_first:j_last,a_comp(c_s),b_comp(c_s))=W_k
      WS(target_rows,j_first:j_last,b_comp(c_s),a_comp(c_s))=W_k
    else if (k<=k_Sn) then
      WD(target_rows,j_first:j_last,a_comp(c_s),b_comp(c_s))=W_k
      WD(target_rows,j_first:j_last,b_comp(c_s),a_comp(c_s))=W_k
    else if (k<=k_Dn) then
      WSn(target_rows,j_first:j_last,a_comp(c_s),b_comp(c_s))=W_k
      WSn(target_rows,j_first:j_last,b_comp(c_s),a_comp(c_s))=W_k
    else
      c_s=k-k_Dn
      WDn(target_rows,j_first:j_last,mod(c_s-1_8,3_8)+1_8,(c_s-1_8)/3_8+1_8)=W_k
    endif
  end subroutine store_stok_weights

end subroutine quadSD_stok_mid_w16

