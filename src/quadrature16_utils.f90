module surfquadr
  ! Shared quadrature steps for the eight near/self drivers.
  implicit none
  private
  public :: compute_legendre_values, compute_gauss_legendre, compute_legendre_translation
  public :: source_geometry, compute_source_geometry
  public :: target_geometry, compute_target_geometry
  public :: edge_tables, compute_edge_rule
  public :: radial_tables, radial_expansion_tables, compute_radial_expansion
  public :: displacement_expansion, compute_displacement_expansion
  public :: jacobian_expansion, compute_jacobian_expansion
  public :: distance_factor, compute_distance_factor
  public :: kernel_weights, compute_kernelradial_S_on, compute_kernelradial_D_on, compute_kernelradial_Sn_on
  public :: compute_kernelradial_S_helm_on, compute_kernelradial_D_helm_on, compute_kernelradial_Sn_helm_on, &
    compute_kernelradial_Dn_helm_on
  public :: compute_kernelradial_Dn_on, compute_kernelradial_D_off, compute_kernelradial_Sn_off, compute_kernelradial_Dn_off
  public :: compute_kernelradial_S_helm_off, compute_kernelradial_D_helm_off, compute_kernelradial_Sn_helm_off, &
    compute_kernelradial_Dn_helm_off, compute_singular_correction_helm_off
  public :: compute_kernelradial_S_stok_on, compute_kernelradial_D_stok_on, compute_kernelradial_Sn_stok_on, &
    compute_kernelradial_Dn_stok_on
  public :: compute_kernelradial_S_stok_off, compute_kernelradial_Dn_stok_off, compute_singular_correction_stok_off
  public :: compute_singular_correction_on, compute_singular_correction_off, compute_muhat_weights, compute_final_weights

  integer(8), parameter :: nr=24_8      ! radial Gauss order (fixed)
  real(8), parameter :: pi=3.141592653589793238462643383279502884d0
  ! Legendre recurrence coefficients, DLMF 18.9.2: x*p_n = a_n*p_{n+1}+b_n*p_n+c_n*p_{n-1}, a_n = (n+1)/(2n+1),
  ! b_n = 0, c_n = n/(2n+1); a_n_inv = 1/a_n. n_index: implied-do index of these constructors only.
  integer(8) :: n_index
  real(8), parameter :: a_n(0:256)=[(real(n_index+1,8)/real(2*n_index+1,8), n_index=0,256)]
  real(8), parameter :: b_n(0:256)=0.0d0
  real(8), parameter :: c_n(0:256)=[(real(n_index,8)/real(2*n_index+1,8), n_index=0,256)]
  real(8), parameter :: a_n_inv(0:256)=[(1.0d0/a_n(n_index), n_index=0,256)]

  ! Edge Gauss rules with order(r) nodes on [-1,1], r=1..32.
  type :: edge_tables
    real(8) :: budget
    integer(8) :: order(32)
    real(8) :: nodes(256,32), weights(256,32)
  end type edge_tables

  ! Radial Gauss rule on [0,1] and the tables computed by prepare_tables_w16.
  type :: radial_tables
    real(8) :: lambda(nr), lambda_inv(nr), w(nr)   ! nodes, 1/nodes, weights
    real(8) :: T1(nr,nr), T2(nr,nr)        ! c1r = T1*V1r (order 0, 1); cr = T2*(d.Gr) (order 2)
    ! Off surface (step 8): Lj = lambda^(0:nr+11) (p=3 reads to nr+11); Lj(:,1:nr)^T = U_V*diag(1/S_inv)*W_V^T;
    ! ell0 = Lagrange basis at lambda=0 and its eps-expansion weights; binom: Pascal; c_y(m+1,p+1) = binomial(-(2p+1)/2,m).
    real(8) :: Lj(nr,nr+12), U_V(nr,nr), S_inv(nr), W_V(nr,nr)
    real(8) :: ell0(nr), ell0_coeff(nr,4), ell0_over_lambda(nr), radial_inv_sums(6), ell0_inv_sums(2,4)
    real(8) :: binom(7,7), c_y(7,4)
  end type radial_tables

  ! Step 0. Source-patch geometry shared by all targets.
  type :: source_geometry
    real(8) :: C_Xi_fixed(16,6,16,2)  ! (fixed-axis mode, r_u 1:3 | r_v 4:6, moving-axis mode, fixed axis)
    logical :: if_har                 ! patch-center aspect ratio > 8
  end type source_geometry

  ! Step 1. Surface geometry at the projected target.
  type :: target_geometry
    real(8) :: t0(2), h
    real(8) :: P0(16,2)               ! Legendre values at u' (column 1) and v' (column 2)
    real(8) :: dP0(16,2)              ! their derivatives
    real(8) :: r_u0(3), r_v0(3)
    real(8) :: metric(2,2)
    real(8) :: N0(3), J0              ! N0 = r_u0 x r_v0, J0 = ||N0||
    real(8) :: x_E(4)                 ! (1+v', 1-u', 1-v', 1+u')
    real(8) :: r_uu0(3), r_uv0(3), r_vv0(3)   ! off surface only
    integer(8) :: K                   ! off surface only: order of step 8's expansion in d_curv
  end type target_geometry

  ! Step 4. Displacement expansion along the rays, with (nr,3,n_E) = (node, component, ray).
  type :: displacement_expansion
    real(8), allocatable :: C0(:,:), q0(:), c1r(:,:,:)   ! (3,n_E), (n_E), (nr,3,n_E)
    real(8), allocatable :: V1r(:,:,:)                   ! (nr,3,n_E): order 0, 1
    real(8), allocatable :: c1(:,:,:), cr(:,:,:)         ! (1,3,n_E), (nr,3,n_E): order 2
  end type displacement_expansion

  ! Step 5. Surface normal and Jacobian expansions along the rays.
  type :: jacobian_expansion
    real(8) :: N0(1,3), J0
    real(8), allocatable :: N(:,:,:), J(:,:)              ! (nr,3,n_E) if if_N, otherwise (nr,n_E)
    real(8), allocatable :: N1(:,:,:), J1(:)              ! (1,3,n_E), (n_E): order >= 1
    real(8), allocatable :: Nr(:,:,:), J1r(:,:), Jr(:,:)  ! (nr,3,n_E), (nr,n_E), (nr,n_E): order >= 1
  end type jacobian_expansion

  ! Step 6. Distance factors and their radial expansions.
  type :: distance_factor
    real(8), allocatable :: g0(:), g(:,:)                ! (n_E), (nr,n_E)
    real(8), allocatable :: g1r(:,:)                     ! (nr,n_E): order >= 1
    real(8), allocatable :: g1(:), gr(:,:), beta1(:)     ! (n_E), (nr,n_E), (n_E): order 2
    real(8), allocatable :: d_curv(:,:)                  ! (nr,n_E): off surface, for step 8
  end type distance_factor

  ! Steps 7-8. Weights on the density samples of one edge.
  type :: kernel_weights
    real(8), allocatable :: weight_lambda(:,:,:)   ! (nr,n_E,ncomp): on mu(t'+lambda*d) (SLP, DLP, SLPn)
    real(8), allocatable :: weight_r(:,:,:)        ! (nr,n_E,ncomp): on mu's order-th remainder (SLPn mu1r, DLPn mur)
  end type kernel_weights

  ! Step 3 tables for the current edge, reused in step 9.
  type :: radial_expansion_tables
    real(8), allocatable :: T(:,:,:), Pfixed(:,:), Pedge(:,:)         ! (nr,16,16), (nr,16), (n_E,16)
    real(8) :: Pf0(16), Pm0(16), Pm1(16), Pf1(16)                     ! order>=1 or if_har
    real(8), allocatable :: T_dd(:,:,:), Pfr(:,:), Pf1r(:,:)         ! (nr,16,16), (nr,16), (nr,16)
  end type radial_expansion_tables

  ! Step 8, off surface. compute_moment_correction_off stores the corrected rays (closed-form eps first)
  ! and their V^-1 solves. For ray rays(e), term m occupies column (K_m+1)*(e-1)+m+1.
  type :: moment_correction
    integer(8) :: K_m
    real(8) :: epsilon_big
    real(8), allocatable :: epsilon(:)                                           ! (n_E)
    integer(8), allocatable :: rays(:)                                           ! (corrected rays)
    real(8), allocatable :: sigma(:,:)                                           ! (nr,ncol): p=0
    real(8), allocatable :: omega(:,:), z(:,:), alpha(:), rest0(:), h_scale(:)   ! (nr,ncol), (ncol): p=1
    real(8), allocatable :: omega_shift(:,:,:,:)                     ! (nr,K_m+1,3,n_E): p=2; (nr,K_m+1,K_m+7,n_E): p=3
    real(8), allocatable :: omega_shift_gradient(:,:,:,:,:)           ! (nr,K_m+1,3,3,n_E): p=2
    real(8), allocatable :: DeltaM01(:,:,:), DeltaM01_gradient(:,:,:,:)            ! (2,K_m+1,n_E), (2,K_m+1,3,n_E): p=2
  end type moment_correction

contains

  subroutine compute_legendre_values(deg,x,d,p,dp)
    ! Legendre polynomials p_n = P_n and their derivatives along the scalar direction d at x, n = 0..deg (1 <= deg <= 256).
    !   p(i,n) = p_n(x_i); dp(i,n) = d*p_n'(x_i), formed only when present; d defaults to 1. Rows 1..size(x) are set;
    !   x, p, dp contiguous (pass whole arrays: no copy).
    ! DLMF 18.9.2 (a_n, b_n, c_n, a_n_inv: module parameters) solved for p_{n+1}:
    !   p_{n+1}  = ((x-b_n)*p_n-c_n*p_{n-1})*(1/a_n),           p_0 = 1, p_{-1} = 0,
    !   dp_{n+1} = ((x-b_n)*dp_n+d*p_n-c_n*dp_{n-1})*(1/a_n)    (its derivative along d).
    integer(8), intent(in) :: deg
    real(8), intent(in), contiguous :: x(:)
    real(8), intent(in), optional :: d
    real(8), intent(inout), contiguous :: p(:,0:)
    real(8), intent(inout), optional, contiguous :: dp(:,0:)
    real(8) :: d_value
    integer(8) :: n, n_x

    n_x=size(x,kind=8)
    p(1:n_x,0)=1.0d0
    p(1:n_x,1)=((x-b_n(0))*p(1:n_x,0))*a_n_inv(0)                                ! p_{-1} = 0
    do n=1,deg-1
      p(1:n_x,n+1)=((x-b_n(n))*p(1:n_x,n)-c_n(n)*p(1:n_x,n-1))*a_n_inv(n)
    enddo
    if (.not.present(dp)) return
    d_value=1.0d0
    if (present(d)) d_value=d
    dp(1:n_x,0)=0.0d0
    dp(1:n_x,1)=((x-b_n(0))*dp(1:n_x,0)+d_value*p(1:n_x,0))*a_n_inv(0)          ! dp_{-1} = 0
    do n=1,deg-1
      dp(1:n_x,n+1)=(((x-b_n(n))*dp(1:n_x,n)+d_value*p(1:n_x,n))-c_n(n)*dp(1:n_x,n-1))*a_n_inv(n)
    enddo
  end subroutine compute_legendre_values

  subroutine compute_gauss_legendre(n,x,w)
    ! n-point Gauss-Legendre rule on [-1,1], 2 <= n <= 256, with ascending nodes x_k and weights w_k.
    ! DLMF 3.5(v): the x_k are the zeros of p_n = P_n, and the rule is exact for polynomials of degree <= 2n-1.
    !   Symmetry: x_{n+1-k} = -x_k, w_{n+1-k} = w_k; for odd n the middle node is 0.
    !   Nodes k <= n/2: Newton's method, x <- x-p_n(x)/p_n'(x), from x_k = -cos((k-1/4)*pi/(n+1/2)), until no node
    !   moves by more than eps; p_{n-1}, p_n by DLMF 18.9.2 as in compute_legendre_values, and
    !   p_n'(x) = n*(p_{n-1}(x)-x*p_n(x))/(1-x^2) (DLMF 14.10.5, mu = 0).
    !   Weights: by DLMF (3.5.20_2), sum_k p_l(x_k)*p_m(x_k)*w_k = h_m*delta_lm for l, m <= n-1, with h_m = 2/(2m+1)
    !   (DLMF Table 18.3.1), so 1/w_k = sum_{m=0}^{n-1} (m+1/2)*p_m(x_k)^2.
    integer(8), intent(in) :: n
    real(8), intent(inout) :: x(n), w(n)
    real(8) :: p_mm1((n+1)/2), p_m((n+1)/2), p_mp1((n+1)/2), x_new(n/2), step
    integer(8) :: k, m, iteration, n_half, n_w

    n_half=n/2
    n_w=(n+1)/2                    ! nodes x <= 0
    do k=1,n_half
      x(k)=-cos((real(k,8)-0.25d0)*pi/(real(n,8)+0.5d0))
    enddo
    do iteration=1,10
      p_mm1(1:n_half)=0.0d0        ! p_{-1}
      p_m(1:n_half)=1.0d0          ! p_0
      do m=0,n-1
        p_mp1(1:n_half)=((x(1:n_half)-b_n(m))*p_m(1:n_half)-c_n(m)*p_mm1(1:n_half))*a_n_inv(m)
        p_mm1(1:n_half)=p_m(1:n_half)
        p_m(1:n_half)=p_mp1(1:n_half)
      enddo
      ! p_m = p_n, p_mm1 = p_{n-1}
      x_new=x(1:n_half)-p_m(1:n_half)/(real(n,8)*(p_mm1(1:n_half)-x(1:n_half)*p_m(1:n_half)) &
        /((1.0d0-x(1:n_half))*(1.0d0+x(1:n_half))))
      step=maxval(abs(x_new-x(1:n_half)))
      x(1:n_half)=x_new
      if (step<=epsilon(1.0d0)) exit
    enddo
    if (n_w>n_half) x(n_w)=0.0d0
    p_mm1=0.0d0
    p_m=1.0d0
    w(1:n_w)=0.0d0
    do m=0,n-1
      w(1:n_w)=w(1:n_w)+(real(m,8)+0.5d0)*p_m**2
      p_mp1=((x(1:n_w)-b_n(m))*p_m-c_n(m)*p_mm1)*a_n_inv(m)
      p_mm1=p_m
      p_m=p_mp1
    enddo
    w(1:n_w)=1.0d0/w(1:n_w)
    do k=1,n_half
      x(n+1-k)=-x(k)
      w(n+1-k)=w(k)
    enddo
  end subroutine compute_gauss_legendre

  subroutine compute_legendre_translation(deg,x0,lambda,T,T_dd)
    ! Legendre coefficients of p_n(y), y = x0+lambda*(x-x0), n = 0..deg (2 <= deg <= 256).
    !   p_n(y) = sum_k T(i,k,n)*p_k(x);   with T_dd:  p_n(y) = p_n(x0)+(y-x0)*p_n'(x0)+(y-x0)^2*d_n(y),
    !   d_n(y) = sum_k T_dd(i,k,n)*p_k(x).   p_n = P_n (Legendre), i = 1..size(lambda); T_dd formed only when present.
    ! DLMF 18.9.2 (a_n, b_n, c_n, a_n_inv: module parameters), x*p_n = a_n*p_{n+1}+b_n*p_n+c_n*p_{n-1}, used twice:
    !   at y, p_{n+1}(y) = ((y-b_n)*p_n(y)-c_n*p_{n-1}(y))*(1/a_n), on coefficient vectors in p_k(x);
    !   and for x times a coefficient vector v (x_v; y_v = x0*v+lambda*(x_v-x0*v), coefficient k).
    !   Terms containing b_n vanish because b_n = 0.
    ! T_dd: the same recurrence after substituting the expansion above,
    !   d_{n+1} = ((y-b_n)*d_n+p_n'(x0)-c_n*d_{n-1})*(1/a_n), d_0 = d_1 = 0; no p_n(y)-p_n(x0) is formed,
    !   so it stays accurate as lambda -> 0.
    ! T(:,k,n) = 0 for k > n and T_dd(:,k,n) = 0 for k > n-2: only T(:,n+1,n) and T_dd(:,n-1,n), which the
    !   recurrence reads, are set (to 0); the rest is not referenced. T, T_dd explicit-shape: in step 3 the
    !   assumed-shape version cost about 1 us per target.
    integer(8), intent(in) :: deg
    real(8), intent(in) :: x0, lambda(:)
    real(8), intent(inout) :: T(size(lambda),0:deg,0:deg)
    real(8), intent(inout), optional :: T_dd(size(lambda),0:deg,0:deg)
    real(8) :: p(1,0:256), dp(1,0:256), x_v, y_v
    integer(8) :: n, n_prev, k, i, n_lambda

    n_lambda=size(lambda,kind=8)
    T(:,0:1,0)=0.0d0
    T(:,0,0)=1.0d0                 ! p_0 = 1
    ! p_{n+1} from p_n = T(:,:,n) (degree n) and p_{n-1} = T(:,:,n_prev); p_{-1} = 0, and c_0 = 0 drops the term.
    ! x_v: coefficient k of x*p_n, from x*p_j = a_j*p_{j+1}+c_j*p_{j-1}, j = k-1, k+1 (no j = -1).
    do n=0,deg-1
      n_prev=max(n-1,0_8)
      do i=1,n_lambda              ! k = 0
        x_v=T(i,1,n)*c_n(1)
        y_v=x0*T(i,0,n)+lambda(i)*(x_v-x0*T(i,0,n))
        T(i,0,n+1)=(y_v-c_n(n)*T(i,0,n_prev))*a_n_inv(n)
      enddo
      do k=1,n
        do i=1,n_lambda
          x_v=T(i,k-1,n)*a_n(k-1)+T(i,k+1,n)*c_n(k+1)
          y_v=x0*T(i,k,n)+lambda(i)*(x_v-x0*T(i,k,n))
          T(i,k,n+1)=(y_v-c_n(n)*T(i,k,n_prev))*a_n_inv(n)
        enddo
      enddo
      do i=1,n_lambda              ! k = n+1: x_v = a_n*(coefficient n), y_v = lambda*x_v
        T(i,n+1,n+1)=(lambda(i)*(T(i,n,n)*a_n(n)))*a_n_inv(n)
      enddo
      if (n+2<=deg) T(:,n+2,n+1)=0.0d0
    enddo
    if (.not.present(T_dd)) return
    call compute_legendre_values(deg,[x0],p=p,dp=dp)   ! dp(1,n) = p_n'(x0)
    T_dd(:,0:1,0:1)=0.0d0          ! d_0 = d_1 = 0
    ! d_{n+1} (degree n-1) from d_n (degree n-2) and d_{n-1}, as for T; p_n'(x0) enters coefficient 0.
    do n=1,deg-1
      do i=1,n_lambda              ! k = 0
        x_v=T_dd(i,1,n)*c_n(1)
        y_v=x0*T_dd(i,0,n)+lambda(i)*(x_v-x0*T_dd(i,0,n))
        T_dd(i,0,n+1)=((y_v-c_n(n)*T_dd(i,0,n-1))+dp(1,n))*a_n_inv(n)
      enddo
      do k=1,n-2
        do i=1,n_lambda
          x_v=T_dd(i,k-1,n)*a_n(k-1)+T_dd(i,k+1,n)*c_n(k+1)
          y_v=x0*T_dd(i,k,n)+lambda(i)*(x_v-x0*T_dd(i,k,n))
          T_dd(i,k,n+1)=(y_v-c_n(n)*T_dd(i,k,n-1))*a_n_inv(n)
        enddo
      enddo
      if (n>=2) then
        do i=1,n_lambda            ! k = n-1: x_v = a_{n-2}*(coefficient n-2), y_v = lambda*x_v
          T_dd(i,n-1,n+1)=(lambda(i)*(T_dd(i,n-2,n)*a_n(n-2)))*a_n_inv(n)
        enddo
      endif
      T_dd(:,n,n+1)=0.0d0
    enddo
  end subroutine compute_legendre_translation

  subroutine compute_source_geometry(C_Xi,src)
    ! Step 0: source-patch data, the same for all targets and kernels.
    ! C_Xi_fixed: C_Xi with the fixed-axis mode first and the moving-axis mode third;
    !   (:,:,:,1) for u fixed (edges 2, 4), (:,:,:,2) for v fixed (edges 1, 3).
    real(8), intent(in) :: C_Xi(16,16,3,2)
    type(source_geometry), intent(inout) :: src
    real(8) :: P_center(1,16), r_u_center(3), r_v_center(3), aspect
    integer(8) :: ifield, ic, ib, i

    do ifield=1,2
      do ic=1,3
        do ib=1,16
          do i=1,16
            src%C_Xi_fixed(i,ic+3*(ifield-1),ib,1)=C_Xi(i,ib,ic,ifield)
            src%C_Xi_fixed(i,ic+3*(ifield-1),ib,2)=C_Xi(ib,i,ic,ifield)
          enddo
        enddo
      enddo
    enddo
    call compute_legendre_values(15_8,[0.0d0],p=P_center)   ! Legendre values at 0
    do ic=1,3
      r_u_center(ic)=dot_product(matmul(P_center(1,:),C_Xi(:,:,ic,1)),P_center(1,:))
      r_v_center(ic)=dot_product(matmul(P_center(1,:),C_Xi(:,:,ic,2)),P_center(1,:))
    enddo
    aspect=norm2(r_v_center)/norm2(r_u_center)
    src%if_har=max(aspect,1.0d0/aspect)>8.0d0
  end subroutine compute_source_geometry

  subroutine compute_target_geometry(t0,h,C_Xi,C_Xij,geo)
    ! Step 1: target geometry, the same for all four kernels.
    ! C_Xij: second derivatives, used off surface only.
    ! Off surface (h/=0) also r_uu0, r_uv0, r_vv0 and K, the order of step 8's expansion in d_curv:
    !   h_eff = |h|*(largest |principal curvature|), from the shape operator in an orthonormal tangent frame;
    !   K = 1, 2, 3, 4 for h_eff <= 1e-8, 1e-6, 1e-4, 1e-3, otherwise 5 (truncation ~ h_eff^(K+1)).
    real(8), intent(in) :: t0(2), h, C_Xi(16,16,3,2), C_Xij(16,16,3,3)
    type(target_geometry), intent(inout) :: geo
    real(8) :: P0(2,16), dP0(2,16), A0(16,6), A_ab(16,9), n_r0(3), tangent_norm(2), tangent_unit(3,2), tangent_cos
    real(8) :: tangent_cross(3), tangent_sin, curvature_h(3), shape_uv, shape_vv, curvature_mid, h_eff
    integer(8) :: ic
    external :: dgemv

    geo%t0=t0
    geo%h=h
    call compute_legendre_values(15_8,t0,p=P0,dp=dP0)
    geo%P0=transpose(P0)              ! (16,2): u' in column 1, v' in column 2
    geo%dP0=transpose(dP0)
    ! u first (A0 = C_Xi^T P0(:,1), C_Xi as 16 x 96), then v.
    call dgemv('T',16_8,96_8,1.0d0,C_Xi,16_8,geo%P0(1,1),1_8,0.0d0,A0,1_8)
    do ic=1,3
      geo%r_u0(ic)=dot_product(A0(:,ic),geo%P0(:,2))
      geo%r_v0(ic)=dot_product(A0(:,ic+3),geo%P0(:,2))
    enddo
    geo%metric(1,1)=dot_product(geo%r_u0,geo%r_u0)
    geo%metric(1,2)=dot_product(geo%r_u0,geo%r_v0)
    geo%metric(2,1)=geo%metric(1,2)
    geo%metric(2,2)=dot_product(geo%r_v0,geo%r_v0)
    geo%N0(1)=geo%r_u0(2)*geo%r_v0(3)-geo%r_u0(3)*geo%r_v0(2)
    geo%N0(2)=geo%r_u0(3)*geo%r_v0(1)-geo%r_u0(1)*geo%r_v0(3)
    geo%N0(3)=geo%r_u0(1)*geo%r_v0(2)-geo%r_u0(2)*geo%r_v0(1)
    geo%J0=sqrt(dot_product(geo%N0,geo%N0))   ! from N0 directly, not sqrt(g11*g22-g12^2)
    geo%x_E=[1.0d0+t0(2),1.0d0-t0(1),1.0d0-t0(2),1.0d0+t0(1)]
    if (h/=0.0d0) then   ! off surface
      call dgemv('T',16_8,144_8,1.0d0,C_Xij,16_8,geo%P0(1,1),1_8,0.0d0,A_ab,1_8)
      do ic=1,3
        geo%r_uu0(ic)=dot_product(A_ab(:,ic),geo%P0(:,2))
        geo%r_uv0(ic)=dot_product(A_ab(:,ic+3),geo%P0(:,2))
        geo%r_vv0(ic)=dot_product(A_ab(:,ic+6),geo%P0(:,2))
      enddo
      n_r0=geo%N0/geo%J0
      ! h-scaled second fundamental form in the frame of the normalized tangent QR
      tangent_norm=[norm2(geo%r_u0),norm2(geo%r_v0)]
      tangent_unit(:,1)=geo%r_u0/tangent_norm(1)
      tangent_unit(:,2)=geo%r_v0/tangent_norm(2)
      tangent_cos=dot_product(tangent_unit(:,1),tangent_unit(:,2))
      tangent_cross(1)=tangent_unit(2,1)*tangent_unit(3,2)-tangent_unit(3,1)*tangent_unit(2,2)
      tangent_cross(2)=tangent_unit(3,1)*tangent_unit(1,2)-tangent_unit(1,1)*tangent_unit(3,2)
      tangent_cross(3)=tangent_unit(1,1)*tangent_unit(2,2)-tangent_unit(2,1)*tangent_unit(1,2)
      tangent_sin=norm2(tangent_cross)
      curvature_h(1)=(abs(h)/tangent_norm(1))*dot_product(n_r0,geo%r_uu0/tangent_norm(1))
      curvature_h(2)=(abs(h)/tangent_norm(1))*dot_product(n_r0,geo%r_uv0/tangent_norm(2))
      curvature_h(3)=(abs(h)/tangent_norm(2))*dot_product(n_r0,geo%r_vv0/tangent_norm(2))
      shape_uv=(curvature_h(2)-tangent_cos*curvature_h(1))/tangent_sin
      shape_vv=((curvature_h(3)-tangent_cos*curvature_h(2))/tangent_sin-tangent_cos*shape_uv)/tangent_sin
      curvature_mid=0.5d0*curvature_h(1)+0.5d0*shape_vv
      h_eff=abs(curvature_mid)+hypot(0.5d0*curvature_h(1)-0.5d0*shape_vv,shape_uv)
      geo%K=5_8
      if (h_eff<=1.0d-8) then
        geo%K=1_8
      elseif (h_eff<=1.0d-6) then
        geo%K=2_8
      elseif (h_eff<=1.0d-4) then
        geo%K=3_8
      elseif (h_eff<=1.0d-3) then
        geo%K=4_8
      endif
    endif
  end subroutine compute_target_geometry

  subroutine compute_edge_rule(iedge,geo,E,n_E,d,w_s)
    ! Step 2: sinh-mapped Gauss rule on one edge, the same for all four kernels.
    ! Edges, counterclockwise: 1 v=-1, 2 u=1, 3 v=1, 4 u=-1.
    ! On the edge, q0+h^2 = g_mm*[(s-s*)^2+delta^2]  (s: moving coordinate, s* = tm+shift).
    ! Nodes s = s* + sgn*delta*sinh(eta). Gauss order per panel from the Bernstein ellipse through eta=i*pi/2:
    ! n = smallest E%order >= max(8, B/log(rho)), B = E%budget, at most maxval(E%order).
    ! When s* lies on the edge, two panels split at eta=0 are used unless one panel needs fewer nodes
    ! (one panel is cheaper for a short eta-interval, the split for a long one; ties keep the split).
    ! A single panel whose B/log(rho) exceeds maxval(E%order) is not eligible: its capped rule misses the budget.
    ! d(:,1:n_E): displacement t-t0 at the nodes, moving component formed directly (HAR form).
    ! w_s(1:n_E): Gauss weight times ds/deta = delta*cosh(eta).
    ! A target on the edge (d_fixed=0) gives n_E=0: the edge factor (u-u') or (v-v') vanishes.
    ! d, w_s: allocated on first use to 2*maxval(E%order) nodes, then reused.
    integer(8), intent(in) :: iedge
    type(target_geometry), intent(in) :: geo
    type(edge_tables), intent(in) :: E
    integer(8), intent(inout) :: n_E
    real(8), allocatable, intent(inout) :: d(:,:), w_s(:)
    real(8), parameter :: edge_fixed(4)=[-1.0d0,1.0d0,1.0d0,-1.0d0]
    integer(8) :: fixed_axis, moving_axis, npanel, ipanel, n, irule(2), irule_split(2), i, j
    logical :: over_budget, half_over_budget
    real(8) :: d_fixed, g_mm, tm, sgn, shift, delta, eta_lo, eta_hi, panel(2,2)
    real(8) :: panel_half, eta, w_eta

    if (.not.allocated(d)) allocate(d(2,2*maxval(E%order)),w_s(2*maxval(E%order)))
    n_E=0_8
    fixed_axis=1_8+mod(iedge,2_8)
    moving_axis=3_8-fixed_axis
    d_fixed=edge_fixed(iedge)-geo%t0(fixed_axis)
    if (d_fixed==0.0d0) return
    g_mm=geo%metric(moving_axis,moving_axis)
    tm=geo%t0(moving_axis)
    sgn=merge(-1.0d0,1.0d0,iedge>2_8)
    shift=-geo%metric(1,2)*d_fixed/g_mm
    delta=abs(d_fixed)*geo%J0/g_mm
    delta=hypot(delta,abs(geo%h)/sqrt(g_mm))
    eta_lo=asinh(sgn*((-sgn-tm)-shift)/delta)
    eta_hi=asinh(sgn*((sgn-tm)-shift)/delta)
    npanel=1_8
    panel(1,:)=[eta_lo,eta_hi]
    call compute_panel_rule(eta_lo,eta_hi,irule(1),over_budget)
    if (eta_lo<0.0d0 .and. eta_hi>0.0d0) then   ! split at eta=0 unless one panel within the rule table needs fewer nodes
      call compute_panel_rule(eta_lo,0.0d0,irule_split(1),half_over_budget)
      call compute_panel_rule(0.0d0,eta_hi,irule_split(2),half_over_budget)
      if (over_budget .or. E%order(irule_split(1))+E%order(irule_split(2))<=E%order(irule(1))) then
        npanel=2_8
        panel(1,:)=[eta_lo,0.0d0]
        panel(2,:)=[0.0d0,eta_hi]
        irule=irule_split
      endif
    endif
    do ipanel=1,npanel
      panel_half=(panel(ipanel,2)-panel(ipanel,1))/2.0d0
      n=E%order(irule(ipanel))
      do i=1,n
        j=n_E+i
        eta=(E%nodes(i,irule(ipanel))+1.0d0)*panel_half+panel(ipanel,1)
        w_eta=E%weights(i,irule(ipanel))*panel_half
        w_s(j)=w_eta*delta*cosh(eta)
        d(fixed_axis,j)=d_fixed
        d(moving_axis,j)=shift+sgn*delta*sinh(eta)
      enddo
      n_E=n_E+n
    enddo
  contains
    subroutine compute_panel_rule(eta_a,eta_b,irule_panel,over_budget_panel)
      ! Rule index for the panel [eta_a,eta_b]: smallest E%order >= max(8, B/log(rho)), rho from the Bernstein
      ! ellipse through eta=i*pi/2, at most maxval(E%order); over_budget_panel: B/log(rho) > maxval(E%order).
      real(8), intent(in) :: eta_a, eta_b
      integer(8), intent(out) :: irule_panel
      logical, intent(out) :: over_budget_panel
      real(8) :: panel_mid, panel_half, n_need
      complex(8) :: z
      integer(8) :: n_p
      panel_mid=(eta_a+eta_b)/2.0d0
      panel_half=(eta_b-eta_a)/2.0d0
      z=cmplx(-panel_mid,pi/2.0d0,kind=8)/panel_half
      n_need=E%budget/abs(real(acosh(z),8))
      over_budget_panel=n_need>real(maxval(E%order),8)
      n_p=max(8_8,ceiling(min(n_need,real(maxval(E%order),8)),kind=8))
      irule_panel=minloc(E%order,1,mask=E%order>=n_p)
    end subroutine compute_panel_rule
  end subroutine compute_edge_rule

  subroutine compute_radial_expansion(order,iedge,geo,src,R,n_E,d,G,G0,G1,Gr,T_radial)
    ! Step 3: r_u, r_v on the rays of one edge, the same for all four kernels.
    ! Ray j, node i: t = t0 + lambda_i*d(:,j), lambda = R%lambda.
    ! G(i,1:6,j) = (r_u, r_v)(t) = G0 + lambda*G1(1,1:6,j) + lambda^2*Gr(i,1:6,j), G0(1,1:6) = (r_u0, r_v0).
    ! G1, Gr: formed only when order>=1 or src%if_har (not needed by SLP, DLP); otherwise left as they are.
    ! Legendre p_n = P_n, n = 0..15; the argument says the axis. Fixed axis (constant on the edge), y = t_fixed+lambda*d_fixed:
    !   p_n(y) = Pfixed = Pf0+lambda*Pf1+lambda^2*Pfr,  Pfr = d_fixed^2*p_n[t_fixed,t_fixed,y].
    ! Moving axis, y = tm+lambda*dm: p_n(y) = sum_k T(:,k,n)*Pedge(j,k) = Pm0+lambda*dm*Pm1+lambda^2*dm^2*(T_dd part),
    !   Pedge = p_k at the edge nodes tm+dm; T, T_dd: compute_legendre_translation centred at tm.
    !   Pfixed, Pedge, Pf1: compute_legendre_values; Pf0, Pm0, Pm1: geo%P0, geo%dP0.
    ! T_radial: tables, also used by step 9. The driver runs the edges as 1, 3, 2, 4; edges 1, 2 build what
    !   depends on t0 only (G0, T; Pf0, Pm0, Pm1, T_dd), and opposite edges 3, 4 reuse it.
    ! G, G1, Gr and the T_radial arrays: allocated on first use (size(d,2) nodes), then reused.
    integer(8), intent(in) :: order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(source_geometry), intent(in) :: src
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:)
    real(8), intent(inout) :: G0(1,6)
    real(8), allocatable, intent(inout) :: G(:,:,:), G1(:,:,:), Gr(:,:,:)
    type(radial_expansion_tables), intent(inout) :: T_radial
    real(8), parameter :: edge_fixed(4)=[-1.0d0,1.0d0,1.0d0,-1.0d0]
    integer(8) :: fixed_axis, moving_axis, ib, i, j, ell, ell2, n
    real(8) :: d_fixed, t_fixed, tm, dm, Cvalue, y_fixed(nr), P_t_fixed(1,16), dP_t_fixed(1,16)
    real(8) :: x_edge(512)                 ! tm+dm; n_E <= 2*256 (two panels of an edge_tables rule)
    real(8) :: A_ux, A_uy, A_uz, A_vx, A_vy, A_vz
    real(8) :: A(nr,6,16), Ar(nr,6,16)
    real(8) :: Cm0(16,6), Cm1(16,6), G1_m0(6), G1_m1(6), Gr_m0(nr,6), Gr_m1(nr,6)
    external :: dgemm, dgemv

    fixed_axis=1_8+mod(iedge,2_8)
    moving_axis=3_8-fixed_axis
    d_fixed=edge_fixed(iedge)-geo%t0(fixed_axis)
    t_fixed=geo%t0(fixed_axis)
    tm=geo%t0(moving_axis)
    if (.not.allocated(G)) allocate(G(nr,6,size(d,2)))
    if (.not.allocated(T_radial%T)) allocate(T_radial%T(nr,16,16),T_radial%Pfixed(nr,16), &
      T_radial%Pedge(size(d,2),16))
    if (iedge<=2_8) then
      G0(1,1:3)=geo%r_u0
      G0(1,4:6)=geo%r_v0
      if (order>=1_8 .or. src%if_har) then
        if (.not.allocated(T_radial%T_dd)) allocate(T_radial%T_dd(nr,16,16),T_radial%Pfr(nr,16),T_radial%Pf1r(nr,16))
        T_radial%Pf0=geo%P0(:,fixed_axis)
        T_radial%Pm0=geo%P0(:,moving_axis)
        T_radial%Pm1=geo%dP0(:,moving_axis)
        call compute_legendre_translation(15_8,tm,R%lambda,T_radial%T,T_radial%T_dd)
      else
        call compute_legendre_translation(15_8,tm,R%lambda,T_radial%T)
      endif
    endif
    ! Legendre values: Pfixed along the rays (fixed axis), Pedge at the edge nodes (moving axis).
    y_fixed=t_fixed+R%lambda*d_fixed
    call compute_legendre_values(15_8,y_fixed,p=T_radial%Pfixed)
    x_edge(1:n_E)=tm+d(moving_axis,1:n_E)
    call compute_legendre_values(15_8,x_edge(1:n_E),p=T_radial%Pedge)
    ! Fixed axis by DGEMM: A(i,:,b) = sum_a Pfixed(i,a)*C_Xi_fixed(a,:,b).
    call dgemm('N','N',nr,96_8,16_8,1.0d0,T_radial%Pfixed,nr,src%C_Xi_fixed(1,1,1,fixed_axis),16_8,0.0d0,A,nr)
    if (order>=1_8 .or. src%if_har) then
      if (.not.allocated(G1)) allocate(G1(1,6,size(d,2)),Gr(nr,6,size(d,2)))
      call compute_legendre_values(15_8,[t_fixed],d_fixed,P_t_fixed,dP_t_fixed)
      T_radial%Pf1=dP_t_fixed(1,:)   ! d_fixed*P'(t_fixed), the order-1 term along the ray
      ! Pfr = d_fixed^2*p_n[t_fixed,t_fixed,y], the one-point case of T_dd (compute_legendre_translation), DLMF 18.9.2:
      ! Pfr_{n+1} = ((y-b_n)*Pfr_n+d_fixed*Pf1_n-c_n*Pfr_{n-1})*(1/a_n), Pfr_0 = Pfr_1 = 0.
      T_radial%Pfr(:,1:2)=0.0d0
      do n=1,14
        T_radial%Pfr(:,n+2)=(((y_fixed-b_n(n))*T_radial%Pfr(:,n+1)+d_fixed*T_radial%Pf1(n+1)) &
          -c_n(n)*T_radial%Pfr(:,n))*a_n_inv(n)
      enddo
      do ib=1,16
        T_radial%Pf1r(:,ib)=T_radial%Pf1(ib)+R%lambda*T_radial%Pfr(:,ib)
      enddo
      call dgemv('N',96_8,16_8,1.0d0,src%C_Xi_fixed(1,1,1,fixed_axis),96_8,T_radial%Pm0,1_8,0.0d0,Cm0,1_8)
      call dgemv('N',96_8,16_8,1.0d0,src%C_Xi_fixed(1,1,1,fixed_axis),96_8,T_radial%Pm1,1_8,0.0d0,Cm1,1_8)
      call dgemv('T',16_8,6_8,1.0d0,Cm0,16_8,T_radial%Pf1,1_8,0.0d0,G1_m0,1_8)
      call dgemv('T',16_8,6_8,1.0d0,Cm1,16_8,T_radial%Pf0,1_8,0.0d0,G1_m1,1_8)
      ! The second remainder starts two moving modes above the output mode (T_dd(:,k,n) = 0 for k > n-2).
      !GCC$ unroll 14
      do ell=1,14
        !GCC$ ivdep
        do i=1,nr
          Cvalue=T_radial%T_dd(i,ell,ell+2)
          A_ux=Cvalue*A(i,1,ell+2)
          A_uy=Cvalue*A(i,2,ell+2)
          A_uz=Cvalue*A(i,3,ell+2)
          A_vx=Cvalue*A(i,4,ell+2)
          A_vy=Cvalue*A(i,5,ell+2)
          A_vz=Cvalue*A(i,6,ell+2)
          !GCC$ unroll 13
          do ell2=ell+3,16
            Cvalue=T_radial%T_dd(i,ell,ell2)
            A_ux=A_ux+Cvalue*A(i,1,ell2)
            A_uy=A_uy+Cvalue*A(i,2,ell2)
            A_uz=A_uz+Cvalue*A(i,3,ell2)
            A_vx=A_vx+Cvalue*A(i,4,ell2)
            A_vy=A_vy+Cvalue*A(i,5,ell2)
            A_vz=A_vz+Cvalue*A(i,6,ell2)
          enddo
          Ar(i,1,ell)=A_ux
          Ar(i,2,ell)=A_uy
          Ar(i,3,ell)=A_uz
          Ar(i,4,ell)=A_vx
          Ar(i,5,ell)=A_vy
          Ar(i,6,ell)=A_vz
        enddo
      enddo
      Ar(:,:,15:16)=0.0d0
      call dgemm('N','T',6*nr,n_E,16_8,1.0d0,Ar,6*nr,T_radial%Pedge,size(T_radial%Pedge,1,kind=8),0.0d0,Gr,6*nr)
      call dgemm('N','N',nr,6_8,16_8,1.0d0,T_radial%Pfr,nr,Cm0,16_8,0.0d0,Gr_m0,nr)
      call dgemm('N','N',nr,6_8,16_8,1.0d0,T_radial%Pf1r,nr,Cm1,16_8,0.0d0,Gr_m1,nr)
      do j=1,n_E
        dm=d(moving_axis,j)
        G1(1,:,j)=G1_m0+G1_m1*dm
        Gr(:,:,j)=Gr(:,:,j)*dm**2+Gr_m0+Gr_m1*dm
      enddo
    endif
    ! Moving axis by T, in place (ascending ell keeps the unused columns), then the edge nodes.
    !GCC$ unroll 16
    do ell=1,16
      !GCC$ ivdep
      do i=1,nr
        Cvalue=T_radial%T(i,ell,ell)
        A_ux=Cvalue*A(i,1,ell)
        A_uy=Cvalue*A(i,2,ell)
        A_uz=Cvalue*A(i,3,ell)
        A_vx=Cvalue*A(i,4,ell)
        A_vy=Cvalue*A(i,5,ell)
        A_vz=Cvalue*A(i,6,ell)
        !GCC$ unroll 15
        do ell2=ell+1,16
          Cvalue=T_radial%T(i,ell,ell2)
          A_ux=A_ux+Cvalue*A(i,1,ell2)
          A_uy=A_uy+Cvalue*A(i,2,ell2)
          A_uz=A_uz+Cvalue*A(i,3,ell2)
          A_vx=A_vx+Cvalue*A(i,4,ell2)
          A_vy=A_vy+Cvalue*A(i,5,ell2)
          A_vz=A_vz+Cvalue*A(i,6,ell2)
        enddo
        A(i,1,ell)=A_ux
        A(i,2,ell)=A_uy
        A(i,3,ell)=A_uz
        A(i,4,ell)=A_vx
        A(i,5,ell)=A_vy
        A(i,6,ell)=A_vz
      enddo
    enddo
    call dgemm('N','T',6*nr,n_E,16_8,1.0d0,A,6*nr,T_radial%Pedge,size(T_radial%Pedge,1,kind=8),0.0d0,G,6*nr)
  end subroutine compute_radial_expansion

  subroutine compute_displacement_expansion(order,src,R,n_E,d,G,G0,G1,Gr,X)
    ! Step 4: chord and tangent along the rays of one edge, the same for all four kernels.
    ! r(t0+lambda*d)-r(t0) = lambda*(C0+lambda*c1r), c1r = c1+lambda*cr;  dr/dlambda = V = C0+lambda*V1r.
    ! d.G = d_u*G(:,1:3)+d_v*G(:,4:6), likewise for G0, G1, Gr.
    ! Always: C0(3,n_E) = d.G0, q0(n_E) = |C0|^2, c1r(nr,3,n_E).
    ! Order 0, 1: V1r(nr,3,n_E), c1r = T1*V1r.
    !   (a) order 0, not HAR (SLP, DLP, SLPn off; G1, Gr not formed): V1r = (d.G-C0)/lambda.
    !   (b) order 1 or HAR: V1r = d.G1+lambda*d.Gr, no cancellation as lambda -> 0.
    ! Order 2: c1(1,3,n_E) = d.G1/2, cr(nr,3,n_E) = T2*(d.Gr), c1r = c1+lambda*cr.
    ! Fields not formed are left as they are. X arrays: allocated on first use (size(d,2) rays), then reused.
    integer(8), intent(in) :: order, n_E
    type(source_geometry), intent(in) :: src
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), G0(1,6)
    real(8), allocatable, intent(in) :: G(:,:,:), G1(:,:,:), Gr(:,:,:)
    type(displacement_expansion), intent(inout) :: X
    integer(8) :: i, j, ic
    real(8) :: d_u, d_v, C0_x, C0_y, C0_z
    real(8), allocatable :: Gr_d(:,:,:)
    external :: dgemm

    if (.not.allocated(X%C0)) allocate(X%C0(3,size(d,2)),X%q0(size(d,2)),X%c1r(nr,3,size(d,2)))
    do j=1,n_E
      X%C0(:,j)=d(1,j)*G0(1,1:3)+d(2,j)*G0(1,4:6)
      X%q0(j)=dot_product(X%C0(:,j),X%C0(:,j))
    enddo
    if (order==0_8 .and. .not.src%if_har) then
      if (.not.allocated(X%V1r)) allocate(X%V1r(nr,3,size(d,2)))
      do j=1,n_E
        d_u=d(1,j)
        d_v=d(2,j)
        C0_x=X%C0(1,j)
        C0_y=X%C0(2,j)
        C0_z=X%C0(3,j)
        !GCC$ ivdep
        do i=1,nr
          X%V1r(i,1,j)=(d_u*G(i,1,j)+d_v*G(i,4,j)-C0_x)*R%lambda_inv(i)
          X%V1r(i,2,j)=(d_u*G(i,2,j)+d_v*G(i,5,j)-C0_y)*R%lambda_inv(i)
          X%V1r(i,3,j)=(d_u*G(i,3,j)+d_v*G(i,6,j)-C0_z)*R%lambda_inv(i)
        enddo
      enddo
      call dgemm('N','N',nr,3*n_E,nr,1.0d0,R%T1,nr,X%V1r,nr,0.0d0,X%c1r,nr)
    elseif (order<=1_8) then
      if (.not.allocated(X%V1r)) allocate(X%V1r(nr,3,size(d,2)))
      do j=1,n_E
        d_u=d(1,j)
        d_v=d(2,j)
        do ic=1,3
          !GCC$ ivdep
          do i=1,nr
            X%V1r(i,ic,j)=d_u*(G1(1,ic,j)+R%lambda(i)*Gr(i,ic,j))+d_v*(G1(1,ic+3,j)+R%lambda(i)*Gr(i,ic+3,j))
          enddo
        enddo
      enddo
      call dgemm('N','N',nr,3*n_E,nr,1.0d0,R%T1,nr,X%V1r,nr,0.0d0,X%c1r,nr)
    else
      allocate(Gr_d(nr,3,n_E))
      if (.not.allocated(X%c1)) allocate(X%c1(1,3,size(d,2)),X%cr(nr,3,size(d,2)))
      do j=1,n_E
        d_u=d(1,j)
        d_v=d(2,j)
        X%c1(1,:,j)=(d_u*G1(1,1:3,j)+d_v*G1(1,4:6,j))/2.0d0
        do ic=1,3
          !GCC$ ivdep
          do i=1,nr
            Gr_d(i,ic,j)=d_u*Gr(i,ic,j)+d_v*Gr(i,ic+3,j)
          enddo
        enddo
      enddo
      call dgemm('N','N',nr,3*n_E,nr,1.0d0,R%T2,nr,Gr_d,nr,0.0d0,X%cr,nr)
      do j=1,n_E
        do ic=1,3
          X%c1r(:,ic,j)=X%c1(1,ic,j)+R%lambda*X%cr(:,ic,j)
        enddo
      enddo
    endif
  end subroutine compute_displacement_expansion

  subroutine compute_jacobian_expansion(order,if_N,R,n_E,G,G0,G1,Gr,Jx)
    ! Step 5: N = r_u x r_v and J = |N| on the rays of one edge, the same for all four kernels.
    ! N = N0 + lambda*N1 + lambda^2*Nr,  J = J0 + lambda*J1 + lambda^2*Jr = J0 + lambda*J1r,
    ! J^2 = J0^2 + lambda*Q1 + lambda^2*Qr:  J1 = Q1/(2*J0), J1r = (Q1+lambda*Qr)/(J+J0), Jr = (Qr-J1r^2)/(2*J0).
    ! Always: N0(1,3), J0 from G0. From G directly: N(nr,3,n_E) if if_N (DLP), otherwise J(nr,n_E);
    !   DLP reads only N, the other kernels only J.
    ! Order >= 1 (never with if_N): N1(1,3,n_E), J1(n_E), Nr(nr,3,n_E), J1r(nr,n_E), Jr(nr,n_E) from G0, G1, Gr,
    !   with G1r = G1+lambda*Gr.
    ! Fields not formed are left as they are. Jx arrays: allocated on first use (size(G,3) rays), then reused.
    integer(8), intent(in) :: order, n_E
    logical, intent(in) :: if_N
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: G0(1,6)
    real(8), allocatable, intent(in) :: G(:,:,:), G1(:,:,:), Gr(:,:,:)
    type(jacobian_expansion), intent(inout) :: Jx
    integer(8) :: i, j
    real(8) :: N_x, N_y, N_z, N1_x, N1_y, N1_z, Nr_x, Nr_y, Nr_z, Q1, Qr, J1r
    real(8) :: G1r_ux, G1r_uy, G1r_uz, G1r_vx, G1r_vy, G1r_vz

    Jx%N0(1,1)=G0(1,2)*G0(1,6)-G0(1,3)*G0(1,5)
    Jx%N0(1,2)=G0(1,3)*G0(1,4)-G0(1,1)*G0(1,6)
    Jx%N0(1,3)=G0(1,1)*G0(1,5)-G0(1,2)*G0(1,4)
    Jx%J0=sqrt((Jx%N0(1,1)*Jx%N0(1,1)+Jx%N0(1,2)*Jx%N0(1,2))+Jx%N0(1,3)*Jx%N0(1,3))
    if (if_N) then
      if (.not.allocated(Jx%N)) allocate(Jx%N(nr,3,size(G,3)))
      do j=1,n_E
        !GCC$ ivdep
        do i=1,nr
          Jx%N(i,1,j)=G(i,2,j)*G(i,6,j)-G(i,3,j)*G(i,5,j)
          Jx%N(i,2,j)=G(i,3,j)*G(i,4,j)-G(i,1,j)*G(i,6,j)
          Jx%N(i,3,j)=G(i,1,j)*G(i,5,j)-G(i,2,j)*G(i,4,j)
        enddo
      enddo
    else
      if (.not.allocated(Jx%J)) allocate(Jx%J(nr,size(G,3)))
      do j=1,n_E
        !GCC$ ivdep
        do i=1,nr
          N_x=G(i,2,j)*G(i,6,j)-G(i,3,j)*G(i,5,j)
          N_y=G(i,3,j)*G(i,4,j)-G(i,1,j)*G(i,6,j)
          N_z=G(i,1,j)*G(i,5,j)-G(i,2,j)*G(i,4,j)
          Jx%J(i,j)=sqrt((N_x*N_x+N_y*N_y)+N_z*N_z)
        enddo
      enddo
    endif
    if (order>=1_8) then
      if (.not.allocated(Jx%N1)) allocate(Jx%N1(1,3,size(G,3)),Jx%J1(size(G,3)),Jx%Nr(nr,3,size(G,3)), &
        Jx%J1r(nr,size(G,3)),Jx%Jr(nr,size(G,3)))
      do j=1,n_E
        N1_x=(G1(1,2,j)*G0(1,6)-G1(1,3,j)*G0(1,5))+(G0(1,2)*G1(1,6,j)-G0(1,3)*G1(1,5,j))
        N1_y=(G1(1,3,j)*G0(1,4)-G1(1,1,j)*G0(1,6))+(G0(1,3)*G1(1,4,j)-G0(1,1)*G1(1,6,j))
        N1_z=(G1(1,1,j)*G0(1,5)-G1(1,2,j)*G0(1,4))+(G0(1,1)*G1(1,5,j)-G0(1,2)*G1(1,4,j))
        Jx%N1(1,1,j)=N1_x
        Jx%N1(1,2,j)=N1_y
        Jx%N1(1,3,j)=N1_z
        Q1=2.0d0*((Jx%N0(1,1)*N1_x+Jx%N0(1,2)*N1_y)+Jx%N0(1,3)*N1_z)
        Jx%J1(j)=Q1/(2.0d0*Jx%J0)
        !GCC$ ivdep
        do i=1,nr
          G1r_ux=G1(1,1,j)+R%lambda(i)*Gr(i,1,j)
          G1r_uy=G1(1,2,j)+R%lambda(i)*Gr(i,2,j)
          G1r_uz=G1(1,3,j)+R%lambda(i)*Gr(i,3,j)
          G1r_vx=G1(1,4,j)+R%lambda(i)*Gr(i,4,j)
          G1r_vy=G1(1,5,j)+R%lambda(i)*Gr(i,5,j)
          G1r_vz=G1(1,6,j)+R%lambda(i)*Gr(i,6,j)
          Nr_x=((Gr(i,2,j)*G0(1,6)-Gr(i,3,j)*G0(1,5))+(G0(1,2)*Gr(i,6,j)-G0(1,3)*Gr(i,5,j)))+(G1r_uy*G1r_vz-G1r_uz*G1r_vy)
          Nr_y=((Gr(i,3,j)*G0(1,4)-Gr(i,1,j)*G0(1,6))+(G0(1,3)*Gr(i,4,j)-G0(1,1)*Gr(i,6,j)))+(G1r_uz*G1r_vx-G1r_ux*G1r_vz)
          Nr_z=((Gr(i,1,j)*G0(1,5)-Gr(i,2,j)*G0(1,4))+(G0(1,1)*Gr(i,5,j)-G0(1,2)*Gr(i,4,j)))+(G1r_ux*G1r_vy-G1r_uy*G1r_vx)
          Jx%Nr(i,1,j)=Nr_x
          Jx%Nr(i,2,j)=Nr_y
          Jx%Nr(i,3,j)=Nr_z
          Qr=2.0d0*((Jx%N0(1,1)*Nr_x+Jx%N0(1,2)*Nr_y)+Jx%N0(1,3)*Nr_z) &
            +(((N1_x+R%lambda(i)*Nr_x)**2+(N1_y+R%lambda(i)*Nr_y)**2)+(N1_z+R%lambda(i)*Nr_z)**2)
          J1r=(Q1+R%lambda(i)*Qr)/(Jx%J(i,j)+Jx%J0)
          Jx%J1r(i,j)=J1r
          Jx%Jr(i,j)=(Qr-J1r**2)/(2.0d0*Jx%J0)
        enddo
      enddo
    endif
  end subroutine compute_jacobian_expansion

  subroutine compute_distance_factor(p,order,R,geo,n_E,X,gp)
    ! Step 6: g_p = B^(-p-1/2) on the rays of one edge, one formula for all kernels and h.
    ! B = |C0+lambda*c1r-(h/lambda)*n_r0|^2 = q0*y^2, n_r0 = N0/J0, eps_sq = h^2/q0:
    !   y^2 = 1 + d_curv + eps_sq/lambda^2,  d_curv = (lambda^2*|c1r|^2 - 2*h*n_r0.c1r + 2*lambda*C0.c1r)/q0.
    ! Always: g0(n_E) = q0^(-p-1/2), g(nr,n_E) = g0/y^(2p+1). Per ray, B_ratio = y^2 first (h terms only
    !   off surface), then g for the given p. Off surface also d_curv(nr,n_E), used by step 8.
    ! Order 1 (p=1, SLPn on): g1r = (g-g0)/lambda = -g*beta1r*(y^2+y+1)/(y+1), beta1r = (lambda*|c1r|^2+2*C0.c1r)/q0.
    ! Order 2 (p=2 DLPn, p=3 Stokes DLPn, on surface): beta1 = 2*C0.c1/q0, betar = (2*C0.cr+|c1r|^2)/q0,
    !   g1 = -(p+1/2)*g0*beta1, gr = g0*(-(p+1/2)*betar+(beta1+lambda*betar)^2*r_y(y)), g1r = g1+lambda*gr.
    ! The order-2 split is on surface only: off-surface DLPn reads g (p=2) and d_curv. Fields not formed are left
    ! as they are. gp arrays: allocated on first use (size(X%q0) rays), then reused.
    integer(8), intent(in) :: p, order, n_E
    type(radial_tables), intent(in) :: R
    type(target_geometry), intent(in) :: geo
    type(displacement_expansion), intent(in) :: X
    type(distance_factor), intent(inout) :: gp
    integer(8) :: i, j
    real(8) :: n_r0(3), q0inv, q0inv_sqrt, eps_sq, c1r_x, c1r_y, c1r_z, c1r_sq, C0_c1r, n_c1r
    real(8) :: d_curv, B_ratio(nr), y, beta1r, betar, r_y, g0j

    if (.not.allocated(gp%g)) allocate(gp%g0(size(X%q0)),gp%g(nr,size(X%q0)))
    if (geo%h/=0.0d0 .and. .not.allocated(gp%d_curv)) allocate(gp%d_curv(nr,size(X%q0)))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      q0inv_sqrt=1.0d0/sqrt(X%q0(j))
      q0inv=1.0d0/X%q0(j)
      eps_sq=(abs(geo%h)/sqrt(X%q0(j)))**2
      g0j=q0inv_sqrt/(merge(X%q0(j),1.0d0,p>=1_8)*merge(X%q0(j),1.0d0,p>=2_8)*merge(X%q0(j),1.0d0,p>=3_8))   ! q0^p, p = 0..3
      gp%g0(j)=g0j
      if (geo%h==0.0d0) then   ! on surface: no h terms (1.3 us per target if kept)
        !GCC$ ivdep
        do i=1,nr
          c1r_x=X%c1r(i,1,j)
          c1r_y=X%c1r(i,2,j)
          c1r_z=X%c1r(i,3,j)
          c1r_sq=(c1r_x*c1r_x+c1r_y*c1r_y)+c1r_z*c1r_z
          C0_c1r=(X%C0(1,j)*c1r_x+X%C0(2,j)*c1r_y)+X%C0(3,j)*c1r_z
          B_ratio(i)=1.0d0+(R%lambda(i)**2*c1r_sq+2.0d0*R%lambda(i)*C0_c1r)*q0inv
        enddo
      else
        !GCC$ ivdep
        do i=1,nr
          c1r_x=X%c1r(i,1,j)
          c1r_y=X%c1r(i,2,j)
          c1r_z=X%c1r(i,3,j)
          c1r_sq=(c1r_x*c1r_x+c1r_y*c1r_y)+c1r_z*c1r_z
          C0_c1r=(X%C0(1,j)*c1r_x+X%C0(2,j)*c1r_y)+X%C0(3,j)*c1r_z
          n_c1r=(n_r0(1)*c1r_x+n_r0(2)*c1r_y)+n_r0(3)*c1r_z
          d_curv=(R%lambda(i)**2*c1r_sq-2.0d0*geo%h*n_c1r+2.0d0*R%lambda(i)*C0_c1r)*q0inv
          gp%d_curv(i,j)=d_curv
          B_ratio(i)=(1.0d0+d_curv)+eps_sq*R%lambda_inv(i)*R%lambda_inv(i)
        enddo
      endif
      ! g = g0/(B_ratio^p*sqrt(B_ratio)); p chosen outside the loop, which a runtime merge inside it would not allow.
      select case (p)
      case (0)
        gp%g(:,j)=g0j/sqrt(B_ratio)
      case (1)
        gp%g(:,j)=g0j/(B_ratio*sqrt(B_ratio))
      case (3)
        gp%g(:,j)=g0j/((B_ratio*B_ratio*B_ratio)*sqrt(B_ratio))
      case default
        gp%g(:,j)=g0j/((B_ratio*B_ratio)*sqrt(B_ratio))
      end select
    enddo
    if (order==1_8) then
      if (.not.allocated(gp%g1r)) allocate(gp%g1r(nr,size(X%q0)))
      do j=1,n_E
        q0inv=1.0d0/X%q0(j)
        eps_sq=(abs(geo%h)/sqrt(X%q0(j)))**2
        !GCC$ ivdep
        do i=1,nr
          c1r_x=X%c1r(i,1,j)
          c1r_y=X%c1r(i,2,j)
          c1r_z=X%c1r(i,3,j)
          c1r_sq=(c1r_x*c1r_x+c1r_y*c1r_y)+c1r_z*c1r_z
          C0_c1r=(X%C0(1,j)*c1r_x+X%C0(2,j)*c1r_y)+X%C0(3,j)*c1r_z
          n_c1r=(n_r0(1)*c1r_x+n_r0(2)*c1r_y)+n_r0(3)*c1r_z
          d_curv=(R%lambda(i)**2*c1r_sq-2.0d0*geo%h*n_c1r+2.0d0*R%lambda(i)*C0_c1r)*q0inv
          y=sqrt((1.0d0+d_curv)+eps_sq*R%lambda_inv(i)*R%lambda_inv(i))
          beta1r=(R%lambda(i)*c1r_sq+2.0d0*C0_c1r)*q0inv
          gp%g1r(i,j)=-gp%g(i,j)*beta1r*(y*y+y+1.0d0)/(y+1.0d0)
        enddo
      enddo
    elseif (order==2_8 .and. geo%h==0.0d0) then   ! the order-2 split is on surface only
      if (.not.allocated(gp%g1)) allocate(gp%g1(size(X%q0)),gp%gr(nr,size(X%q0)),gp%beta1(size(X%q0)), &
        gp%g1r(nr,size(X%q0)))
      do j=1,n_E
        q0inv=1.0d0/X%q0(j)
        eps_sq=(abs(geo%h)/sqrt(X%q0(j)))**2
        gp%beta1(j)=2.0d0*((X%C0(1,j)*X%c1(1,1,j)+X%C0(2,j)*X%c1(1,2,j))+X%C0(3,j)*X%c1(1,3,j))*q0inv
        gp%g1(j)=-(real(p,8)+0.5d0)*gp%g0(j)*gp%beta1(j)
        !GCC$ ivdep
        do i=1,nr
          c1r_x=X%c1r(i,1,j)
          c1r_y=X%c1r(i,2,j)
          c1r_z=X%c1r(i,3,j)
          c1r_sq=(c1r_x*c1r_x+c1r_y*c1r_y)+c1r_z*c1r_z
          C0_c1r=(X%C0(1,j)*c1r_x+X%C0(2,j)*c1r_y)+X%C0(3,j)*c1r_z
          n_c1r=(n_r0(1)*c1r_x+n_r0(2)*c1r_y)+n_r0(3)*c1r_z
          d_curv=(R%lambda(i)**2*c1r_sq-2.0d0*geo%h*n_c1r+2.0d0*R%lambda(i)*C0_c1r)*q0inv
          y=sqrt((1.0d0+d_curv)+eps_sq*R%lambda_inv(i)*R%lambda_inv(i))
          betar=(2.0d0*((X%C0(1,j)*X%cr(i,1,j)+X%C0(2,j)*X%cr(i,2,j))+X%C0(3,j)*X%cr(i,3,j))+c1r_sq)*q0inv
          if (p==3_8) then
            r_y=(((((((7.0d0*y+14.0d0)*y+12.0d0)*y+10.0d0)*y+8.0d0)*y+6.0d0)*y+4.0d0)*y+2.0d0)/(2.0d0*y**7*(y+1.0d0)**2)
          else
            r_y=(((((5.0d0*y+10.0d0)*y+8.0d0)*y+6.0d0)*y+4.0d0)*y+2.0d0)/(2.0d0*((y*y)*(y*y)*y)*((y+1.0d0)*(y+1.0d0)))
          endif
          gp%gr(i,j)=gp%g0(j)*(-(real(p,8)+0.5d0)*betar+(gp%beta1(j)+R%lambda(i)*betar)**2*r_y)
          gp%g1r(i,j)=gp%g1(j)+R%lambda(i)*gp%gr(i,j)
        enddo
      enddo
    endif
  end subroutine compute_distance_factor

  subroutine compute_kernelradial_S_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, SLP, on and off surface: lambda*F = J*g*mu. On surface it is smooth and step 8 adds nothing;
    ! off surface the h terms are in g (step 6) and compute_singular_correction_off corrects the Gauss sum.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,1) = f*w_R*J*g, the weights on mu(t'+lambda*d).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),1))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      !GCC$ ivdep
      do i=1,nr
        wk%weight_lambda(i,j,1)=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
      enddo
    enddo
  end subroutine compute_kernelradial_S_on

  subroutine compute_kernelradial_D_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, DLP on surface: r(t')-r(t'+lambda*d) = -lambda*V+lambda^2*H1, V.N = 0, so
    ! lambda*F = (H1.N)*g*mu is smooth (H1 = V1r-c1r) and step 8 adds nothing.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,1) = f*w_R*JP_D*g, JP_D = (V1r-c1r).N.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, JP_D

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),1))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      !GCC$ ivdep
      do i=1,nr
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        wk%weight_lambda(i,j,1)=factor*R%w(i)*(JP_D*gp%g(i,j))
      enddo
    enddo
  end subroutine compute_kernelradial_D_on

  subroutine compute_kernelradial_D_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, DLP off surface (h/=0): x = r(t')+h*n_r0, n_r0 = N0/J0, so
    ! r(t'+lambda*d)-x = lambda*V-lambda^2*H1-h*n_r0 with V.N = 0, and lambda*F = (JP_D+(h/lambda^2)*JP_Dh)*g*mu,
    ! JP_D = (V1r-c1r).N as on surface, JP_Dh = n_r0.N. The h terms of the distance are in g (step 6);
    ! compute_singular_correction_off corrects the Gauss sum near lambda = 0.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,1) = f*w_R*(JP_D+h*lambda^-2*JP_Dh)*g.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, JP_D, JP_Dh, n_r0(3)

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),1))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      !GCC$ ivdep
      do i=1,nr
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        JP_Dh=(n_r0(1)*Jx%N(i,1,j)+n_r0(2)*Jx%N(i,2,j))+n_r0(3)*Jx%N(i,3,j)
        wk%weight_lambda(i,j,1)=factor*R%w(i)*((JP_D+geo%h*JP_Dh*R%lambda_inv(i)*R%lambda_inv(i))*gp%g(i,j))
      enddo
    enddo
  end subroutine compute_kernelradial_D_off

  subroutine compute_kernelradial_Sn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, SLPn off surface (h/=0): x = r(t')+h*n_r0, n_r0 = N0/J0, so
    ! r(t'+lambda*d)-x = lambda^2*(C0/lambda+c1r-(h/lambda^2)*n_r0) and lambda*F = J*g*(C0/lambda+c1r-(h/lambda^2)*n_r0)*mu.
    ! mu is not split (order 0): h keeps C0/lambda finite; compute_singular_correction_off corrects the Gauss sum near 0.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,3) = f*w_R*J*g*(C0/lambda+c1r-h*lambda^-2*n_r0).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, w_Jg_Sn, n_r0(3)

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),3))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      !GCC$ ivdep
      do i=1,nr
        w_Jg_Sn=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        wk%weight_lambda(i,j,1)=w_Jg_Sn*((X%C0(1,j)*R%lambda_inv(i)+X%c1r(i,1,j)) &
          -geo%h*n_r0(1)*(R%lambda_inv(i)*R%lambda_inv(i)))
        wk%weight_lambda(i,j,2)=w_Jg_Sn*((X%C0(2,j)*R%lambda_inv(i)+X%c1r(i,2,j)) &
          -geo%h*n_r0(2)*(R%lambda_inv(i)*R%lambda_inv(i)))
        wk%weight_lambda(i,j,3)=w_Jg_Sn*((X%C0(3,j)*R%lambda_inv(i)+X%c1r(i,3,j)) &
          -geo%h*n_r0(3)*(R%lambda_inv(i)*R%lambda_inv(i)))
      enddo
    enddo
  end subroutine compute_kernelradial_Sn_off

  subroutine compute_kernelradial_Dn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, DLPn off surface (h/=0): the Hessian of the single layer (the target gradient of SLPn), all 9
    ! components, (a,k) stored at a+3*(k-1). x = r(t')+h*n_r0, n_r0 = N0/J0:
    !   r(t'+lambda*d)-x = lambda*c,  c = C0+lambda*c1r-(h/lambda)*n_r0,  |c|^2 = q0*y^2,
    !   lambda*F_ak = J*g*(3*c_a*c_k-|c|^2*delta_ak)/lambda^2,  g = q0^(-5/2)*y^(-5) (step 6, p=2).
    ! With mu = mu0+lambda*mu1+lambda^2*mur, f = x_E*w_s/(2*pi):
    !   weight_r(nr,n_E,9) = f*w_R*J*g*(3*c*c'-|c|^2*I)   (on mur).
    ! The endpoint terms on mu0, mu1 are formed for every ray in compute_singular_correction_off.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, ia, ik
    real(8) :: factor, n_r0(3), c(3), c_sq, w_Jg

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),9))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        c(1)=(X%C0(1,j)+R%lambda(i)*X%c1r(i,1,j))-(geo%h*R%lambda_inv(i))*n_r0(1)
        c(2)=(X%C0(2,j)+R%lambda(i)*X%c1r(i,2,j))-(geo%h*R%lambda_inv(i))*n_r0(2)
        c(3)=(X%C0(3,j)+R%lambda(i)*X%c1r(i,3,j))-(geo%h*R%lambda_inv(i))*n_r0(3)
        c_sq=(c(1)*c(1)+c(2)*c(2))+c(3)*c(3)
        w_Jg=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        do ik=1,3
          do ia=1,3
            if (ia==ik) then
              wk%weight_r(i,j,ia+3*(ik-1))=w_Jg*(3.0d0*c(ia)*c(ik)-c_sq)
            else
              wk%weight_r(i,j,ia+3*(ik-1))=w_Jg*(3.0d0*c(ia)*c(ik))
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Dn_off

  subroutine compute_kernelradial_Sn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, SLPn on surface, the Gauss sums. With mu = mu0+lambda*mu1r and J*g = J0*g0+lambda*Jg_Sn1r:
    !   lambda*F = (C0/lambda+c1r)*J*g*mu = c1r*J*g*mu + C0*J*g*mu1r + C0*Jg_Sn1r*mu0 + C0*J0*g0*mu0/lambda,
    ! the last term being step 8. f = x_E*w_s/(2*pi):
    !   weight_lambda(nr,n_E,3) = f*w_R*J*g*c1r (on mu),  weight_r(nr,n_E,3) = f*w_R*J*g*C0 (on mu1r),
    !   weight_0(3) += sum_j f*(sum_i w_R*Jg_Sn1r)*C0 (on mu(t')).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, w_Jg_Sn, Jg_Sn1r(nr), s

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),3),wk%weight_r(nr,size(d,2),3))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      !GCC$ ivdep
      do i=1,nr
        w_Jg_Sn=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        wk%weight_lambda(i,j,1)=w_Jg_Sn*X%c1r(i,1,j)
        wk%weight_lambda(i,j,2)=w_Jg_Sn*X%c1r(i,2,j)
        wk%weight_lambda(i,j,3)=w_Jg_Sn*X%c1r(i,3,j)
        wk%weight_r(i,j,1)=w_Jg_Sn*X%C0(1,j)
        wk%weight_r(i,j,2)=w_Jg_Sn*X%C0(2,j)
        wk%weight_r(i,j,3)=w_Jg_Sn*X%C0(3,j)
        Jg_Sn1r(i)=Jx%J1r(i,j)*gp%g(i,j)+Jx%J0*gp%g1r(i,j)
      enddo
      s=0.0d0
      do i=1,nr
        s=s+R%w(i)*Jg_Sn1r(i)
      enddo
      weight_0=weight_0+(factor*s)*X%C0(:,j)
    enddo
  end subroutine compute_kernelradial_Sn_on

  subroutine compute_kernelradial_Dn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, DLPn on surface, the Gauss sums. K = J*P*g = K0+lambda*K1+lambda^2*Kr, P = 3*c*c'-|c|^2*I
    ! (upper triangle (1,1),(1,2),(1,3),(2,2),(2,3),(3,3)), mu = mu0+lambda*mu1+lambda^2*mur:
    !   lambda*F = K*mu/lambda^2 = K*mur + Kr*mu0 + lambda*Kr*mu1 + [(K0/lambda^2+K1/lambda)*mu0+(K0/lambda+K1)*mu1],
    ! the bracket being step 8; mu1 = d.grad(mu(t')). f = x_E*w_s/(2*pi):
    !   weight_r(nr,n_E,6) = f*w_R*K (on mur),
    !   weight_0 += sum_j sum_i f*w_R*Kr,  (weight_u, weight_v) += sum_j (sum_i f*w_R*lambda*Kr)*d(:,j).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, C0(3), c1(3), B0, B1, Br(nr), dij, P_Dn0, P_Dn1, P_Dnr, JP_Dn0, JP_Dn1, JP_Dnr, JP_Dn
    real(8) :: K_Dn0, K_Dn1, K_Dnr, regular_weight(nr), regular_0, regular_1

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),6))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      C0=X%C0(:,j)
      c1=X%c1(1,:,j)
      B0=X%q0(j)
      B1=2.0d0*((C0(1)*c1(1)+C0(2)*c1(2))+C0(3)*c1(3))
      !GCC$ ivdep
      do i=1,nr
        Br(i)=2.0d0*((C0(1)*X%cr(i,1,j)+C0(2)*X%cr(i,2,j))+C0(3)*X%cr(i,3,j)) &
          +((X%c1r(i,1,j)*X%c1r(i,1,j)+X%c1r(i,2,j)*X%c1r(i,2,j))+X%c1r(i,3,j)*X%c1r(i,3,j))
      enddo
      k=0
      do ic=1,3
        do jc=ic,3
          k=k+1
          dij=merge(1.0d0,0.0d0,ic==jc)
          P_Dn0=3.0d0*C0(ic)*C0(jc)-dij*B0
          P_Dn1=3.0d0*(c1(ic)*C0(jc)+C0(ic)*c1(jc))-dij*B1
          JP_Dn0=Jx%J0*P_Dn0
          JP_Dn1=Jx%J1(j)*P_Dn0+Jx%J0*P_Dn1
          K_Dn0=JP_Dn0*gp%g0(j)
          K_Dn1=JP_Dn1*gp%g0(j)+JP_Dn0*gp%g1(j)
          !GCC$ ivdep
          do i=1,nr
            P_Dnr=3.0d0*((X%cr(i,ic,j)*(C0(jc)+R%lambda(i)*c1(jc))+c1(ic)*c1(jc)) &
              +(C0(ic)+R%lambda(i)*X%c1r(i,ic,j))*X%cr(i,jc,j))-dij*Br(i)
            JP_Dnr=Jx%Jr(i,j)*(P_Dn0+R%lambda(i)*P_Dn1)+Jx%J1(j)*P_Dn1+Jx%J(i,j)*P_Dnr
            JP_Dn=Jx%J(i,j)*(P_Dn0+R%lambda(i)*(P_Dn1+R%lambda(i)*P_Dnr))
            K_Dnr=JP_Dnr*gp%g0(j)+(JP_Dn1+R%lambda(i)*JP_Dnr)*gp%g1(j)+JP_Dn*gp%gr(i,j)
            regular_weight(i)=factor*R%w(i)*K_Dnr
            wk%weight_r(i,j,k)=factor*R%w(i)*(K_Dn0+R%lambda(i)*K_Dn1)+R%lambda(i)**2*regular_weight(i)
          enddo
          regular_0=0.0d0
          regular_1=0.0d0
          do i=1,nr
            regular_0=regular_0+regular_weight(i)
            regular_1=regular_1+R%lambda(i)*regular_weight(i)
          enddo
          weight_0(k)=weight_0(k)+regular_0
          weight_u(k)=weight_u(k)+regular_1*d(1,j)
          weight_v(k)=weight_v(k)+regular_1*d(2,j)
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Dn_on

  subroutine compute_kernelradial_S_helm_on(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz SLP on surface: compute_kernelradial_S_on times e = exp(1i*z), z = zk*lambda*|c|,
    ! c = C0+lambda*c1r (any complex zk); smooth, step 8 adds nothing. weight_lambda(nr,n_E,2): real part, imaginary part.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, w_Jg, c_abs
    complex(8) :: e_z

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),2))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        c_abs=sqrt((X%C0(1,j)+R%lambda(i)*X%c1r(i,1,j))**2+(X%C0(2,j)+R%lambda(i)*X%c1r(i,2,j))**2 &
          +(X%C0(3,j)+R%lambda(i)*X%c1r(i,3,j))**2)
        e_z=exp((0.0d0,1.0d0)*(zk*(R%lambda(i)*c_abs)))
        w_Jg=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        wk%weight_lambda(i,j,1)=w_Jg*real(e_z,8)
        wk%weight_lambda(i,j,2)=w_Jg*aimag(e_z)
      enddo
    enddo
  end subroutine compute_kernelradial_S_helm_on

  subroutine compute_kernelradial_D_helm_on(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz DLP on surface: compute_kernelradial_D_on times c = e*(1-1i*z), e = exp(1i*z), z = zk*lambda*|c|
    ! (any complex zk); smooth, step 8 adds nothing. weight_lambda(nr,n_E,2): real part, imaginary part.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, w_D, JP_D, c_abs
    complex(8) :: z, c_z

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),2))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        c_abs=sqrt((X%C0(1,j)+R%lambda(i)*X%c1r(i,1,j))**2+(X%C0(2,j)+R%lambda(i)*X%c1r(i,2,j))**2 &
          +(X%C0(3,j)+R%lambda(i)*X%c1r(i,3,j))**2)
        z=zk*(R%lambda(i)*c_abs)
        c_z=exp((0.0d0,1.0d0)*z)*(1.0d0-(0.0d0,1.0d0)*z)
        w_D=factor*R%w(i)*(JP_D*gp%g(i,j))
        wk%weight_lambda(i,j,1)=w_D*real(c_z,8)
        wk%weight_lambda(i,j,2)=w_D*aimag(c_z)
      enddo
    enddo
  end subroutine compute_kernelradial_D_helm_on

  subroutine compute_kernelradial_Sn_helm_on(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz SLPn on surface: compute_kernelradial_Sn_on with J*g replaced by J*g*c_z, c_z = e*(1-1i*z),
    ! e = exp(1i*z), z = zk*lambda*|c| (any complex zk). c_z = 1+O(lambda^2): the lambda^-1 term and step 8 are
    ! Laplace's (driver: step 8 on components 1:3), and (J*g*c_z)_1r = (J*g)_1r*c_z+J0*g0*zk^2*lambda*|c|^2*phi.
    ! Components 1:3 real parts, 4:6 imaginary parts of weight_lambda, weight_r, weight_0.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, n_term
    real(8) :: factor, c(3), c_sq
    complex(8) :: z, c_z, phi, term, w_Jg, Jg1r, s

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),6),wk%weight_r(nr,size(d,2),6))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      s=(0.0d0,0.0d0)
      do i=1,nr
        c=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)
        c_sq=(c(1)*c(1)+c(2)*c(2))+c(3)*c(3)
        z=zk*(R%lambda(i)*sqrt(c_sq))
        c_z=exp((0.0d0,1.0d0)*z)*(1.0d0-(0.0d0,1.0d0)*z)
        ! phi = (c_z-1)/z^2 = sum_{n>=2} (n-1)*(1i*z)^(n-2)/n!: the series for |z| < 1
        if (abs(z)<1.0d0) then
          term=(0.5d0,0.0d0); phi=term
          do n_term=1,20
            term=term*(0.0d0,1.0d0)*z/real(n_term+2,8)
            phi=phi+real(n_term+1,8)*term
          enddo
        else
          phi=(c_z-1.0d0)/z**2
        endif
        w_Jg=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))*c_z
        wk%weight_lambda(i,j,1:3)=real(w_Jg,8)*X%c1r(i,:,j); wk%weight_lambda(i,j,4:6)=aimag(w_Jg)*X%c1r(i,:,j)
        wk%weight_r(i,j,1:3)=real(w_Jg,8)*X%C0(:,j); wk%weight_r(i,j,4:6)=aimag(w_Jg)*X%C0(:,j)
        Jg1r=(Jx%J1r(i,j)*gp%g(i,j)+Jx%J0*gp%g1r(i,j))*c_z+Jx%J0*gp%g0(j)*(zk**2*(R%lambda(i)*c_sq))*phi
        s=s+R%w(i)*Jg1r
      enddo
      weight_0(1:3)=weight_0(1:3)+(factor*real(s,8))*X%C0(:,j)
      weight_0(4:6)=weight_0(4:6)+(factor*aimag(s))*X%C0(:,j)
    enddo
  end subroutine compute_kernelradial_Sn_helm_on

  subroutine compute_kernelradial_Dn_helm_on(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz DLPn on surface: the target Hessian of exp(1i*z)/|r| is c_z*(3*r*r'-|r|^2*I)/|r|^5
    ! -zk^2*e*r*r'/|r|^3, e = exp(1i*z), c_z = e*(1-1i*z), z = zk*lambda*|c|, c = C0+lambda*c1r (any complex zk).
    ! As compute_kernelradial_Dn_on with K = c_z*J*P*g+lambda^2*E = K0+lambda*K1+lambda^2*Kr, E = -zk^2*e*J*c*c'*|c|^2*g:
    ! c_z = 1+O(lambda^2), so K0, K1 and step 8 are Laplace's (driver: step 8 on components 1:6), and
    !   Kr = c_z*Kr_Laplace+zk^2*|c|^2*phi*(K0+lambda*K1)+E,  phi = (c_z-1)/z^2.
    ! Components 1:6 real parts, 7:12 imaginary parts of weight_r, weight_0, weight_u, weight_v.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc, n_term
    real(8) :: factor, C0(3), c1(3), B0, B1, Br(nr), dij, P_Dn0, P_Dn1, P_Dnr, JP_Dn0, JP_Dn1, JP_Dnr, JP_Dn
    real(8) :: K_Dn0, K_Dn1, K_Dnr, c(nr,3), c_sq(nr)
    complex(8) :: z, e_z(nr), c_z(nr), phi(nr), term, regular_weight(nr), regular_0, regular_1

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),12))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      C0=X%C0(:,j)
      c1=X%c1(1,:,j)
      B0=X%q0(j)
      B1=2.0d0*((C0(1)*c1(1)+C0(2)*c1(2))+C0(3)*c1(3))
      do i=1,nr
        Br(i)=2.0d0*((C0(1)*X%cr(i,1,j)+C0(2)*X%cr(i,2,j))+C0(3)*X%cr(i,3,j)) &
          +((X%c1r(i,1,j)*X%c1r(i,1,j)+X%c1r(i,2,j)*X%c1r(i,2,j))+X%c1r(i,3,j)*X%c1r(i,3,j))
        c(i,:)=C0+R%lambda(i)*X%c1r(i,:,j)
        c_sq(i)=B0+R%lambda(i)*(B1+R%lambda(i)*Br(i))
        z=zk*(R%lambda(i)*sqrt(c_sq(i)))
        e_z(i)=exp((0.0d0,1.0d0)*z); c_z(i)=e_z(i)*(1.0d0-(0.0d0,1.0d0)*z)
        ! phi(i) = (c_z(i)-1)/z^2 = sum_{n>=2} (n-1)*(1i*z)^(n-2)/n!: the series for |z| < 1
        if (abs(z)<1.0d0) then
          term=(0.5d0,0.0d0); phi(i)=term
          do n_term=1,20
            term=term*(0.0d0,1.0d0)*z/real(n_term+2,8)
            phi(i)=phi(i)+real(n_term+1,8)*term
          enddo
        else
          phi(i)=(c_z(i)-1.0d0)/z**2
        endif
      enddo
      k=0
      do ic=1,3
        do jc=ic,3
          k=k+1
          dij=merge(1.0d0,0.0d0,ic==jc)
          P_Dn0=3.0d0*C0(ic)*C0(jc)-dij*B0
          P_Dn1=3.0d0*(c1(ic)*C0(jc)+C0(ic)*c1(jc))-dij*B1
          JP_Dn0=Jx%J0*P_Dn0
          JP_Dn1=Jx%J1(j)*P_Dn0+Jx%J0*P_Dn1
          K_Dn0=JP_Dn0*gp%g0(j)
          K_Dn1=JP_Dn1*gp%g0(j)+JP_Dn0*gp%g1(j)
          do i=1,nr
            P_Dnr=3.0d0*((X%cr(i,ic,j)*(C0(jc)+R%lambda(i)*c1(jc))+c1(ic)*c1(jc)) &
              +(C0(ic)+R%lambda(i)*X%c1r(i,ic,j))*X%cr(i,jc,j))-dij*Br(i)
            JP_Dnr=Jx%Jr(i,j)*(P_Dn0+R%lambda(i)*P_Dn1)+Jx%J1(j)*P_Dn1+Jx%J(i,j)*P_Dnr
            JP_Dn=Jx%J(i,j)*(P_Dn0+R%lambda(i)*(P_Dn1+R%lambda(i)*P_Dnr))
            K_Dnr=JP_Dnr*gp%g0(j)+(JP_Dn1+R%lambda(i)*JP_Dnr)*gp%g1(j)+JP_Dn*gp%gr(i,j)
            regular_weight(i)=factor*R%w(i)*(c_z(i)*K_Dnr+(zk**2*(c_sq(i)*phi(i)))*(K_Dn0+R%lambda(i)*K_Dn1) &
              -zk**2*e_z(i)*(Jx%J(i,j)*c(i,ic)*c(i,jc)*c_sq(i)*gp%g(i,j)))
            wk%weight_r(i,j,k)=factor*R%w(i)*(K_Dn0+R%lambda(i)*K_Dn1)+R%lambda(i)**2*real(regular_weight(i),8)
            wk%weight_r(i,j,6+k)=R%lambda(i)**2*aimag(regular_weight(i))
          enddo
          regular_0=(0.0d0,0.0d0)
          regular_1=(0.0d0,0.0d0)
          do i=1,nr
            regular_0=regular_0+regular_weight(i)
            regular_1=regular_1+R%lambda(i)*regular_weight(i)
          enddo
          weight_0(k)=weight_0(k)+real(regular_0,8); weight_0(6+k)=weight_0(6+k)+aimag(regular_0)
          weight_u(k)=weight_u(k)+real(regular_1,8)*d(1,j); weight_u(6+k)=weight_u(6+k)+aimag(regular_1)*d(1,j)
          weight_v(k)=weight_v(k)+real(regular_1,8)*d(2,j); weight_v(6+k)=weight_v(6+k)+aimag(regular_1)*d(2,j)
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Dn_helm_on

  subroutine compute_kernelradial_S_helm_off(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz SLP off surface (h/=0), any complex zk. With the chord r(t'+lambda*d)-x = lambda*(C0+lambda*c1r)-h*n_r0
    ! of length dist and z = zk*dist, exp(1i*z)/dist = cos(z)/dist+1i*sin(z)/dist: cos(z) is a function of dist^2, smooth
    ! in lambda, and multiplies the Laplace SLP weight of compute_kernelradial_S_on (compute_singular_correction_helm_off
    ! corrects it with the same factor); sin(z)/dist is smooth and is summed by Gauss. Off surface dist is not smooth at
    ! the scale eps = |h|/sqrt(q0), so exp(1i*z) itself cannot multiply the Laplace weight.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,2) = f*w_R*(J*g*cos(z)+1i*lambda*J*sin(z)/dist): real part, imaginary part.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j
    real(8) :: factor, n_r0(3), chord(3), dist
    complex(8) :: z, weight_S

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),2))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        chord=R%lambda(i)*(X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j))-geo%h*n_r0
        dist=sqrt((chord(1)*chord(1)+chord(2)*chord(2))+chord(3)*chord(3))
        z=zk*dist
        weight_S=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j)*cos(z)+(0.0d0,1.0d0)*R%lambda(i)*Jx%J(i,j)*sin(z)/dist)
        wk%weight_lambda(i,j,1)=real(weight_S,8)
        wk%weight_lambda(i,j,2)=aimag(weight_S)
      enddo
    enddo
  end subroutine compute_kernelradial_S_helm_off

  subroutine compute_kernelradial_D_helm_off(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz DLP off surface (h/=0), any complex zk. The Laplace DLP integrand of compute_kernelradial_D_off,
    ! lambda*F = (JP_D+(h/lambda^2)*JP_Dh)*g = lambda*T/dist^3, T = lambda^2*JP_D+h*JP_Dh, times
    ! c = exp(1i*z)*(1-1i*z) = c_e+1i*dist^3*q, z = zk*dist, dist = |lambda*(C0+lambda*c1r)-h*n_r0|:
    !   c_e = cos(z)+z*sin(z) (a function of dist^2, smooth in lambda) multiplies the Laplace weight, corrected with the
    !   same factor in compute_singular_correction_helm_off; q = zk^3*j1(z)/z (series for |z|<1) gives the smooth 1i*q*lambda*T.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,2) = f*w_R*(c_e*(JP_D+h*lambda^-2*JP_Dh)*g+1i*q*lambda*T): real, imaginary part.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, n_term
    real(8) :: factor, n_r0(3), chord(3), JP_D, JP_Dh
    complex(8) :: z, c_e, q, term, w_D

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),2))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        chord=R%lambda(i)*(X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j))-geo%h*n_r0
        z=zk*sqrt((chord(1)*chord(1)+chord(2)*chord(2))+chord(3)*chord(3))
        c_e=cos(z)+z*sin(z)
        ! q = zk^3*j1(z)/z = zk^3*sum_m (-z^2/2)^m/(m!*(2m+3)!!): the series for |z|<1
        if (abs(z)<1.0d0) then
          term=(1.0d0,0.0d0)/3.0d0; q=term
          do n_term=1,20
            term=term*(-z*z/2.0d0)/real(n_term*(2*n_term+3),8)
            q=q+term
          enddo
        else
          q=(sin(z)-z*cos(z))/z**3
        endif
        q=zk**3*q
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        JP_Dh=(n_r0(1)*Jx%N(i,1,j)+n_r0(2)*Jx%N(i,2,j))+n_r0(3)*Jx%N(i,3,j)
        w_D=factor*R%w(i)*(c_e*((JP_D+geo%h*JP_Dh*R%lambda_inv(i)*R%lambda_inv(i))*gp%g(i,j)) &
          +(0.0d0,1.0d0)*q*R%lambda(i)*(R%lambda(i)**2*JP_D+geo%h*JP_Dh))
        wk%weight_lambda(i,j,1)=real(w_D,8)
        wk%weight_lambda(i,j,2)=aimag(w_D)
      enddo
    enddo
  end subroutine compute_kernelradial_D_helm_off

  subroutine compute_kernelradial_Sn_helm_off(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz SLPn off surface (h/=0), any complex zk. The Laplace SLPn integrand of compute_kernelradial_Sn_off,
    ! lambda*F = J*g*(C0/lambda+c1r-(h/lambda^2)*n_r0) = lambda*J*chord/dist^3, chord = lambda*(C0+lambda*c1r)-h*n_r0,
    ! dist = |chord|, times c = exp(1i*z)*(1-1i*z) = c_e+1i*dist^3*q, z = zk*dist:
    !   c_e = cos(z)+z*sin(z) multiplies the Laplace weight (corrected with the same factor in
    !   compute_singular_correction_helm_off); q = zk^3*j1(z)/z (series for |z|<1) gives the smooth 1i*q*lambda*J*chord.
    ! f = x_E*w_s/(2*pi);  weight_lambda(nr,n_E,6) = f*w_R*(c_e*J*g*(C0/lambda+c1r-h*lambda^-2*n_r0)+1i*q*lambda*J*chord):
    ! components 1:3 real parts, 4:6 imaginary parts.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, ic, n_term
    real(8) :: factor, n_r0(3), chord(3), w_Jg_Sn
    complex(8) :: z, c_e, q, term, w_J_smooth, w_Sn

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),6))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        chord=R%lambda(i)*(X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j))-geo%h*n_r0
        z=zk*sqrt((chord(1)*chord(1)+chord(2)*chord(2))+chord(3)*chord(3))
        c_e=cos(z)+z*sin(z)
        ! q = zk^3*j1(z)/z = zk^3*sum_m (-z^2/2)^m/(m!*(2m+3)!!): the series for |z|<1
        if (abs(z)<1.0d0) then
          term=(1.0d0,0.0d0)/3.0d0; q=term
          do n_term=1,20
            term=term*(-z*z/2.0d0)/real(n_term*(2*n_term+3),8)
            q=q+term
          enddo
        else
          q=(sin(z)-z*cos(z))/z**3
        endif
        q=zk**3*q
        w_Jg_Sn=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        w_J_smooth=factor*R%w(i)*((0.0d0,1.0d0)*q*R%lambda(i)*Jx%J(i,j))
        do ic=1,3
          w_Sn=c_e*(w_Jg_Sn*((X%C0(ic,j)*R%lambda_inv(i)+X%c1r(i,ic,j)) &
            -geo%h*n_r0(ic)*(R%lambda_inv(i)*R%lambda_inv(i))))+w_J_smooth*chord(ic)
          wk%weight_lambda(i,j,ic)=real(w_Sn,8)
          wk%weight_lambda(i,j,3+ic)=aimag(w_Sn)
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Sn_helm_off

  subroutine compute_kernelradial_Dn_helm_off(zk,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Helmholtz DLPn off surface (h/=0): compute_kernelradial_Dn_off, the Laplace Hessian H on mur, in
    ! weight_r(:,:,1:9); compute_singular_correction_helm_off integrates H*psi, psi = A*mu, turns it into weights on mu,
    ! fills the imaginary parts weight_r(:,:,10:18) and adds the other Helmholtz terms. zk is not used here.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),18))
    call compute_kernelradial_Dn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
  end subroutine compute_kernelradial_Dn_helm_off

  subroutine compute_kernelradial_S_stok_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Stokes SLP on surface: (I/|r|+r*r'/|r|^3)/(8*pi), r = lambda*c, c = C0+lambda*c1r, so
    ! lambda*F = J*g*(I+c*c'/|c|^2)*mu with g = 1/|c| (p=0) is smooth and step 8 adds nothing.
    ! f = x_E*w_s/(8*pi);  weight_lambda(nr,n_E,6) = f*w_R*J*g*(delta_ij+c_i*c_j/|c|^2), the weights on mu_j(t'+lambda*d),
    ! (i,j) upper triangle (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, c(3), c_sq, w_Jg, dij

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),6))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(8.0d0*pi)
      do i=1,nr
        c=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)
        c_sq=(c(1)*c(1)+c(2)*c(2))+c(3)*c(3)
        w_Jg=factor*R%w(i)*(Jx%J(i,j)*gp%g(i,j))
        k=0
        do ic=1,3
          do jc=ic,3
            k=k+1
            dij=merge(1.0d0,0.0d0,ic==jc)
            wk%weight_lambda(i,j,k)=w_Jg*(dij+c(ic)*c(jc)/c_sq)
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_S_stok_on

  subroutine compute_kernelradial_D_stok_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Stokes DLP on surface: 3*(r.n)*r*r'/(4*pi*|r|^5), r = lambda*c, c = C0+lambda*c1r, n = -N/J. As in
    ! compute_kernelradial_D_on, J*(r.n) = lambda^2*JP_D, JP_D = (V1r-c1r).N, so lambda*F = 6/(8*pi)*JP_D*g*(c*c'/|c|^2)*mu
    ! with g = 1/|c|^3 (p=1) is smooth and step 8 adds nothing.
    ! f = x_E*w_s/(8*pi);  weight_lambda(nr,n_E,6) = 6*f*w_R*JP_D*g*c_i*c_j/|c|^2, the weights on mu_j(t'+lambda*d),
    ! (i,j) upper triangle (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, c(3), c_sq, JP_D, w_D

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),6))
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(8.0d0*pi)
      do i=1,nr
        c=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)
        c_sq=(c(1)*c(1)+c(2)*c(2))+c(3)*c(3)
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        w_D=6.0d0*factor*R%w(i)*(JP_D*gp%g(i,j))
        k=0
        do ic=1,3
          do jc=ic,3
            k=k+1
            wk%weight_lambda(i,j,k)=w_D*(c(ic)*c(jc))/c_sq
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_D_stok_on

  subroutine compute_kernelradial_Sn_stok_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Stokes SLPn (traction of the SLP) on surface, target normal n_r0 = N0/J0:
    ! -3*((x-y).n_r0)*r*r'/(4*pi*|r|^5), r = x-y = -lambda*c, c = C0+lambda*c1r. C0.n_r0 = 0, so (x-y).n_r0 = -lambda^2*c1r.n_r0
    ! and lambda*F = 6/(8*pi)*JP_Sn*g*(c*c'/|c|^2)*mu, JP_Sn = J*(c1r.n_r0), g = 1/|c|^3 (p=1), is smooth: step 8 adds nothing.
    ! compute_kernelradial_D_stok_on with JP_D -> JP_Sn.
    ! f = x_E*w_s/(8*pi);  weight_lambda(nr,n_E,6) = 6*f*w_R*JP_Sn*g*c_i*c_j/|c|^2, the weights on mu_j(t'+lambda*d),
    ! (i,j) upper triangle (1,1),(1,2),(1,3),(2,2),(2,3),(3,3).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, n_r0(3), c(3), c_sq, JP_Sn, w_Sn

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),6))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(8.0d0*pi)
      do i=1,nr
        c=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)
        c_sq=(c(1)*c(1)+c(2)*c(2))+c(3)*c(3)
        JP_Sn=Jx%J(i,j)*((X%c1r(i,1,j)*n_r0(1)+X%c1r(i,2,j)*n_r0(2))+X%c1r(i,3,j)*n_r0(3))
        w_Sn=6.0d0*factor*R%w(i)*(JP_Sn*gp%g(i,j))
        k=0
        do ic=1,3
          do jc=ic,3
            k=k+1
            wk%weight_lambda(i,j,k)=w_Sn*(c(ic)*c(jc))/c_sq
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Sn_stok_on

  subroutine compute_kernelradial_Dn_stok_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7, Stokes DLPn (finite-part traction of the DLP) on surface, target normal n_r0 = N0/J0, viscosity 1, the
    ! Gauss sums. Sign: the negative of Sto3dDLPmat's second output, N = r_u x r_v, source normal -N/J.
    ! K = JP_Dn*g = K_Dn0+lambda*K_Dn1+lambda^2*K_Dnr, g = B^(-7/2) (p=3), B = |c|^2 = B0+lambda*B1r, c = C0+lambda*c1r,
    !   JP_Dn_ij = 2*n_r0_i*B^2*N_j+3*B*(N.n_r0)*c_i*c_j+3*lambda*B*(cn1r*N_i*c_j+Nc1r*c_i*n_r0_j)
    !              +lambda^2*cn1r*Nc1r*(3*B*delta_ij-30*c_i*c_j),
    ! compute_kernelradial_Dn_stok_off's numerator at h = 0, tnx = n_r0: c.n_r0 = lambda*cn1r, N.c = lambda*Nc1r
    ! (C0.n_r0 = N0.C0 = 0), cn1r = cn1+lambda*cnr, Nc1r = Nc1+lambda*Ncr. With N1r = N1+lambda*Nr and N.V = 0:
    ! Nc1 = -N0.c1, Ncr = Nr.C0+N1r.c1+N.cr, no division by lambda. mu = mu0+lambda*mu1+lambda^2*mur:
    !   lambda*F = K*mur + K_Dnr*mu0 + lambda*K_Dnr*mu1 + [(K_Dn0/lambda^2+K_Dn1/lambda)*mu0+(K_Dn0/lambda+K_Dn1)*mu1],
    ! the bracket being step 8 (compute_singular_correction_on, compute_singular_correction_Dnn_on). f = x_E*w_s/(4*pi):
    !   weight_r(nr,n_E,9) = f*w_R*K (on mur),
    !   weight_0 += sum_j sum_i f*w_R*K_Dnr,  (weight_u, weight_v) += sum_j (sum_i f*w_R*lambda*K_Dnr)*d(:,j).
    ! Components (ic,jc) at ic+3*(jc-1), without symmetry.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, n_r0(3), N0(3), C0(3), c1(3), N1(3), B0, B1, cn1, Nc1, Nn0, Nn1, B_sq0, B_sq1, dij
    real(8) :: c(nr,3), N1r(nr,3), N(nr,3), Br(nr), B1r(nr), B(nr), cnr(nr), cn1r(nr), Ncr(nr), Nc1r(nr)
    real(8) :: Nnr(nr), Nn(nr), B_sqr(nr), B_sq(nr), regular_weight(nr), regular_0, regular_1
    real(8) :: U0, U1, Ur, U, V0, V1, Vr, P0, P1r, L0, L1r, JP_Dn0, JP_Dn1, JP_Dnr, JP_Dn, K_Dn0, K_Dn1, K_Dnr

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),9))
    n_r0=geo%N0/geo%J0
    N0=Jx%N0(1,:)
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(4.0d0*pi)
      C0=X%C0(:,j)
      c1=X%c1(1,:,j)
      N1=Jx%N1(1,:,j)
      B0=X%q0(j)
      B1=2.0d0*((C0(1)*c1(1)+C0(2)*c1(2))+C0(3)*c1(3))
      cn1=(c1(1)*n_r0(1)+c1(2)*n_r0(2))+c1(3)*n_r0(3)
      Nc1=-((N0(1)*c1(1)+N0(2)*c1(2))+N0(3)*c1(3))
      Nn0=(N0(1)*n_r0(1)+N0(2)*n_r0(2))+N0(3)*n_r0(3)
      Nn1=(N1(1)*n_r0(1)+N1(2)*n_r0(2))+N1(3)*n_r0(3)
      B_sq0=B0*B0
      B_sq1=2.0d0*B0*B1
      do i=1,nr
        c(i,:)=C0+R%lambda(i)*X%c1r(i,:,j)
        Br(i)=2.0d0*((X%cr(i,1,j)*C0(1)+X%cr(i,2,j)*C0(2))+X%cr(i,3,j)*C0(3)) &
          +((X%c1r(i,1,j)*X%c1r(i,1,j)+X%c1r(i,2,j)*X%c1r(i,2,j))+X%c1r(i,3,j)*X%c1r(i,3,j))
        B1r(i)=B1+R%lambda(i)*Br(i)
        B(i)=B0+R%lambda(i)*B1r(i)
        N1r(i,:)=N1+R%lambda(i)*Jx%Nr(i,:,j)
        N(i,:)=N0+R%lambda(i)*N1r(i,:)
        cnr(i)=(X%cr(i,1,j)*n_r0(1)+X%cr(i,2,j)*n_r0(2))+X%cr(i,3,j)*n_r0(3)
        cn1r(i)=cn1+R%lambda(i)*cnr(i)
        Ncr(i)=(((Jx%Nr(i,1,j)*C0(1)+Jx%Nr(i,2,j)*C0(2))+Jx%Nr(i,3,j)*C0(3)) &
          +((N1r(i,1)*c1(1)+N1r(i,2)*c1(2))+N1r(i,3)*c1(3))) &
          +((N(i,1)*X%cr(i,1,j)+N(i,2)*X%cr(i,2,j))+N(i,3)*X%cr(i,3,j))
        Nc1r(i)=Nc1+R%lambda(i)*Ncr(i)
        Nnr(i)=(Jx%Nr(i,1,j)*n_r0(1)+Jx%Nr(i,2,j)*n_r0(2))+Jx%Nr(i,3,j)*n_r0(3)
        Nn(i)=Nn0+R%lambda(i)*(Nn1+R%lambda(i)*Nnr(i))
        B_sqr(i)=2.0d0*B0*Br(i)+B1r(i)*B1r(i)
        B_sq(i)=B(i)*B(i)
      enddo
      do jc=1,3
        do ic=1,3
          k=ic+3*(jc-1)
          dij=merge(1.0d0,0.0d0,ic==jc)
          U0=C0(ic)*C0(jc)
          U1=c1(ic)*C0(jc)+C0(ic)*c1(jc)
          V0=Nn0*U0
          V1=Nn1*U0+Nn0*U1
          P0=N0(ic)*C0(jc)
          L0=cn1*P0+Nc1*C0(ic)*n_r0(jc)
          JP_Dn0=2.0d0*n_r0(ic)*B_sq0*N0(jc)+3.0d0*B0*V0
          JP_Dn1=2.0d0*n_r0(ic)*(B_sq1*N0(jc)+B_sq0*N1(jc))+3.0d0*((B1*V0+B0*V1)+B0*L0)
          K_Dn0=JP_Dn0*gp%g0(j)
          K_Dn1=JP_Dn1*gp%g0(j)+JP_Dn0*gp%g1(j)
          !GCC$ ivdep
          do i=1,nr
            Ur=(X%cr(i,ic,j)*C0(jc)+X%c1r(i,ic,j)*c1(jc))+c(i,ic)*X%cr(i,jc,j)
            U=c(i,ic)*c(i,jc)
            Vr=(Nnr(i)*U0+(Nn1+R%lambda(i)*Nnr(i))*U1)+Nn(i)*Ur
            P1r=N1r(i,ic)*C0(jc)+N(i,ic)*X%c1r(i,jc,j)
            L1r=(cnr(i)*P0+cn1r(i)*P1r)+(Ncr(i)*C0(ic)+Nc1r(i)*X%c1r(i,ic,j))*n_r0(jc)
            JP_Dnr=2.0d0*n_r0(ic)*((B_sqr(i)*N0(jc)+(B_sq1+R%lambda(i)*B_sqr(i))*N1(jc))+B_sq(i)*Jx%Nr(i,jc,j))
            JP_Dnr=(JP_Dnr+3.0d0*((((Br(i)*V0+B1r(i)*V1)+B(i)*Vr)+B1r(i)*L0)+B(i)*L1r)) &
              +cn1r(i)*Nc1r(i)*(3.0d0*B(i)*dij-30.0d0*U)
            JP_Dn=(JP_Dn0+R%lambda(i)*JP_Dn1)+R%lambda(i)**2*JP_Dnr
            K_Dnr=(JP_Dnr*gp%g0(j)+(JP_Dn1+R%lambda(i)*JP_Dnr)*gp%g1(j))+JP_Dn*gp%gr(i,j)
            regular_weight(i)=factor*R%w(i)*K_Dnr
            wk%weight_r(i,j,k)=factor*R%w(i)*(K_Dn0+R%lambda(i)*K_Dn1)+R%lambda(i)**2*regular_weight(i)
          enddo
          regular_0=0.0d0
          regular_1=0.0d0
          do i=1,nr
            regular_0=regular_0+regular_weight(i)
            regular_1=regular_1+R%lambda(i)*regular_weight(i)
          enddo
          weight_0(k)=weight_0(k)+regular_0
          weight_u(k)=weight_u(k)+regular_1*d(1,j)
          weight_v(k)=weight_v(k)+regular_1*d(2,j)
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Dn_stok_on

  subroutine compute_kernelradial_S_stok_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 7 of quadS_stok_off_w16: compute_kernelradial_Sn_off (the Laplace SLPn weights in components 1:3) with
    ! weight_lambda allocated for the 9 components (a,b) that compute_singular_correction_stok_off (kind 'S') forms.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)

    if (.not.allocated(wk%weight_lambda)) allocate(wk%weight_lambda(nr,size(d,2),9))
    call compute_kernelradial_Sn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
  end subroutine compute_kernelradial_S_stok_off

  subroutine compute_kernelradial_Dn_stok_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v,tnx)
    ! Step 7, Stokes DLPn off surface (h/=0), target traction normal tnx, viscosity 1.
    ! Sign: negative of Sto3dDLPmat's second output, source normal -N/J. With chord y-x=lambda*c,
    ! c=C0+lambda*c1r-(h/lambda)*n_r0, B=c.c, N=N0+lambda*N1r, g=B^(-7/2):
    ! JP_Dn_ij = 2*tnx_i*B^2*N_j+3*B*(N.tnx)*c_i*c_j
    !   +3*B*((c.tnx)*N_i*c_j+(N.c)*c_i*tnx_j)+(c.tnx)*(N.c)*(3*B*delta_ij-30*c_i*c_j).
    ! weight_r=f*w_R*JP_Dn*g is on mur; step 8 forms the endpoints and corrects the radial moments.
    ! Components (ic,jc) at ic+3*(jc-1), without symmetry; f=x_E*w_s/(4*pi).
    integer(8), intent(in) :: iedge,n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:),w_s(:),tnx(3)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:),weight_u(:),weight_v(:)
    integer(8) :: i,j,ic,jc,k
    real(8) :: n_r0(3),factor,c(nr,3),N1r(3),N(nr,3),B(nr),ct(nr),Nc(nr),Nt(nr),w_g(nr),U,JP_Dn,dij

    if (.not.allocated(wk%weight_r)) allocate(wk%weight_r(nr,size(d,2),9))
    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(4.0d0*pi)
      do i=1,nr
        c(i,:)=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)-(geo%h*R%lambda_inv(i))*n_r0
        N1r=Jx%N1(1,:,j)+R%lambda(i)*Jx%Nr(i,:,j)
        N(i,:)=Jx%N0(1,:)+R%lambda(i)*N1r
        B(i)=(c(i,1)*c(i,1)+c(i,2)*c(i,2))+c(i,3)*c(i,3)
        ct(i)=(c(i,1)*tnx(1)+c(i,2)*tnx(2))+c(i,3)*tnx(3)
        Nc(i)=(N(i,1)*c(i,1)+N(i,2)*c(i,2))+N(i,3)*c(i,3)
        Nt(i)=(N(i,1)*tnx(1)+N(i,2)*tnx(2))+N(i,3)*tnx(3)
        w_g(i)=factor*R%w(i)*gp%g(i,j)
      enddo
      do jc=1,3
        do ic=1,3
          k=ic+3*(jc-1)
          dij=merge(1.0d0,0.0d0,ic==jc)
          !GCC$ ivdep
          do i=1,nr
            U=c(i,ic)*c(i,jc)
            JP_Dn=2.0d0*tnx(ic)*B(i)**2*N(i,jc)+3.0d0*B(i)*Nt(i)*U
            JP_Dn=JP_Dn+3.0d0*B(i)*(ct(i)*N(i,ic)*c(i,jc)+Nc(i)*c(i,ic)*tnx(jc))
            JP_Dn=JP_Dn+ct(i)*Nc(i)*(3.0d0*B(i)*dij-30.0d0*U)
            wk%weight_r(i,j,k)=w_g(i)*JP_Dn
          enddo
        enddo
      enddo
    enddo
  end subroutine compute_kernelradial_Dn_stok_off

  subroutine compute_singular_correction_on(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 8 on surface: sum_j f*reg-int_0^1 Fhat, from the lambda^0, lambda^1 data at t',
    !   reg-int lambda^-1 = log(sqrt(q0)) = log_gamma,  reg-int lambda^-2 = kappa-1,  kappa = beta1/2.
    ! By order: 0 (SLP, DLP) no Fhat; 1 compute_singular_correction_Sn_on; 2 compute_singular_correction_Dn_on (p=2),
    ! compute_singular_correction_Dnn_on (p=3, Stokes DLPn). wk is unused; the arguments are those of every driver's step 8.
    integer(8), intent(in) :: p, order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)

    if (order==1_8) then
      call compute_singular_correction_Sn_on(iedge,geo,n_E,w_s,X,Jx,gp,weight_0)
    elseif (order==2_8 .and. p==3_8) then
      call compute_singular_correction_Dnn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,weight_0,weight_u,weight_v)
    elseif (order==2_8) then
      call compute_singular_correction_Dn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,weight_0,weight_u,weight_v)
    endif
  end subroutine compute_singular_correction_on

  subroutine compute_singular_correction_Sn_on(iedge,geo,n_E,w_s,X,Jx,gp,weight_0)
    ! Order 1 (SLPn): Fhat = C0*J0*g0*mu0/lambda;  weight_0 += f*J0*g0*log_gamma*C0.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    real(8), intent(in) :: w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    real(8), intent(inout) :: weight_0(:)
    integer(8) :: j
    real(8) :: factor, log_gamma

    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      log_gamma=log(X%q0(j))/2.0d0
      weight_0=weight_0+factor*Jx%J0*gp%g0(j)*log_gamma*X%C0(:,j)
    enddo
  end subroutine compute_singular_correction_Sn_on

  subroutine compute_singular_correction_Dn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,weight_0,weight_u,weight_v)
    ! Order 2 (DLPn): Fhat = (K0/lambda^2+K1/lambda)*mu0+(K0/lambda+K1)*mu1, K0, K1 as in compute_kernelradial_Dn_on;
    !   e0 = (kappa-1)*K0+log_gamma*K1, e1 = sum(w_R)*K1+log_gamma*K0;
    !   weight_0 += f*e0,  (weight_u, weight_v) += f*e1*d(:,j).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, log_gamma, kappa, sum_w, C0(3), c1(3), B0, B1, dij, P_Dn0, P_Dn1, JP_Dn0, JP_Dn1
    real(8) :: K_Dn0, K_Dn1, e0, e1

    sum_w=0.0d0
    do i=1,nr
      sum_w=sum_w+R%w(i)
    enddo
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      C0=X%C0(:,j)
      c1=X%c1(1,:,j)
      B0=X%q0(j)
      B1=2.0d0*((C0(1)*c1(1)+C0(2)*c1(2))+C0(3)*c1(3))
      log_gamma=log(B0)/2.0d0
      kappa=0.5d0*gp%beta1(j)
      k=0
      do ic=1,3
        do jc=ic,3
          k=k+1
          dij=merge(1.0d0,0.0d0,ic==jc)
          P_Dn0=3.0d0*C0(ic)*C0(jc)-dij*B0
          P_Dn1=3.0d0*(c1(ic)*C0(jc)+C0(ic)*c1(jc))-dij*B1
          JP_Dn0=Jx%J0*P_Dn0
          JP_Dn1=Jx%J1(j)*P_Dn0+Jx%J0*P_Dn1
          K_Dn0=JP_Dn0*gp%g0(j)
          K_Dn1=JP_Dn1*gp%g0(j)+JP_Dn0*gp%g1(j)
          e0=(kappa-1.0d0)*K_Dn0+log_gamma*K_Dn1
          e1=sum_w*K_Dn1+log_gamma*K_Dn0
          weight_0(k)=weight_0(k)+factor*e0
          weight_u(k)=weight_u(k)+factor*e1*d(1,j)
          weight_v(k)=weight_v(k)+factor*e1*d(2,j)
        enddo
      enddo
    enddo
  end subroutine compute_singular_correction_Dn_on

  subroutine compute_singular_correction_Dnn_on(iedge,geo,R,n_E,d,w_s,X,Jx,gp,weight_0,weight_u,weight_v)
    ! Order 2, p=3 (Stokes DLPn): Fhat = (K_Dn0/lambda^2+K_Dn1/lambda)*mu0+(K_Dn0/lambda+K_Dn1)*mu1, K_Dn0, K_Dn1 as in
    !   compute_kernelradial_Dn_stok_on; e0 = (kappa-1)*K_Dn0+log_gamma*K_Dn1, e1 = sum(w_R)*K_Dn1+log_gamma*K_Dn0
    !   (physical-distance finite part, reg-int lambda^-1 = log_gamma = log(sqrt(q0)), reg-int lambda^-2 = kappa-1,
    !   kappa = beta1/2); f = x_E*w_s/(4*pi), weight_0 += f*e0, (weight_u, weight_v) += f*e1*d(:,j).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ic, jc
    real(8) :: factor, sum_w, n_r0(3), N0(3), C0(3), c1(3), N1(3), B0, B1, cn1, Nc1, Nn0, Nn1, B_sq0, B_sq1
    real(8) :: log_gamma, kappa, U0, U1, V0, V1, P0, L0, JP_Dn0, JP_Dn1, K_Dn0, K_Dn1, e0, e1

    sum_w=0.0d0
    do i=1,nr
      sum_w=sum_w+R%w(i)
    enddo
    n_r0=geo%N0/geo%J0
    N0=Jx%N0(1,:)
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(4.0d0*pi)
      C0=X%C0(:,j)
      c1=X%c1(1,:,j)
      N1=Jx%N1(1,:,j)
      B0=X%q0(j)
      B1=2.0d0*((C0(1)*c1(1)+C0(2)*c1(2))+C0(3)*c1(3))
      cn1=(c1(1)*n_r0(1)+c1(2)*n_r0(2))+c1(3)*n_r0(3)
      Nc1=-((N0(1)*c1(1)+N0(2)*c1(2))+N0(3)*c1(3))
      Nn0=(N0(1)*n_r0(1)+N0(2)*n_r0(2))+N0(3)*n_r0(3)
      Nn1=(N1(1)*n_r0(1)+N1(2)*n_r0(2))+N1(3)*n_r0(3)
      B_sq0=B0*B0
      B_sq1=2.0d0*B0*B1
      log_gamma=log(B0)/2.0d0
      kappa=0.5d0*gp%beta1(j)
      do jc=1,3
        do ic=1,3
          k=ic+3*(jc-1)
          U0=C0(ic)*C0(jc)
          U1=c1(ic)*C0(jc)+C0(ic)*c1(jc)
          V0=Nn0*U0
          V1=Nn1*U0+Nn0*U1
          P0=N0(ic)*C0(jc)
          L0=cn1*P0+Nc1*C0(ic)*n_r0(jc)
          JP_Dn0=2.0d0*n_r0(ic)*B_sq0*N0(jc)+3.0d0*B0*V0
          JP_Dn1=2.0d0*n_r0(ic)*(B_sq1*N0(jc)+B_sq0*N1(jc))+3.0d0*((B1*V0+B0*V1)+B0*L0)
          K_Dn0=JP_Dn0*gp%g0(j)
          K_Dn1=JP_Dn1*gp%g0(j)+JP_Dn0*gp%g1(j)
          e0=(kappa-1.0d0)*K_Dn0+log_gamma*K_Dn1
          e1=sum_w*K_Dn1+log_gamma*K_Dn0
          weight_0(k)=weight_0(k)+factor*e0
          weight_u(k)=weight_u(k)+factor*e1*d(1,j)
          weight_v(k)=weight_v(k)+factor*e1*d(2,j)
        enddo
      enddo
    enddo
  end subroutine compute_singular_correction_Dnn_on

  subroutine compute_singular_correction_off(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v,tnx)
    ! Step 8. Off-surface epsilon-moment correction (h/=0).
    ! eps = |h|/sqrt(q0) per ray. y^-(2p+1) = sum_m c_y(m+1,p+1)*d_curv^m*(lambda^2/(lambda^2+eps^2))^(m+p+1/2).
    ! For m = 0..geo%K the Gauss sum of step 7 is corrected to the exact integral of P(lambda)*(...)^(m+p+1/2),
    ! P of degree < nr; m > K stays with step 7.
    ! compute_moment_correction_off (shared, by p): the moment deficits of the corrected rays and their V^-1 solves, mc.
    ! Then by kernel: compute_singular_correction_S_off (p=0), _D_off (p=1, ncomp=1), _Sn_off (p=1, ncomp=3),
    ! _Dn_off (p=2), _Dnn_off (p=3, R^-7), ncomp = size(wk%weight_lambda,3).
    ! order is unused; tnx is the target traction normal, required only for p=3.
    integer(8), intent(in) :: p, order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    real(8), intent(in), optional :: tnx(3)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    type(moment_correction) :: mc

    call compute_moment_correction_off(p,geo,R,n_E,X,mc)
    if (p==0_8) then
      call compute_singular_correction_S_off(iedge,geo,R,n_E,w_s,Jx,gp,mc,wk)
    elseif (p==1_8) then
      if (size(wk%weight_lambda,3)==1) then
        call compute_singular_correction_D_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0)
      else
        call compute_singular_correction_Sn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0)
      endif
    elseif (p==2_8) then
      call compute_singular_correction_Dn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0,weight_u,weight_v)
    elseif (p==3_8) then
      call compute_singular_correction_Dnn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0,weight_u,weight_v,tnx)
    endif
  end subroutine compute_singular_correction_off

  subroutine compute_moment_correction_off(p,geo,R,n_E,X,mc)
    ! Shared part of step 8 off: eps per ray and, for every corrected ray (eps < epsilon_big) and m = 0..K_m, the V^-1
    ! solve of its moment deficit; K_m = K, K+1 (at most 5) for p=2, and K+1 for p=3. Moments M_m(k) of
    ! k_m = lambda^(2m+1)/(lambda^2+eps^2)^(m+3/2), k = 0..nr+1, and DeltaM = M - Lj^T*kernel (kernel = w*k_m).
    ! Rays: eps <= epsilon_small (1e-4*min(lambda)): DeltaM from its closed form through eps^4, V^-1 through ell0_coeff;
    !   epsilon_small < eps < epsilon_big (0.5; 0.3 if |h| > 1e-2): analytic moments, one solve through U_V, S_inv, W_V.
    ! mc%rays: the corrected rays, closed-form ones first; ray mc%rays(e), term m is column (K_m+1)*(e-1)+m+1. By p:
    ! p=0: mc%sigma = V^-1*(DeltaM_m(2:nr+1)+eps^2*DeltaM_m(0:nr-1)).
    ! p=1: mc%omega = V^-1*DeltaM_m(2:nr+1), mc%alpha = DeltaM_m(1)-sum(omega_m/lambda), mc%z = omega_m/lambda^2
    !   +alpha_m*ell0/lambda, and h*(DeltaM_m(0)-sum(z_m)) = mc%h_scale*mc%rest0: the closed form gives h*DeltaM_m(0)
    !   itself (h_scale = 1), the middle range DeltaM_m(0) (h_scale = h).
    ! p=2, per ray j: the p=1 moments to k = nr+3, the p=2 moments N_m to k = nr+4, grad_x DeltaM_m(k)
    !   = -(a/q0)*(h*n_r0*DeltaN_m(k)-C0*DeltaN_m(k+1)), a = 2m+3; mc%omega_shift(:,m+1,s+1,j) = V^-1*DeltaM_m(2+s:nr+1+s),
    !   s = 0,1,2, mc%omega_shift_gradient, mc%DeltaM01 = DeltaM_m(0:1), mc%DeltaM01_gradient.
    ! p=3, per ray j: moments P_m of k_m/(lambda^2+eps^2)^2 to k = nr+K_m+5, from their closed form for tiny eps
    !   and P_m(k+2) = N_m(k)-eps^2*P_m(k) otherwise. mc%omega_shift(:,m+1,s+1,j) = V^-1*DeltaP_m(s:nr-1+s),
    !   s = 0..K_m+6; step 8 supplies the h powers. Both ranges use the same SVD solve.
    ! The middle-range moment arrays are automatic, sized (K_m+1)*(middle rays); mc is allocated per call.
    integer(8), intent(in) :: p, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    type(displacement_expansion), intent(in) :: X
    type(moment_correction), intent(inout) :: mc
    real(8), parameter :: beta_m(6)=[1.0d0,2.0d0/3,8.0d0/15,16.0d0/35,128.0d0/315,256.0d0/693]
    real(8), parameter :: Hodd(6)=[1.0d0,4.0d0/3,23.0d0/15,176.0d0/105,563.0d0/315,6508.0d0/3465]
    real(8), parameter :: odd_inverse(7)=[1.0d0,1.0d0/3,1.0d0/5,1.0d0/7,1.0d0/9,1.0d0/11,1.0d0/13]
    integer(8) :: i, j, e, im, jexp, l, k, icol, col, ncol, n_tiny, n_sel, tiny(n_E), rays(n_E), ic, K_m, n_mom, ishift
    integer(8) :: offset, n_N, n_P
    real(8) :: epsilon(n_E), epsilon_small, epsilon_big, eps_sq, eps_sq_sq, log_inverse, n_r0(3)
    real(8) :: a, t, gamma_m, dS(4), dM0_h, dM1, dMg(4), alpha, zsum, omega_i, z_i
    real(8) :: Dlog, dMg_eps(4), DeltaN(9), grad_M(8,3), dMg_shift(6)
    real(8) :: beta_P
    external :: dgemm

    K_m=geo%K
    if (p==2_8) K_m=min(geo%K+1_8,5_8)
    if (p==3_8) K_m=geo%K+1_8
    epsilon_small=1.0d-4*minval(R%lambda)
    epsilon_big=0.5d0
    if (abs(geo%h)>1.0d-2) epsilon_big=0.3d0
    n_r0=geo%N0/geo%J0
    n_tiny=0_8
    n_sel=0_8
    do j=1,n_E
      epsilon(j)=abs(geo%h)/sqrt(X%q0(j))
      if (epsilon(j)<=epsilon_small) then
        n_tiny=n_tiny+1_8
        tiny(n_tiny)=j
      elseif (epsilon(j)<epsilon_big) then
        n_sel=n_sel+1_8
        rays(n_sel)=j
      endif
    enddo
    mc%K_m=K_m
    mc%epsilon_big=epsilon_big
    mc%epsilon=epsilon(1:n_E)
    mc%rays=[tiny(1:n_tiny),rays(1:n_sel)]
    ncol=(K_m+1_8)*(n_tiny+n_sel)
    if (p==0_8) then
      allocate(mc%sigma(nr,ncol))
    elseif (p==1_8) then
      allocate(mc%omega(nr,ncol),mc%z(nr,ncol),mc%alpha(ncol),mc%rest0(ncol),mc%h_scale(ncol))
    elseif (p==2_8) then   ! per ray: V^-1 of DeltaM shifted by s = 0,1,2 (s = 0: omega), DeltaM(0:1), and their gradients
      allocate(mc%omega_shift(nr,K_m+1,3,n_E),mc%omega_shift_gradient(nr,K_m+1,3,3,n_E))
      allocate(mc%DeltaM01(2,K_m+1,n_E),mc%DeltaM01_gradient(2,K_m+1,3,n_E))
    else   ! p=3, per ray: V^-1 of DeltaP shifted by s = 0..K_m+6
      allocate(mc%omega_shift(nr,K_m+1,K_m+7,n_E))
      mc%omega_shift=0.0d0
    endif

    ! eps <= epsilon_small: closed form
    do e=1,n_tiny
      j=tiny(e)
      eps_sq=epsilon(j)*epsilon(j)
      eps_sq_sq=eps_sq*eps_sq
      log_inverse=log(2.0d0)-log(epsilon(j))
      do im=0,K_m
        col=(K_m+1_8)*(e-1_8)+im+1_8
        a=real(im+min(p,1_8),8)+0.5d0   ! m+p+1/2; DLPn differentiates the p=1 moments
        t=0.5d0*a*(a+1.0d0)
        if (p==0_8) then
          gamma_m=real(2*im+1,8)*beta_m(im+1)
          dS(1)=-gamma_m*epsilon(j)+a*eps_sq*(1.0d0+R%radial_inv_sums(2))-t*eps_sq_sq*(1.0d0/3+R%radial_inv_sums(4))
          dS(2)=a*eps_sq*(R%radial_inv_sums(1)-log_inverse+Hodd(im+1)-0.5d0)-t*eps_sq_sq*(0.5d0+R%radial_inv_sums(3))
          dS(3)=real(2*im+2,8)*gamma_m*(eps_sq*epsilon(j))/3.0d0-t*eps_sq_sq*(1.0d0+R%radial_inv_sums(2))
          dS(4)=eps_sq_sq*(t*(log_inverse-Hodd(im+1)-R%radial_inv_sums(1))+a*(3.0d0*a+1.0d0)/8.0d0)
          !GCC$ ivdep
          do i=1,nr
            mc%sigma(i,col)=((dS(1)*R%ell0_coeff(i,1)+dS(2)*R%ell0_coeff(i,2))+dS(3)*R%ell0_coeff(i,3))+dS(4)*R%ell0_coeff(i,4)
          enddo
        elseif (p<=2_8) then
          gamma_m=real(2*im+2,8)*beta_m(im+1)
          dM0_h=sign(sqrt(X%q0(j)),geo%h)*beta_m(im+1)+geo%h*(-(1.0d0+R%radial_inv_sums(2)) &
            +a*eps_sq*(1.0d0/3+R%radial_inv_sums(4))-t*eps_sq_sq*(1.0d0/5+R%radial_inv_sums(6)))   ! h*DeltaM(0)
          dM1=log_inverse-Hodd(im+1)-R%radial_inv_sums(1)+a*eps_sq*(0.5d0+R%radial_inv_sums(3)) &
            -t*eps_sq_sq*(0.25d0+R%radial_inv_sums(5))                                                 ! DeltaM(1)
          dMg(1)=-gamma_m*epsilon(j)+a*eps_sq*(1.0d0+R%radial_inv_sums(2))-t*eps_sq_sq*(1.0d0/3+R%radial_inv_sums(4))
          dMg(2)=eps_sq*(a*(R%radial_inv_sums(1)-log_inverse+Hodd(im+1))-0.5d0*(a-1.0d0)) &
            -t*eps_sq_sq*(0.5d0+R%radial_inv_sums(3))
          dMg(3)=real(2*im+4,8)*gamma_m*(eps_sq*epsilon(j))/3.0d0-t*eps_sq_sq*(1.0d0+R%radial_inv_sums(2))
          dMg(4)=eps_sq_sq*(t*(log_inverse-Hodd(im+1)-R%radial_inv_sums(1))+(a-1.0d0)*(3.0d0*a+2.0d0)/8.0d0)
          if (p==1_8) then   ! DLP, SLPn
            alpha=dM1-(((dMg(1)*R%ell0_inv_sums(1,1)+dMg(2)*R%ell0_inv_sums(1,2))+dMg(3)*R%ell0_inv_sums(1,3)) &
              +dMg(4)*R%ell0_inv_sums(1,4))
            zsum=(((dMg(1)*R%ell0_inv_sums(2,1)+dMg(2)*R%ell0_inv_sums(2,2))+dMg(3)*R%ell0_inv_sums(2,3)) &
              +dMg(4)*R%ell0_inv_sums(2,4))+alpha*R%ell0_inv_sums(1,1)
            !GCC$ ivdep
            do i=1,nr
              omega_i=((dMg(1)*R%ell0_coeff(i,1)+dMg(2)*R%ell0_coeff(i,2))+dMg(3)*R%ell0_coeff(i,3))+dMg(4)*R%ell0_coeff(i,4)
              mc%omega(i,col)=omega_i
              mc%z(i,col)=omega_i*R%lambda_inv(i)*R%lambda_inv(i)+alpha*R%ell0_over_lambda(i)
            enddo
            mc%alpha(col)=alpha
            mc%rest0(col)=dM0_h-geo%h*zsum
            mc%h_scale(col)=1.0d0
          else               ! DLPn: differentiating the series term by term in eps gives the p=2 moments
            Dlog=log_inverse-Hodd(im+1)-R%radial_inv_sums(1)
            dMg_eps(1)=-gamma_m+2.0d0*a*epsilon(j)*(1.0d0+R%radial_inv_sums(2)) &
              -4.0d0*t*(eps_sq*epsilon(j))*(1.0d0/3+R%radial_inv_sums(4))
            dMg_eps(2)=epsilon(j)*(1.0d0-2.0d0*a*Dlog)-4.0d0*t*(eps_sq*epsilon(j))*(0.5d0+R%radial_inv_sums(3))
            dMg_eps(3)=real(2*im+4,8)*gamma_m*eps_sq-4.0d0*t*(eps_sq*epsilon(j))*(1.0d0+R%radial_inv_sums(2))
            dMg_eps(4)=(eps_sq*epsilon(j))*(4.0d0*t*Dlog+(a-1.0d0)*(3.0d0*a+2.0d0)/2.0d0-t)
            DeltaN=0.0d0   ! p=2 moments, k = 0..8 (k = 8 is zero at this order)
            DeltaN(1)=beta_m(im+1)/(real(2*im+3,8)*(eps_sq*epsilon(j)))-(1.0d0/3+R%radial_inv_sums(4)) &
              +(a+1.0d0)*eps_sq*(1.0d0/5+R%radial_inv_sums(6))
            DeltaN(2)=1.0d0/(real(2*im+3,8)*eps_sq)-(0.5d0+R%radial_inv_sums(3))+(a+1.0d0)*eps_sq*(0.25d0+R%radial_inv_sums(5))
            DeltaN(3:6)=-dMg_eps/(real(2*im+3,8)*epsilon(j))
            DeltaN(7)=dMg(3)-eps_sq*DeltaN(5)
            DeltaN(8)=dMg(4)-eps_sq*DeltaN(6)
            do k=1,8   ! grad_x DeltaM(k-1)
              do ic=1,3
                grad_M(k,ic)=-real(2*im+3,8)/X%q0(j)*(geo%h*DeltaN(k)*n_r0(ic)-DeltaN(k+1)*X%C0(ic,j))
              enddo
            enddo
            dMg_shift=0.0d0
            dMg_shift(1:4)=dMg
            do ishift=1,3
              do i=1,nr
                mc%omega_shift(i,im+1,ishift,j)=((dMg_shift(ishift)*R%ell0_coeff(i,1)+dMg_shift(ishift+1)*R%ell0_coeff(i,2)) &
                  +dMg_shift(ishift+2)*R%ell0_coeff(i,3))+dMg_shift(ishift+3)*R%ell0_coeff(i,4)
                do ic=1,3
                  mc%omega_shift_gradient(i,im+1,ishift,ic,j)=((R%ell0_coeff(i,1)*grad_M(ishift+2,ic) &
                    +R%ell0_coeff(i,2)*grad_M(ishift+3,ic))+R%ell0_coeff(i,3)*grad_M(ishift+4,ic)) &
                    +R%ell0_coeff(i,4)*grad_M(ishift+5,ic)
                enddo
              enddo
            enddo
            mc%DeltaM01(1,im+1,j)=dM0_h/geo%h
            mc%DeltaM01(2,im+1,j)=dM1
            mc%DeltaM01_gradient(1:2,im+1,1:3,j)=grad_M(1:2,1:3)
          endif
        else   ! p=3: DeltaP_m(k), k = 0..nr+K_m+5, then V^-1 of its K_m+7 shifts
          block
            integer(8) :: i_log, i_series, k_gauss
            real(8) :: DeltaP(nr+K_m+6), rhs(nr,K_m+7), projection(nr,K_m+7), sol(nr,K_m+7)
            real(8) :: binom_a, harmonic, odd_sum, constant
            a=real(im,8)+3.5d0
            beta_P=0.2d0   ! B(5/2,m+1)/2
            do l=1,im
              beta_P=beta_P*real(2*l,8)/real(2*l+5,8)
            enddo
            do k=0,nr+K_m+5
              if (k>=5_8 .and. mod(k,2_8)==1_8) then   ! log(1/eps) terms
                i_log=(k-5_8)/2_8
                binom_a=1.0d0
                do l=1,i_log
                  binom_a=-binom_a*(a+real(l-1,8))/real(l,8)
                enddo
                harmonic=0.0d0
                do l=1,i_log
                  harmonic=harmonic+1.0d0/real(l,8)
                enddo
                odd_sum=0.0d0
                do l=1,2*(im+i_log+2)+1,2
                  odd_sum=odd_sum+1.0d0/real(l,8)
                enddo
                DeltaP(k+1)=epsilon(j)**(2*i_log)*binom_a &
                  *(log_inverse+harmonic/2.0d0-odd_sum-R%radial_inv_sums(1))
              elseif (mod(k,2_8)==0_8) then
                constant=beta_P
                do l=0,k-2,2
                  constant=constant*real(2*im+l+2,8)/real(3-l,8)
                enddo
                DeltaP(k+1)=epsilon(j)**(k-5)*constant
              elseif (k==1_8) then
                DeltaP(k+1)=2.0d0/(real((2*im+3)*(2*im+5),8)*eps_sq_sq)
              else   ! k=3
                DeltaP(k+1)=1.0d0/(real(2*im+5,8)*eps_sq)
              endif
              binom_a=1.0d0
              do i_series=0,4   ! Gauss sums of the singular powers lambda^(k-6-2*i_series)
                k_gauss=k-6_8-2_8*i_series
                if (k_gauss<-1_8) then
                  DeltaP(k+1)=DeltaP(k+1)+epsilon(j)**(2*i_series)*binom_a &
                    *(1.0d0/real(k_gauss+1,8)-sum(R%w*R%lambda**k_gauss))
                endif
                binom_a=-binom_a*(a+real(i_series,8))/real(i_series+1,8)
              enddo
            enddo
            do ishift=0,K_m+6
              rhs(:,ishift+1)=DeltaP(ishift+1:ishift+nr)
            enddo
            call dgemm('T','N',nr,K_m+7_8,nr,1.0d0,R%U_V,nr,rhs,nr,0.0d0,projection,nr)
            do ishift=1,K_m+7
              do i=1,nr
                projection(i,ishift)=R%S_inv(i)*projection(i,ishift)
              enddo
            enddo
            call dgemm('N','N',nr,K_m+7_8,nr,1.0d0,R%W_V,nr,projection,nr,0.0d0,sol,nr)
            do ishift=1,K_m+7
              mc%omega_shift(:,im+1,ishift,j)=sol(:,ishift)
            enddo
          end block
        endif
      enddo
    enddo

    ! epsilon_small < eps < epsilon_big: analytic moments, one solve for all of them
    if (n_sel>0_8) then
      ncol=(K_m+1_8)*n_sel
      offset=(K_m+1_8)*n_tiny   ! mc columns before the middle range
      n_mom=nr+2_8
      if (p==2_8) then
        n_mom=nr+4_8
        n_N=nr+5_8
      elseif (p==3_8) then
        n_mom=nr+K_m+2_8
        n_N=nr+K_m+4_8
        n_P=nr+K_m+6_8
      endif
      block
        real(8) :: M(n_mom,ncol), kernel(nr,ncol), Gauss_sum(n_mom,ncol), DeltaS(nr,ncol), projection(nr,ncol)
        real(8) :: lambda_sq(nr), ratio(nr), denom, eps_e, s_sq, s, s_eps, s_odd(8), asinh_inverse, base, horner_at0
        real(8) :: horner, moment, a_N, grad_M(nr+4), binom_sum, power
        real(8), allocatable :: N(:,:), kernel_N(:,:), Gauss_N(:,:), rhs(:,:,:), projection_rhs(:,:,:), sol(:,:,:)
        real(8), allocatable :: P_mom(:,:), kernel_P(:,:), Gauss_P(:,:)
        if (p>=2_8) allocate(N(n_N,ncol),kernel_N(nr,ncol))
        if (p==3_8) allocate(P_mom(n_P,ncol),kernel_P(nr,ncol))
        do i=1,nr
          lambda_sq(i)=R%lambda(i)*R%lambda(i)
        enddo
        do e=1,n_sel
          eps_e=epsilon(rays(e))
          eps_sq=eps_e*eps_e
          s_sq=1.0d0/(1.0d0+eps_sq)
          s=sqrt(s_sq)
          s_eps=eps_e*s
          s_odd(1)=s
          do l=1,K_m+merge(1_8,0_8,p>=2_8)
            s_odd(l+1)=s_odd(l)*s_sq
          enddo
          asinh_inverse=asinh(1.0d0/eps_e)
          icol=1_8+(K_m+1_8)*(e-1_8)
          do i=1,nr
            denom=lambda_sq(i)+eps_sq
            ratio(i)=lambda_sq(i)/denom
            kernel(i,icol)=R%w(i)*R%lambda(i)/(denom*sqrt(denom))
          enddo
          base=s
          horner_at0=1.0d0
          do jexp=0,K_m
            icol=jexp+1_8+(K_m+1_8)*(e-1_8)
            if (jexp>0) horner_at0=horner_at0*real(2*jexp,8)*odd_inverse(jexp+1)
            horner=(-1.0d0)**jexp*R%binom(jexp+1,jexp+1)*odd_inverse(jexp+1)
            do l=jexp-1,0,-1
              horner=horner*s_eps*s_eps+(-1.0d0)**l*R%binom(l+1,jexp+1)*odd_inverse(l+1)
            enddo
            moment=asinh_inverse
            do l=0,jexp
              moment=moment-s_odd(l+1)*odd_inverse(l+1)
            enddo
            M(1,icol)=(horner_at0-s_eps*horner)/eps_e
            M(2,icol)=moment
            do k=0,n_mom-3
              M(k+3,icol)=(base-real(2*jexp+2+k,8)*eps_sq*M(k+1,icol))*(1.0d0/real(k+1,8))
            enddo
            base=base*s_sq
            if (jexp<K_m) then
              do i=1,nr
                kernel(i,icol+1)=kernel(i,icol)*ratio(i)
              enddo
            endif
            if (p>=2_8) then   ! p=2 moments of k_m/(lambda^2+eps^2)
              a_N=real(2*jexp+3,8)
              N(1,icol)=(M(1,icol)+s_odd(jexp+2))/(a_N*eps_sq)
              N(2,icol)=s_odd(jexp+2)/(a_N*eps_sq)
              do k=1,n_N-2
                N(k+2,icol)=M(k,icol)-eps_sq*N(k,icol)
              enddo
              do i=1,nr
                kernel_N(i,icol)=kernel(i,icol)/(lambda_sq(i)+eps_sq)
              enddo
            endif
            if (p==3_8) then   ! p=3 moments of k_m/(lambda^2+eps^2)^2
              beta_P=0.2d0
              do l=1,jexp
                beta_P=beta_P*real(2*l,8)/real(2*l+5,8)
              enddo
              binom_sum=0.0d0
              power=s_eps**5
              do l=0,jexp
                binom_sum=binom_sum+(-1.0d0)**l*R%binom(l+1,jexp+1)*power/real(2*l+5,8)
                power=power*s_eps*s_eps
              enddo
              P_mom(1,icol)=(beta_P-binom_sum)/(eps_sq*eps_sq*eps_e)
              P_mom(2,icol)=s_odd(jexp+2)*(2.0d0+real(2*jexp+3,8)*s_eps*s_eps) &
                /(real((2*jexp+3)*(2*jexp+5),8)*eps_sq*eps_sq)
              do k=1,n_P-2
                P_mom(k+2,icol)=N(k,icol)-eps_sq*P_mom(k,icol)
              enddo
              do i=1,nr
                kernel_P(i,icol)=kernel_N(i,icol)/(lambda_sq(i)+eps_sq)
              enddo
            endif
          enddo
        enddo
        call dgemm('T','N',n_mom,ncol,nr,1.0d0,R%Lj,nr,kernel,nr,0.0d0,Gauss_sum,n_mom)
        if (p==0_8) then
          ! M becomes DeltaM; DeltaS = DeltaM(3:nr+2)+eps^2*DeltaM(1:nr)
          do e=1,n_sel
            eps_sq=epsilon(rays(e))*epsilon(rays(e))
            do icol=(geo%K+1_8)*(e-1_8)+1_8,(geo%K+1_8)*e
              do i=1,nr+2
                M(i,icol)=M(i,icol)-Gauss_sum(i,icol)
              enddo
              do i=1,nr
                DeltaS(i,icol)=M(i+2,icol)+eps_sq*M(i,icol)
              enddo
            enddo
          enddo
          call dgemm('T','N',nr,ncol,nr,1.0d0,R%U_V,nr,DeltaS,nr,0.0d0,projection,nr)
          do icol=1,ncol
            do i=1,nr
              projection(i,icol)=R%S_inv(i)*projection(i,icol)
            enddo
          enddo
          call dgemm('N','N',nr,ncol,nr,1.0d0,R%W_V,nr,projection,nr,0.0d0,mc%sigma(1,offset+1),nr)
        elseif (p==1_8) then
          ! M becomes DeltaM; omega = V^-1*DeltaM(3:nr+2)
          do icol=1,ncol
            do i=1,nr+2
              M(i,icol)=M(i,icol)-Gauss_sum(i,icol)
            enddo
          enddo
          call dgemm('T','N',nr,ncol,nr,1.0d0,R%U_V,nr,M(3,1),n_mom,0.0d0,projection,nr)
          do icol=1,ncol
            do i=1,nr
              projection(i,icol)=R%S_inv(i)*projection(i,icol)
            enddo
          enddo
          call dgemm('N','N',nr,ncol,nr,1.0d0,R%W_V,nr,projection,nr,0.0d0,mc%omega(1,offset+1),nr)
          do icol=1,ncol
            col=offset+icol
            alpha=M(2,icol)-sum(mc%omega(:,col)*R%lambda_inv)
            zsum=0.0d0
            do i=1,nr
              z_i=(mc%omega(i,col)*R%lambda_inv(i)+alpha*R%ell0(i))*R%lambda_inv(i)
              mc%z(i,col)=z_i
              zsum=zsum+z_i
            enddo
            mc%alpha(col)=alpha
            mc%rest0(col)=M(1,icol)-zsum
            mc%h_scale(col)=geo%h
          enddo
        elseif (p==2_8) then
          ! p=2: M becomes DeltaM, N becomes DeltaN; one solve for the 3 shifts and the 3x3 shifted gradients
          do icol=1,ncol
            do i=1,n_mom
              M(i,icol)=M(i,icol)-Gauss_sum(i,icol)
            enddo
          enddo
          allocate(Gauss_N(nr+5,ncol))
          call dgemm('T','N',nr+5_8,ncol,nr,1.0d0,R%Lj,nr,kernel_N,nr,0.0d0,Gauss_N,nr+5_8)
          do icol=1,ncol
            do i=1,nr+5
              N(i,icol)=N(i,icol)-Gauss_N(i,icol)
            enddo
          enddo
          allocate(rhs(nr,ncol,12),projection_rhs(nr,ncol,12),sol(nr,ncol,12))
          do e=1,n_sel
            j=rays(e)
            do jexp=0,K_m
              icol=jexp+1_8+(K_m+1_8)*(e-1_8)
              a_N=real(2*jexp+3,8)
              do ishift=1,3
                do i=1,nr
                  rhs(i,icol,ishift)=M(i+ishift+1,icol)
                enddo
              enddo
              do ic=1,3
                do k=1,nr+4   ! grad_x DeltaM(k-1)
                  grad_M(k)=-(geo%h*n_r0(ic)*N(k,icol)-X%C0(ic,j)*N(k+1,icol))*(a_N/X%q0(j))
                enddo
                mc%DeltaM01_gradient(1:2,jexp+1,ic,j)=grad_M(1:2)
                do ishift=1,3
                  do i=1,nr
                    rhs(i,icol,3+3*(ic-1)+ishift)=grad_M(i+ishift+1)
                  enddo
                enddo
              enddo
              mc%DeltaM01(1:2,jexp+1,j)=M(1:2,icol)
            enddo
          enddo
          call dgemm('T','N',nr,12_8*ncol,nr,1.0d0,R%U_V,nr,rhs,nr,0.0d0,projection_rhs,nr)
          do ishift=1,12
            do icol=1,ncol
              do i=1,nr
                projection_rhs(i,icol,ishift)=R%S_inv(i)*projection_rhs(i,icol,ishift)
              enddo
            enddo
          enddo
          call dgemm('N','N',nr,12_8*ncol,nr,1.0d0,R%W_V,nr,projection_rhs,nr,0.0d0,sol,nr)
          do e=1,n_sel
            j=rays(e)
            do jexp=0,K_m
              icol=jexp+1_8+(K_m+1_8)*(e-1_8)
              do ishift=1,3
                mc%omega_shift(:,jexp+1,ishift,j)=sol(:,icol,ishift)
                do ic=1,3
                  mc%omega_shift_gradient(:,jexp+1,ishift,ic,j)=sol(:,icol,3+3*(ic-1)+ishift)
                enddo
              enddo
            enddo
          enddo
        else   ! p=3: P_mom becomes DeltaP, one solve for all K_m+7 shifts
          allocate(Gauss_P(n_P,ncol),rhs(nr,ncol,K_m+7),projection_rhs(nr,ncol,K_m+7),sol(nr,ncol,K_m+7))
          call dgemm('T','N',n_P,ncol,nr,1.0d0,R%Lj,nr,kernel_P,nr,0.0d0,Gauss_P,n_P)
          do icol=1,ncol
            do i=1,n_P
              P_mom(i,icol)=P_mom(i,icol)-Gauss_P(i,icol)
            enddo
          enddo
          do ishift=0,K_m+6
            do icol=1,ncol
              rhs(:,icol,ishift+1)=P_mom(ishift+1:ishift+nr,icol)
            enddo
          enddo
          call dgemm('T','N',nr,(K_m+7_8)*ncol,nr,1.0d0,R%U_V,nr,rhs,nr,0.0d0,projection_rhs,nr)
          do ishift=1,K_m+7
            do icol=1,ncol
              do i=1,nr
                projection_rhs(i,icol,ishift)=R%S_inv(i)*projection_rhs(i,icol,ishift)
              enddo
            enddo
          enddo
          call dgemm('N','N',nr,(K_m+7_8)*ncol,nr,1.0d0,R%W_V,nr,projection_rhs,nr,0.0d0,sol,nr)
          do e=1,n_sel
            j=rays(e)
            do jexp=0,K_m
              icol=jexp+1_8+(K_m+1_8)*(e-1_8)
              do ishift=1,K_m+7
                mc%omega_shift(:,jexp+1,ishift,j)=sol(:,icol,ishift)
              enddo
            enddo
          enddo
        endif
      end block
    endif
  end subroutine compute_moment_correction_off

  subroutine compute_singular_correction_S_off(iedge,geo,R,n_E,w_s,Jx,gp,mc,wk)
    ! SLP (p=0): weight_lambda += f*c_y(m+1,1)*J*g0*d_curv^m*sigma_m.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: w_s(:)
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(moment_correction), intent(in) :: mc
    type(kernel_weights), intent(inout) :: wk
    integer(8) :: i, j, e, im, col
    real(8) :: factor, Jg(nr,n_E)

    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        Jg(i,j)=factor*Jx%J(i,j)*gp%g0(j)   ! f*J*g0, times d_curv after each m
      enddo
    enddo
    do e=1,size(mc%rays,kind=8)
      j=mc%rays(e)
      do im=0,mc%K_m
        col=(mc%K_m+1_8)*(e-1_8)+im+1_8
        !GCC$ ivdep
        do i=1,nr
          wk%weight_lambda(i,j,1)=wk%weight_lambda(i,j,1)+R%c_y(im+1,1)*Jg(i,j)*mc%sigma(i,col)
          Jg(i,j)=Jg(i,j)*gp%d_curv(i,j)
        enddo
      enddo
    enddo
  end subroutine compute_singular_correction_S_off

  subroutine compute_singular_correction_D_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0)
    ! DLP (p=1): omega_m corrects JP_D, z_m corrects (h/lambda^2)*JP_Dh;
    !   weight_lambda += f*c_y(m+1,2)*g0*d_curv^m*(JP_D*omega_m+h*JP_Dh*z_m);
    !   weight_0 += f*c_y(m+1,2)*h*J0*g0*d_curv0^m*(DeltaM_m(0)-sum(z_m)), d_curv0 = -2*h*n_r0.c1/q0,
    !   c1 = (d_u^2*r_uu0+2*d_u*d_v*r_uv0+d_v^2*r_vv0)/2.
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(moment_correction), intent(in) :: mc
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:)
    integer(8) :: i, j, e, im, col
    real(8) :: factor, n_r0(3), c1j(3), JP_D, JP_Dh, JPg_D(nr,n_E), JPg_Dh(nr,n_E), Jg0(n_E), d_curv0(n_E)

    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      Jg0(j)=factor*geo%J0*gp%g0(j)         ! f*J0*g0, times d_curv0 after each m
      c1j=0.5d0*((geo%r_uu0*(d(1,j)*d(1,j))+geo%r_uv0*(2.0d0*d(1,j)*d(2,j)))+geo%r_vv0*(d(2,j)*d(2,j)))
      d_curv0(j)=-2.0d0*geo%h*((n_r0(1)*c1j(1)+n_r0(2)*c1j(2))+n_r0(3)*c1j(3))/X%q0(j)
      !GCC$ ivdep
      do i=1,nr
        JP_D=((X%V1r(i,1,j)-X%c1r(i,1,j))*Jx%N(i,1,j)+(X%V1r(i,2,j)-X%c1r(i,2,j))*Jx%N(i,2,j)) &
          +(X%V1r(i,3,j)-X%c1r(i,3,j))*Jx%N(i,3,j)
        JP_Dh=(n_r0(1)*Jx%N(i,1,j)+n_r0(2)*Jx%N(i,2,j))+n_r0(3)*Jx%N(i,3,j)
        JPg_D(i,j)=factor*JP_D*gp%g0(j)     ! f*JP_D*g0, times d_curv after each m
        JPg_Dh(i,j)=factor*JP_Dh*gp%g0(j)   ! f*JP_Dh*g0, times d_curv after each m
      enddo
    enddo
    do e=1,size(mc%rays,kind=8)
      j=mc%rays(e)
      do im=0,mc%K_m
        col=(mc%K_m+1_8)*(e-1_8)+im+1_8
        !GCC$ ivdep
        do i=1,nr
          wk%weight_lambda(i,j,1)=wk%weight_lambda(i,j,1)+R%c_y(im+1,2)*(geo%h*JPg_Dh(i,j)*mc%z(i,col)+JPg_D(i,j)*mc%omega(i,col))
          JPg_D(i,j)=JPg_D(i,j)*gp%d_curv(i,j)
          JPg_Dh(i,j)=JPg_Dh(i,j)*gp%d_curv(i,j)
        enddo
        weight_0(1)=weight_0(1)+R%c_y(im+1,2)*mc%h_scale(col)*Jg0(j)*mc%rest0(col)
        Jg0(j)=Jg0(j)*d_curv0(j)
      enddo
    enddo
  end subroutine compute_singular_correction_D_off

  subroutine compute_singular_correction_Sn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0)
    ! SLPn (p=1, ncomp=3): the same omega_m, alpha_m, z_m for lambda*F = J*g*(C0/lambda+c1r-(h/lambda^2)*n_r0)*mu:
    !   weight_lambda(:,:,ic) += f*c_y(m+1,2)*J*g0*d_curv^m*(C0(ic)*omega_m/lambda+c1r(:,ic)*omega_m-h*n_r0(ic)*z_m);
    !   weight_0(ic) += f*c_y(m+1,2)*J0*g0*d_curv0^m*(alpha_m*C0(ic)-h*n_r0(ic)*(DeltaM_m(0)-sum(z_m))).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(moment_correction), intent(in) :: mc
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:)
    integer(8) :: i, j, e, im, col, ic
    real(8) :: factor, n_r0(3), c1j(3), Jg(nr,n_E), Jg0(n_E), d_curv0(n_E), omega_i, omega_inv_i, z_i

    n_r0=geo%N0/geo%J0
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        Jg(i,j)=factor*Jx%J(i,j)*gp%g0(j)   ! f*J*g0, times d_curv after each m
      enddo
      Jg0(j)=factor*geo%J0*gp%g0(j)         ! f*J0*g0, times d_curv0 after each m
      c1j=0.5d0*((geo%r_uu0*(d(1,j)*d(1,j))+geo%r_uv0*(2.0d0*d(1,j)*d(2,j)))+geo%r_vv0*(d(2,j)*d(2,j)))
      d_curv0(j)=-2.0d0*geo%h*((n_r0(1)*c1j(1)+n_r0(2)*c1j(2))+n_r0(3)*c1j(3))/X%q0(j)
    enddo
    do e=1,size(mc%rays,kind=8)
      j=mc%rays(e)
      do im=0,mc%K_m
        col=(mc%K_m+1_8)*(e-1_8)+im+1_8
        !GCC$ ivdep
        do i=1,nr
          omega_i=mc%omega(i,col)
          z_i=mc%z(i,col)
          omega_inv_i=omega_i*R%lambda_inv(i)
          do ic=1,3
            wk%weight_lambda(i,j,ic)=wk%weight_lambda(i,j,ic)+R%c_y(im+1,2)*Jg(i,j) &
              *((X%C0(ic,j)*omega_inv_i+X%c1r(i,ic,j)*omega_i)-geo%h*n_r0(ic)*z_i)
          enddo
          Jg(i,j)=Jg(i,j)*gp%d_curv(i,j)
        enddo
        weight_0=weight_0+R%c_y(im+1,2)*Jg0(j)*(mc%alpha(col)*X%C0(:,j)-mc%h_scale(col)*mc%rest0(col)*n_r0)
        Jg0(j)=Jg0(j)*d_curv0(j)
      enddo
    enddo
  end subroutine compute_singular_correction_Sn_off

  subroutine compute_singular_correction_Dn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0,weight_u,weight_v)
    ! DLPn (p=2, the target gradient of SLPn), from mc's shifted solves, ray by ray: the endpoint terms
    ! e0, e1 from the first-jet model G plus the target gradient of the SLPn correction, weight_r += its gradient on
    ! mur; eps >= epsilon_big: the Gauss sums of weight_r; subtract_common leaves out the 1/|h|, log|h| terms that cancel
    ! over a closed fan, a target projecting onto an edge or a corner gets them back once.
    ! weight_0 += f*e0, (weight_u, weight_v) += f*e1*d(:,j).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(moment_correction), intent(in) :: mc
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, im, ic
    real(8) :: factor, n_r0(3), c1j(3), alpha
    real(8) :: eye3(3,3), nn(3,3), abs_h, edge_distance(4), sqrt_q0, C0j(3), crj(nr,3), c1rj(nr,3)
    real(8) :: e_C(3), a_geom(3), ae, an, ee(3,3), en(3,3), e0_ray(3,3), e1_ray(3,3), de_ray(3,3), log0(3,3), log1(3,3)
    real(8) :: rho(nr), s_lambda(nr), unit0(nr,3), Ec1(nr,3), Ecr(nr,3), Ec1r(nr,3), a_E1(nr), a_Er(nr), a_E(nr)
    real(8) :: y_E(nr), g_E(nr), g_Ediff(nr), r_y(nr), g_Erem(nr), w_rho(nr)
    real(8) :: r_L, u_L, v_L, one_minus_v, b35, b57, b47, b37, A3, A5, A7, K0(3,3), K1(3,3), Ke(3,3), dij
    real(8) :: M0(nr), M1(nr), Mr(nr), U0(nr), U1(nr), Ur(nr), sgn_h
    real(8) :: dc_gradient(nr,3), d_curv0_j, d_curv0_gradient(3), lambda_crj(nr,3), Wr(nr,3,3)
    real(8) :: power(nr), power_gradient(nr,3), power0, power0_gradient(3)
    real(8) :: omega_m(nr), omega_m_gradient(nr,3), alpha_gradient(3), z(nr), z_gradient(nr,3), rest0, rest0_gradient(3)
    real(8) :: scale, scale_J(nr), f(nr), f_gradient(nr,3), g_c1r(nr,3), g_gradient(nr,3,3), f0, f0_gradient(3)
    real(8) :: fd(nr), fd_gradient(nr,3), gd(nr,3), dcdiff(nr)
    real(8) :: fx(nr), fx_gradient(nr,3), gx(nr,3), gx_gradient(nr,3,3), fx0, fx0_gradient(3), v_alpha(3), v_z(3), sum_fz
    real(8) :: omega_1(nr), omega_2(nr), omega_1_gradient(nr,3), omega_2_gradient(nr,3), T_mur(nr,3), T_mur_k(nr)
    integer(8) :: ia, ik, icomp, q
    real(8) :: tangents(3,4), b1(3), b2(3), b_last(3), sides(3,2), sine, cosine, one_minus_cosine, span, int_b(3)
    real(8) :: leading(3,3), theta, theta_w(12), dual(3,2), jac_grad(2), model(3,3,3), e_t(3), d_t(2), A_t(3), ea
    real(8) :: ee_t(3,3), L1(3,3), LC(3,3), corners(2,4)
    integer(8) :: it, jm, iboundary, corner, n_zero, i_first
    logical :: subtract_common

    n_r0=geo%N0/geo%J0
    associate(epsilon=>mc%epsilon, epsilon_big=>mc%epsilon_big, K_m=>mc%K_m, omega_shift=>mc%omega_shift, &
      omega_shift_gradient=>mc%omega_shift_gradient, DeltaM01=>mc%DeltaM01, DeltaM01_gradient=>mc%DeltaM01_gradient)
      eye3=0.0d0
      eye3(1,1)=1.0d0; eye3(2,2)=1.0d0; eye3(3,3)=1.0d0
      do ik=1,3
        do ia=1,3
          nn(ia,ik)=n_r0(ia)*n_r0(ik)
        enddo
      enddo
      abs_h=abs(geo%h)
      sgn_h=sign(1.0d0,geo%h)
      edge_distance=geo%J0*abs(geo%x_E)/[norm2(geo%r_u0),norm2(geo%r_v0),norm2(geo%r_u0),norm2(geo%r_v0)]
      subtract_common=.true.
      do k=1,4
        if (geo%x_E(k)/=0.0d0 .and. edge_distance(k)<abs_h) subtract_common=.false.
      enddo
      do j=1,n_E
        factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
        sqrt_q0=sqrt(X%q0(j))
        C0j=X%C0(:,j)
        c1j=X%c1(1,:,j)
        crj=X%cr(:,:,j)
        c1rj=X%c1r(:,:,j)
        e_C=C0j/sqrt_q0
        a_geom=c1j/X%q0(j)
        ae=(e_C(1)*a_geom(1)+e_C(2)*a_geom(2))+e_C(3)*a_geom(3)
        do ik=1,3
          do ia=1,3
            ee(ia,ik)=e_C(ia)*e_C(ik)
            en(ia,ik)=e_C(ia)*n_r0(ik)+n_r0(ia)*e_C(ik)
          enddo
        enddo
        if (epsilon(j)>=epsilon_big) then   ! Gauss sums of step 7 (weight_r includes f)
          e0_ray=0.0d0
          e1_ray=0.0d0
          do icomp=1,9
            ia=1_8+mod(icomp-1_8,3_8)
            ik=1_8+(icomp-1_8)/3_8
            do i=1,nr
              e0_ray(ia,ik)=e0_ray(ia,ik)+wk%weight_r(i,j,icomp)*(R%lambda_inv(i)*R%lambda_inv(i))
              e1_ray(ia,ik)=e1_ray(ia,ik)+wk%weight_r(i,j,icomp)*R%lambda_inv(i)
            enddo
          enddo
          if (subtract_common) then
            e0_ray=e0_ray-factor*(geo%J0/X%q0(j))*((2.0d0*ee-eye3+nn)/abs_h-en/geo%h)
            log1=3.0d0*ee-eye3
            log0=(Jx%J1(j)/(geo%J0*sqrt_q0))*log1-15.0d0*ae*ee
            do ik=1,3
              do ia=1,3
                log0(ia,ik)=log0(ia,ik)+3.0d0*(e_C(ia)*a_geom(ik)+a_geom(ia)*e_C(ik)+ae*eye3(ia,ik))
              enddo
            enddo
            e0_ray=e0_ray+factor*(geo%J0/X%q0(j))*log(abs_h)*log0
            e1_ray=e1_ray+factor*(geo%J0/(X%q0(j)*sqrt_q0))*log(abs_h)*log1
          endif
        else
          ! first-jet model G: exact integral (K0, K1, Ke) plus the Gauss sum of F-G
          do i=1,nr
            rho(i)=hypot(R%lambda(i)*sqrt_q0,geo%h)
            s_lambda(i)=R%lambda(i)/rho(i)
            do ic=1,3
              unit0(i,ic)=s_lambda(i)*C0j(ic)-(geo%h/rho(i))*n_r0(ic)
              Ec1(i,ic)=(R%lambda(i)*s_lambda(i))*c1j(ic)
              Ecr(i,ic)=(R%lambda(i)*R%lambda(i)*s_lambda(i))*crj(i,ic)
              Ec1r(i,ic)=Ec1(i,ic)+Ecr(i,ic)
            enddo
            a_E1(i)=2.0d0*((unit0(i,1)*Ec1(i,1)+unit0(i,2)*Ec1(i,2))+unit0(i,3)*Ec1(i,3))
            a_Er(i)=2.0d0*((unit0(i,1)*Ecr(i,1)+unit0(i,2)*Ecr(i,2))+unit0(i,3)*Ecr(i,3)) &
              +((Ec1r(i,1)**2+Ec1r(i,2)**2)+Ec1r(i,3)**2)
            a_E(i)=a_E1(i)+a_Er(i)
            y_E(i)=sqrt(1.0d0+a_E(i))
            g_E(i)=1.0d0/y_E(i)**5
            g_Ediff(i)=-a_E(i)*((((y_E(i)+1.0d0)*y_E(i)+1.0d0)*y_E(i)+1.0d0)*y_E(i)+1.0d0)/(y_E(i)**5*(y_E(i)+1.0d0))
            r_y(i)=(((((5.0d0*y_E(i)+10.0d0)*y_E(i)+8.0d0)*y_E(i)+6.0d0)*y_E(i)+4.0d0)*y_E(i)+2.0d0) &
              /(2.0d0*y_E(i)**5*(y_E(i)+1.0d0)**2)
            g_Erem(i)=-2.5d0*a_Er(i)+a_E(i)**2*r_y(i)
            w_rho(i)=R%w(i)*s_lambda(i)/rho(i)**2
          enddo
          an=(n_r0(1)*a_geom(1)+n_r0(2)*a_geom(2))+n_r0(3)*a_geom(3)
          r_L=hypot(sqrt_q0,geo%h)
          u_L=sqrt_q0/r_L
          v_L=abs_h/r_L
          one_minus_v=u_L**2/(1.0d0+v_L)
          b35=sgn_h*one_minus_v**2*(v_L+2.0d0)/3.0d0
          b57=sgn_h*one_minus_v**3*(3.0d0*v_L**2+9.0d0*v_L+8.0d0)/15.0d0
          b47=u_L**5/5.0d0
          b37=sgn_h*one_minus_v**2*(2.0d0+4.0d0*v_L+6.0d0*v_L**2+3.0d0*v_L**3)/15.0d0
          if (subtract_common) then
            A3=log(sqrt_q0+r_L)-u_L
            K0=(eye3-3.0d0*ee+v_L**2*(ee-nn))/r_L+(geo%h/r_L**2)*(1.0d0+u_L+u_L**2)/(1.0d0+u_L)*en
          else   ! eps < epsilon_big gives u_L > 0.89, outside the small-u series range u < 0.125
            A3=asinh(sqrt_q0/abs_h)-u_L
            K0=(one_minus_v**2*(v_L+2.0d0)*ee-one_minus_v*(eye3-nn)+v_L*u_L**2*nn-sgn_h*u_L**3*en)/abs_h
          endif
          A5=A3-u_L**3/3.0d0
          A7=A5-u_L**5/5.0d0
          K1=3.0d0*A5*ee-A3*eye3-3.0d0*b35*en+u_L**3*nn
          do ik=1,3
            do ia=1,3
              Ke(ia,ik)=3.0d0*(A5*(e_C(ia)*a_geom(ik)+a_geom(ia)*e_C(ik)+ae*eye3(ia,ik)) &
                -b35*(n_r0(ia)*a_geom(ik)+a_geom(ia)*n_r0(ik)+an*eye3(ia,ik)))
            enddo
          enddo
          Ke=Ke-15.0d0*((ae*A7-an*b57)*ee-(ae*b57-an*b47)*en+(ae*b47-an*b37)*nn)
          e0_ray=(geo%J0/X%q0(j))*(K0+(Jx%J1(j)/(geo%J0*sqrt_q0))*K1+Ke)
          e1_ray=(geo%J0/(X%q0(j)*sqrt_q0))*K1
          do ik=1,3
            do ia=1,3
              dij=merge(1.0d0,0.0d0,ia==ik)
              M0=3.0d0*unit0(:,ia)*unit0(:,ik)-dij
              M1=3.0d0*(unit0(:,ia)*Ec1(:,ik)+Ec1(:,ia)*unit0(:,ik))-dij*a_E1
              Mr=3.0d0*(unit0(:,ia)*Ecr(:,ik)+Ecr(:,ia)*unit0(:,ik)+Ec1r(:,ia)*Ec1r(:,ik))-dij*a_Er
              U0=geo%J0*M0
              U1=R%lambda*Jx%J1(j)*M0+geo%J0*M1
              Ur=R%lambda**2*Jx%Jr(:,j)*M0+R%lambda*Jx%J1r(:,j)*M1+Jx%J(:,j)*Mr
              e0_ray(ia,ik)=e0_ray(ia,ik)+sum(w_rho*(g_E*Ur+g_Ediff*U1+g_Erem*U0))
              e1_ray(ia,ik)=e1_ray(ia,ik)+sum((w_rho*R%lambda)*(g_E*(U1+Ur)+g_Ediff*U0))
            enddo
          enddo
          ! Target gradient of the SLPn correction, term by term, including its contribution to weight_r.
          dc_gradient=-2.0d0*c1rj/X%q0(j)
          d_curv0_j=-2.0d0*geo%h*((n_r0(1)*c1j(1)+n_r0(2)*c1j(2))+n_r0(3)*c1j(3))/X%q0(j)
          d_curv0_gradient=-2.0d0*c1j/X%q0(j)
          do ic=1,3
            lambda_crj(:,ic)=R%lambda*crj(:,ic)
          enddo
          Wr=0.0d0
          power=1.0d0
          power_gradient=0.0d0
          power0=1.0d0
          power0_gradient=0.0d0
          do im=0,K_m
            omega_m=omega_shift(:,im+1,1,j)
            omega_m_gradient=omega_shift_gradient(:,im+1,1,:,j)
            alpha=DeltaM01(2,im+1,j)-sum(omega_m*R%lambda_inv)
            do ic=1,3
              alpha_gradient(ic)=DeltaM01_gradient(2,im+1,ic,j)-sum(omega_m_gradient(:,ic)*R%lambda_inv)
            enddo
            z=omega_m*R%lambda_inv**2+alpha*R%ell0_over_lambda
            do ic=1,3
              z_gradient(:,ic)=omega_m_gradient(:,ic)*R%lambda_inv**2+R%ell0_over_lambda*alpha_gradient(ic)
            enddo
            rest0=DeltaM01(1,im+1,j)-sum(z)   ! the rest of the lambda^-2 moment, at the endpoint
            do ic=1,3
              rest0_gradient(ic)=DeltaM01_gradient(1,im+1,ic,j)-sum(z_gradient(:,ic))
            enddo
            scale=R%c_y(im+1,2)/(X%q0(j)*sqrt_q0)
            f=scale*Jx%J(:,j)*power
            do ic=1,3
              f_gradient(:,ic)=scale*Jx%J(:,j)*power_gradient(:,ic)
              g_c1r(:,ic)=f*c1rj(:,ic)
            enddo
            do ik=1,3
              do ia=1,3
                g_gradient(:,ia,ik)=c1rj(:,ia)*f_gradient(:,ik)
              enddo
            enddo
            f0=scale*geo%J0*power0
            f0_gradient=scale*geo%J0*power0_gradient
            fd=f
            fd_gradient=f_gradient
            gd=g_c1r
            if (im==0) then   ! m = 0, 1: the difference from G
              fd=scale*R%lambda**2*Jx%Jr(:,j)
              do ic=1,3
                gd(:,ic)=scale*(R%lambda*Jx%J1r(:,j)*c1j(ic)+Jx%J(:,j)*lambda_crj(:,ic))
              enddo
              f0=0.0d0
              f0_gradient=0.0d0
            elseif (im==1) then
              dcdiff=(R%lambda**2*((c1rj(:,1)**2+c1rj(:,2)**2)+c1rj(:,3)**2) &
                +2.0d0*R%lambda*((lambda_crj(:,1)*C0j(1)+lambda_crj(:,2)*C0j(2))+lambda_crj(:,3)*C0j(3)) &
                -2.0d0*geo%h*((lambda_crj(:,1)*n_r0(1)+lambda_crj(:,2)*n_r0(2))+lambda_crj(:,3)*n_r0(3)))/X%q0(j)
              fd=scale*(R%lambda*Jx%J1r(:,j)*gp%d_curv(:,j)+geo%J0*dcdiff)
              do ic=1,3
                fd_gradient(:,ic)=scale*(R%lambda*Jx%J1r(:,j)*dc_gradient(:,ic)-2.0d0*geo%J0*lambda_crj(:,ic)/X%q0(j))
              enddo
              f0=0.0d0
              f0_gradient=0.0d0
            endif
            do q=0,1   ! q = 0: on mu0 (e0); q = 1: on mu1 (e1), one more power of lambda and no endpoint
              if (q==0_8) then
                fx=fd
                fx_gradient=fd_gradient
                gx=gd
                gx_gradient=g_gradient
                fx0=f0
                fx0_gradient=f0_gradient
              else
                fx=R%lambda*f
                if (im==0) fx=scale*R%lambda**2*Jx%J1r(:,j)
                do ic=1,3
                  fx_gradient(:,ic)=R%lambda*f_gradient(:,ic)
                  gx(:,ic)=R%lambda*g_c1r(:,ic)
                enddo
                do ik=1,3
                  do ia=1,3
                    gx_gradient(:,ia,ik)=R%lambda*g_gradient(:,ia,ik)
                  enddo
                enddo
                fx0=0.0d0
                fx0_gradient=0.0d0
              endif
              do ic=1,3
                v_alpha(ic)=sum(fx_gradient(:,ic)*(omega_m*R%lambda_inv)+fx*(omega_m_gradient(:,ic)*R%lambda_inv)) &
                  +fx0_gradient(ic)*alpha+fx0*alpha_gradient(ic)
                v_z(ic)=sum(fx_gradient(:,ic)*z+fx*z_gradient(:,ic))+fx0_gradient(ic)*rest0+fx0*rest0_gradient(ic)
              enddo
              sum_fz=sum(fx*z)+fx0*rest0
              do ik=1,3
                do ia=1,3
                  de_ray(ia,ik)=(C0j(ia)*v_alpha(ik)-(geo%h*n_r0(ia))*v_z(ik)-eye3(ia,ik)*sum_fz) &
                    +sum(gx_gradient(:,ia,ik)*omega_m+gx(:,ia)*omega_m_gradient(:,ik))
                enddo
              enddo
              if (q==0_8) then
                e0_ray=e0_ray+de_ray
              else
                e1_ray=e1_ray+de_ray
              endif
            enddo
            ! on mur: lambda^2*(C0/lambda+c1r-h*n_r0/lambda^2) uses the solves shifted by 1, 2, 0
            omega_1=omega_shift(:,im+1,2,j)
            omega_2=omega_shift(:,im+1,3,j)
            omega_1_gradient=omega_shift_gradient(:,im+1,2,:,j)
            omega_2_gradient=omega_shift_gradient(:,im+1,3,:,j)
            do ic=1,3
              T_mur(:,ic)=(omega_1*C0j(ic)+omega_2*c1rj(:,ic))-geo%h*omega_m*n_r0(ic)
            enddo
            scale_J=R%c_y(im+1,2)*Jx%J(:,j)/(X%q0(j)*sqrt_q0)
            do ik=1,3
              do ia=1,3
                T_mur_k=((C0j(ia)*omega_1_gradient(:,ik)+c1rj(:,ia)*omega_2_gradient(:,ik)) &
                  -merge(1.0d0,0.0d0,ia==ik)*omega_m)-geo%h*n_r0(ia)*omega_m_gradient(:,ik)
                Wr(:,ia,ik)=Wr(:,ia,ik)+(scale_J*power_gradient(:,ik))*T_mur(:,ia)+(scale_J*power)*T_mur_k
              enddo
            enddo
            do ic=1,3
              power_gradient(:,ic)=power_gradient(:,ic)*gp%d_curv(:,j)+power*dc_gradient(:,ic)
            enddo
            power0_gradient=power0_gradient*d_curv0_j+power0*d_curv0_gradient
            power=power*gp%d_curv(:,j)
            power0=power0*d_curv0_j
          enddo
          e0_ray=factor*e0_ray
          e1_ray=factor*e1_ray
          do ik=1,3
            do ia=1,3
              wk%weight_r(:,j,ia+3*(ik-1))=wk%weight_r(:,j,ia+3*(ik-1))+factor*Wr(:,ia,ik)
            enddo
          enddo
        endif
        weight_0=weight_0+reshape(e0_ray,[9])
        weight_u=weight_u+reshape(e1_ray,[9])*d(1,j)
        weight_v=weight_v+reshape(e1_ray,[9])*d(2,j)
      enddo

      ! boundary completion: a target projecting onto an edge or a corner gets the 1/|h|, log|h| terms of the missing fan
      i_first=0_8
      n_zero=0_8
      do k=4,1,-1
        if (geo%x_E(k)/=0.0d0) i_first=k
        if (geo%x_E(k)==0.0d0) then
          n_zero=n_zero+1_8
          iboundary=k
        endif
      enddo
      if (subtract_common .and. iedge==i_first .and. n_zero>0_8 .and. all(abs(geo%t0)<=1.0d0)) then
        if (n_zero==1_8) then
          tangents(:,1)=geo%r_u0
          tangents(:,2)=geo%r_v0
          tangents(:,3)=-geo%r_u0
          tangents(:,4)=-geo%r_v0
          b1=tangents(:,iboundary)/norm2(tangents(:,iboundary))
          span=pi
          sine=0.0d0
          cosine=-1.0d0
          one_minus_cosine=2.0d0
        else
          corners=reshape([-1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,1.0d0],[2,4])
          corner=0_8
          do k=4,1,-1
            if (corners(1,k)==geo%t0(1) .and. corners(2,k)==geo%t0(2)) corner=k
          enddo
          sides(:,1)=geo%r_u0*(corners(1,mod(corner,4_8)+1)-geo%t0(1))+geo%r_v0*(corners(2,mod(corner,4_8)+1)-geo%t0(2))
          sides(:,2)=geo%r_u0*(corners(1,mod(corner-2_8+4_8,4_8)+1)-geo%t0(1)) &
            +geo%r_v0*(corners(2,mod(corner-2_8+4_8,4_8)+1)-geo%t0(2))
          b1=sides(:,1)/norm2(sides(:,1))
          b_last=sides(:,2)/norm2(sides(:,2))
          sine=norm2([b1(2)*b_last(3)-b1(3)*b_last(2),b1(3)*b_last(1)-b1(1)*b_last(3),b1(1)*b_last(2)-b1(2)*b_last(1)])
          cosine=(b1(1)*b_last(1)+b1(2)*b_last(2))+b1(3)*b_last(3)
          span=atan2(sine,cosine)
          one_minus_cosine=2.0d0*sin(span/2.0d0)**2
        endif
        b2=[n_r0(2)*b1(3)-n_r0(3)*b1(2),n_r0(3)*b1(1)-n_r0(1)*b1(3),n_r0(1)*b1(2)-n_r0(2)*b1(1)]
        int_b=sine*b1+one_minus_cosine*b2
        do ik=1,3
          do ia=1,3
            leading(ia,ik)=sine*cosine*(b1(ia)*b1(ik)-b2(ia)*b2(ik))+sine**2*(b1(ia)*b2(ik)+b2(ia)*b1(ik))
            leading(ia,ik)=leading(ia,ik)/abs_h-(int_b(ia)*n_r0(ik)+n_r0(ia)*int_b(ik))/geo%h
          enddo
        enddo
        do it=1,12   ! exact for the degree-5 angular polynomial of the log terms
          theta=2.0d0*pi*real(it-1,8)/12.0d0
          theta_w(it)=span
          do jm=1,5
            theta_w(it)=theta_w(it)+2.0d0*(sin(real(jm,8)*span)*cos(real(jm,8)*theta) &
              +2.0d0*sin(real(jm,8)*span/2.0d0)**2*sin(real(jm,8)*theta))/real(jm,8)
          enddo
          theta_w(it)=theta_w(it)/12.0d0
        enddo
        dual(:,1)=[geo%r_v0(2)*n_r0(3)-geo%r_v0(3)*n_r0(2),geo%r_v0(3)*n_r0(1)-geo%r_v0(1)*n_r0(3), &
          geo%r_v0(1)*n_r0(2)-geo%r_v0(2)*n_r0(1)]/geo%J0
        dual(:,2)=[n_r0(2)*geo%r_u0(3)-n_r0(3)*geo%r_u0(2),n_r0(3)*geo%r_u0(1)-n_r0(1)*geo%r_u0(3), &
          n_r0(1)*geo%r_u0(2)-n_r0(2)*geo%r_u0(1)]/geo%J0
        jac_grad(1)=sum(dual(:,1)*geo%r_uu0)+sum(dual(:,2)*geo%r_uv0)
        jac_grad(2)=sum(dual(:,1)*geo%r_uv0)+sum(dual(:,2)*geo%r_vv0)
        model=0.0d0
        do it=1,12
          theta=2.0d0*pi*real(it-1,8)/12.0d0
          e_t=cos(theta)*b1+sin(theta)*b2
          d_t(1)=sum(dual(:,1)*e_t)
          d_t(2)=sum(dual(:,2)*e_t)
          A_t=(geo%r_uu0*(d_t(1)**2)+geo%r_uv0*(2.0d0*d_t(1)*d_t(2))+geo%r_vv0*(d_t(2)**2))/2.0d0
          ea=(e_t(1)*A_t(1)+e_t(2)*A_t(2))+e_t(3)*A_t(3)
          do ik=1,3
            do ia=1,3
              ee_t(ia,ik)=e_t(ia)*e_t(ik)
              L1(ia,ik)=3.0d0*ee_t(ia,ik)-eye3(ia,ik)
              LC(ia,ik)=3.0d0*(e_t(ia)*A_t(ik)+A_t(ia)*e_t(ik)+ea*eye3(ia,ik))-15.0d0*ea*ee_t(ia,ik)
            enddo
          enddo
          model(:,:,1)=model(:,:,1)+theta_w(it)*((jac_grad(1)*d_t(1)+jac_grad(2)*d_t(2))*L1+LC)
          model(:,:,2)=model(:,:,2)+theta_w(it)*d_t(1)*L1
          model(:,:,3)=model(:,:,3)+theta_w(it)*d_t(2)*L1
        enddo
        model=(-log(abs_h)/(2.0d0*pi))*model
        model(:,:,1)=model(:,:,1)+leading/(2.0d0*pi)
        weight_0=weight_0+reshape(model(:,:,1),[9])
        weight_u=weight_u+reshape(model(:,:,2),[9])
        weight_v=weight_v+reshape(model(:,:,3),[9])
      endif
    end associate
  end subroutine compute_singular_correction_Dn_off

  subroutine compute_singular_correction_Dnn_off(iedge,geo,R,n_E,d,w_s,X,Jx,gp,mc,wk,weight_0,weight_u,weight_v,tnx)
    ! p=3: constant and linear models on mu0, the constant model on mu1, the Gauss sum of F-G, and shifted
    ! moment corrections. Terms h^a*lambda^k use h^a*mc%omega_shift(:,m+1,k+1,j). The 1/|h| and log|h| terms
    ! that cancel over a closed fan are omitted ray by ray; an edge or corner gets its missing fan back once.
    ! The leading boundary term is integrated in closed form, the degree-5 log terms by the exact 12-point rule.
    ! f = x_E*w_s/(4*pi): weight_0 += f*e0, (weight_u,weight_v) += f*e1*d(:,j).
    integer(8), intent(in) :: iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:), tnx(3)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(moment_correction), intent(in) :: mc
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i,j,k,ic,jc,ix,iy,ell,im,ib,q,k_extra,k_shift,h_power,imodel,l
    integer(8), parameter :: model_m(4)=[0_8,0_8,1_8,1_8], model_k(4)=[4_8,5_8,5_8,4_8]
    integer(8), parameter :: model_a(4)=[0_8,0_8,0_8,1_8]
    logical :: subtract_common,strip_leading(4),strip_log(4)
    real(8) :: factor,n_r0(3),h,abs_h,sign_h,eye3(3,3),nn(3,3),proj_t(3,3),tnx_n,tnx_t(3),edge_distance(4)
    real(8) :: q0,sqrt_q0,C0j(3),c1j(3),crj(nr,3),c1rj(nr,3),N1j(3),Nrj(nr,3),a(nr,3),N0(3),N1r(nr,3),N(nr,3)
    real(8) :: B_w0(3),B_w1(3),B_wr(nr,3),B_w1r(nr,3),B_w(nr,3),ct0(2),ct1(2),ctr(nr,2),ct1r(nr,2),ct(nr,2)
    real(8) :: Nc0(2),Nc1(2),Ncr(nr,2),Nc1r(nr,2),Nc(nr,2),Nt0,Nt1,Ntr(nr),Nt1r(nr),Nt(nr)
    real(8) :: B_sq0(5),B_sq1(5),B_sqr(nr,5),B_sq1r(nr,5),B_sq(nr,5)
    real(8) :: ctNc0(3),ctNc1(3),ctNcr(nr,3),ctNc1r(nr,3),ctNc(nr,3)
    real(8) :: JP_w0(5,9),JP_w1(5,9),JP_wr(nr,9,5),dij,U0(3),U1(3),Ur(nr,3),P0(2),P1(2),Pr(nr,2)
    real(8) :: L0(2),L1(2),Lr(nr,2),V0(3),V1(3),Vr(nr,3),e0_ray(9),e1_ray(9),e_C(3),ee(3,3),tnx_e_C,v_C(3)
    real(8) :: leading_C(3,3),log0(9),log1(9),r_eps,u_eps,v_eps,M_model(5,4),moment,exponent,power,coefficient,tail,term
    real(8) :: beta_model,C0_c1,n_c1,rho(nr),s_lambda(nr),t_lambda(nr),unit0(nr,3),Ec1(nr,3),Ecr(nr,3),Ec1r(nr,3)
    real(8) :: a_E1(nr),a_Er(nr),a_E(nr),y_E(nr),g_E(nr),g_Ediff(nr),r_y(nr),g_Erem(nr),w_rho(nr),power_E(nr)
    real(8) :: JP_E0(nr,9),JP_E1(nr,9),JP_Er(nr,9),beta1r(nr),betar(nr),gamma_h(nr),gamma_h1r(nr),T_mur(nr,9)
    real(8) :: scale,JP_wr_ell(nr,9),JP_w1r_ell(nr,9),JP_w_ell(nr,9),d_curv_term(nr),smooth(nr,9),de_ray(nr,9)

    n_r0=geo%N0/geo%J0
    h=geo%h; abs_h=abs(h); sign_h=sign(1.0d0,h)
    eye3=0.0d0
    eye3(1,1)=1.0d0; eye3(2,2)=1.0d0; eye3(3,3)=1.0d0
    do jc=1,3
      do ic=1,3
        nn(ic,jc)=n_r0(ic)*n_r0(jc)
      enddo
    enddo
    proj_t=eye3-nn
    tnx_n=dot_product(n_r0,tnx); tnx_t=tnx-tnx_n*n_r0
    edge_distance=geo%J0*abs(geo%x_E)/[norm2(geo%r_u0),norm2(geo%r_v0),norm2(geo%r_u0),norm2(geo%r_v0)]
    subtract_common=.true.
    do k=1,4
      if (geo%x_E(k)/=0.0d0 .and. edge_distance(k)<abs_h) subtract_common=.false.
    enddo
    associate(lambda=>R%lambda,lambda_inv=>R%lambda_inv)
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(4.0d0*pi)
      q0=X%q0(j); sqrt_q0=sqrt(q0)
      C0j=X%C0(:,j); c1j=X%c1(1,:,j)
      crj=X%cr(:,:,j); c1rj=X%c1r(:,:,j)
      N1j=Jx%N1(1,:,j); Nrj=Jx%Nr(:,:,j)
      ! Coefficients of w=h/lambda; suffixes 0,1,r denote the lambda expansion.
      ! B_w=c.c, ct=c.tnx, Nc=N.c, Nt=N.tnx, with c=a-w*n_r0.
      N0=Jx%N0(1,:)
      do ic=1,3
        a(:,ic)=C0j(ic)+lambda*c1rj(:,ic)
        N1r(:,ic)=N1j(ic)+lambda*Nrj(:,ic)
        N(:,ic)=N0(ic)+lambda*N1r(:,ic)
      enddo
      B_w0=[sum(C0j**2),-2.0d0*dot_product(C0j,n_r0),sum(n_r0**2)]
      B_w1=[2.0d0*dot_product(C0j,c1j),-2.0d0*dot_product(c1j,n_r0),0.0d0]
      B_wr(:,1)=2.0d0*matmul(crj,C0j)+sum(c1rj**2,dim=2)
      B_wr(:,2)=-2.0d0*matmul(crj,n_r0); B_wr(:,3)=0.0d0
      do ix=1,3
        B_w1r(:,ix)=B_w1(ix)+lambda*B_wr(:,ix)
        B_w(:,ix)=B_w0(ix)+lambda*B_w1r(:,ix)
      enddo
      ct0=[dot_product(C0j,tnx),-dot_product(n_r0,tnx)]
      ct1=[dot_product(c1j,tnx),0.0d0]
      ctr(:,1)=matmul(crj,tnx); ctr(:,2)=0.0d0
      Nc0=[dot_product(N0,C0j),-dot_product(N0,n_r0)]
      Nc1=[dot_product(N1j,C0j)+dot_product(N0,c1j),-dot_product(N1j,n_r0)]
      Ncr(:,1)=matmul(Nrj,C0j)+matmul(N1r,c1j)+sum(N*crj,dim=2)
      Ncr(:,2)=-matmul(Nrj,n_r0)
      do ix=1,2
        ct1r(:,ix)=ct1(ix)+lambda*ctr(:,ix)
        ct(:,ix)=ct0(ix)+lambda*ct1r(:,ix)
        Nc1r(:,ix)=Nc1(ix)+lambda*Ncr(:,ix)
        Nc(:,ix)=Nc0(ix)+lambda*Nc1r(:,ix)
      enddo
      Nt0=dot_product(N0,tnx); Nt1=dot_product(N1j,tnx); Ntr=matmul(Nrj,tnx)
      Nt1r=Nt1+lambda*Ntr; Nt=Nt0+lambda*Nt1r
      ! B_w^2 and (c.tnx)*(N.c).
      B_sq0=0.0d0; B_sq1=0.0d0; B_sqr=0.0d0
      do ix=1,3
        do iy=1,3
          ell=ix+iy-1
          B_sq0(ell)=B_sq0(ell)+B_w0(ix)*B_w0(iy)
          B_sq1(ell)=B_sq1(ell)+B_w0(ix)*B_w1(iy)+B_w1(ix)*B_w0(iy)
          B_sqr(:,ell)=B_sqr(:,ell)+B_wr(:,ix)*B_w0(iy)+B_w1r(:,ix)*B_w1(iy)+B_w(:,ix)*B_wr(:,iy)
        enddo
      enddo
      do ell=1,5
        B_sq1r(:,ell)=B_sq1(ell)+lambda*B_sqr(:,ell)
        B_sq(:,ell)=B_sq0(ell)+lambda*B_sq1r(:,ell)
      enddo
      ctNc0=0.0d0; ctNc1=0.0d0; ctNcr=0.0d0
      do ix=1,2
        do iy=1,2
          ell=ix+iy-1
          ctNc0(ell)=ctNc0(ell)+ct0(ix)*Nc0(iy)
          ctNc1(ell)=ctNc1(ell)+ct0(ix)*Nc1(iy)+ct1(ix)*Nc0(iy)
          ctNcr(:,ell)=ctNcr(:,ell)+ctr(:,ix)*Nc0(iy)+ct1r(:,ix)*Nc1(iy)+ct(:,ix)*Ncr(:,iy)
        enddo
      enddo
      do ell=1,3
        ctNc1r(:,ell)=ctNc1(ell)+lambda*ctNcr(:,ell)
        ctNc(:,ell)=ctNc0(ell)+lambda*ctNc1r(:,ell)
      enddo
      JP_w0=0.0d0; JP_w1=0.0d0; JP_wr=0.0d0
      do jc=1,3
        do ic=1,3
          k=ic+3*(jc-1); dij=merge(1.0d0,0.0d0,ic==jc)
          ! U=c_i*c_j, P=N_i*c_j, L=tnx_j*c_i+delta_ij*(c.tnx).
          U0=[C0j(ic)*C0j(jc),-C0j(ic)*n_r0(jc)-n_r0(ic)*C0j(jc),n_r0(ic)*n_r0(jc)]
          U1=[c1j(ic)*C0j(jc)+C0j(ic)*c1j(jc),-c1j(ic)*n_r0(jc)-n_r0(ic)*c1j(jc),0.0d0]
          Ur(:,1)=crj(:,ic)*C0j(jc)+c1rj(:,ic)*c1j(jc)+a(:,ic)*crj(:,jc)
          Ur(:,2)=-crj(:,ic)*n_r0(jc)-n_r0(ic)*crj(:,jc); Ur(:,3)=0.0d0
          P0=[N0(ic)*C0j(jc),-N0(ic)*n_r0(jc)]
          P1=[N1j(ic)*C0j(jc)+N0(ic)*c1j(jc),-N1j(ic)*n_r0(jc)]
          Pr(:,1)=Nrj(:,ic)*C0j(jc)+N1r(:,ic)*c1j(jc)+N(:,ic)*crj(:,jc)
          Pr(:,2)=-Nrj(:,ic)*n_r0(jc)
          L0=tnx(jc)*[C0j(ic),-n_r0(ic)]+dij*ct0
          L1=tnx(jc)*[c1j(ic),0.0d0]+dij*ct1
          Lr(:,1)=tnx(jc)*crj(:,ic)+dij*ctr(:,1); Lr(:,2)=dij*ctr(:,2)
          ! V=(N.tnx)*U+(c.tnx)*P+(N.c)*L.
          V0=Nt0*U0; V1=Nt0*U1+Nt1*U0
          do ell=1,3
            Vr(:,ell)=Ntr*U0(ell)+Nt1r*U1(ell)+Nt*Ur(:,ell)
          enddo
          do ix=1,2
            do iy=1,2
              ell=ix+iy-1
              V0(ell)=V0(ell)+ct0(ix)*P0(iy)+Nc0(ix)*L0(iy)
              V1(ell)=V1(ell)+ct0(ix)*P1(iy)+ct1(ix)*P0(iy)+Nc0(ix)*L1(iy)+Nc1(ix)*L0(iy)
              Vr(:,ell)=Vr(:,ell)+ctr(:,ix)*P0(iy)+ct1r(:,ix)*P1(iy)+ct(:,ix)*Pr(:,iy) &
                +Ncr(:,ix)*L0(iy)+Nc1r(:,ix)*L1(iy)+Nc(:,ix)*Lr(:,iy)
            enddo
          enddo
          ! JP_Dn=2*tnx_i*B_w^2*N_j+3*B_w*V-30*(c.tnx)*(N.c)*U.
          JP_w0(:,k)=2.0d0*tnx(ic)*N0(jc)*B_sq0
          JP_w1(:,k)=2.0d0*tnx(ic)*(N0(jc)*B_sq1+N1j(jc)*B_sq0)
          do ell=1,5
            JP_wr(:,k,ell)=2.0d0*tnx(ic)*(B_sqr(:,ell)*N0(jc)+B_sq1r(:,ell)*N1j(jc)+B_sq(:,ell)*Nrj(:,jc))
          enddo
          do ix=1,3
            do iy=1,3
              ell=ix+iy-1
              JP_w0(ell,k)=JP_w0(ell,k)+3.0d0*B_w0(ix)*V0(iy)-30.0d0*ctNc0(ix)*U0(iy)
              JP_w1(ell,k)=JP_w1(ell,k)+3.0d0*(B_w0(ix)*V1(iy)+B_w1(ix)*V0(iy)) &
                -30.0d0*(ctNc0(ix)*U1(iy)+ctNc1(ix)*U0(iy))
              JP_wr(:,k,ell)=JP_wr(:,k,ell)+3.0d0*(B_wr(:,ix)*V0(iy)+B_w1r(:,ix)*V1(iy)+B_w(:,ix)*Vr(:,iy)) &
                -30.0d0*(ctNcr(:,ix)*U0(iy)+ctNc1r(:,ix)*U1(iy)+ctNc(:,ix)*Ur(:,iy))
            enddo
          enddo
        enddo
      enddo
      if (mc%epsilon(j)>=mc%epsilon_big) then   ! Gauss sums of step 7 (weight_r includes f)
        do k=1,9
          e0_ray(k)=sum(wk%weight_r(:,j,k)*lambda_inv**2)
          e1_ray(k)=sum(wk%weight_r(:,j,k)*lambda_inv)
        enddo
        if (subtract_common) then
          e_C=C0j/sqrt_q0; tnx_e_C=dot_product(tnx_t,e_C); v_C=tnx_t-2.0d0*tnx_e_C*e_C
          do jc=1,3
            do ic=1,3
              ee(ic,jc)=e_C(ic)*e_C(jc)
              leading_C(ic,jc)=tnx_n*(proj_t(ic,jc)-2.0d0*ee(ic,jc)) &
                +2.0d0*sign_h*tnx_n*(e_C(ic)*n_r0(jc)+n_r0(ic)*e_C(jc))+2.0d0*v_C(ic)*n_r0(jc)+n_r0(ic)*v_C(jc) &
                +sign_h*(-e_C(ic)*tnx_t(jc)+6.0d0*tnx_e_C*ee(ic,jc)-tnx_e_C*proj_t(ic,jc)+2.0d0*tnx_e_C*nn(ic,jc))
            enddo
          enddo
          log0=gp%g0(j)*(JP_w1(1,:)-7.0d0*dot_product(C0j,c1j)/q0*JP_w0(1,:))
          log1=gp%g0(j)*JP_w0(1,:)
          e0_ray=e0_ray-factor*(geo%J0/q0)*reshape(leading_C,[9])/abs_h+factor*log(abs_h)*log0
          e1_ray=e1_ray+factor*log(abs_h)*log1
        endif
      else
        ! Exact models: h^a*integral lambda^k*k_m. The first model omits its leading beta term when subtract_common;
        ! the second and third omit log|h|. Here k_m=lambda^(2m+1)/(lambda^2+eps^2)^(m+7/2).
        r_eps=hypot(1.0d0,mc%epsilon(j)); u_eps=1.0d0/r_eps; v_eps=mc%epsilon(j)/r_eps
        strip_leading=[subtract_common,.false.,.false.,.false.]
        strip_log=[.false.,subtract_common,subtract_common,.false.]
        M_model=0.0d0
        do ell=0,4
          do imodel=1,4
            im=model_m(imodel); k_shift=model_k(imodel)-ell; h_power=model_a(imodel)+ell
            if (k_shift==5) then   ! asinh(1/eps)=log(1+r_eps)+log(sqrt_q0)-log|h|
              if (strip_log(imodel)) then
                moment=log(1.0d0+r_eps)+0.5d0*log(q0)
              else
                moment=asinh(1.0d0/mc%epsilon(j))
              endif
              do l=0,im+2
                moment=moment-u_eps**(2*l+1)/real(2*l+1,8)
              enddo
              moment=h**h_power*moment
            else   ! lower tail of the beta integral, as a binomial series in v_eps<0.45
              exponent=real(im,8)+real(k_shift,8)/2.0d0; power=1.0d0; coefficient=1.0d0
              tail=1.0d0/real(5-k_shift,8)
              do l=1,40
                coefficient=coefficient*(real(l-1,8)-exponent)/real(l,8); power=power*v_eps*v_eps
                term=coefficient*power/real(5-k_shift+2*l,8); tail=tail+term
                if (abs(term)<=epsilon(1.0d0)*max(1.0d0,abs(tail))) exit
              enddo
              moment=-h**h_power*r_eps**(k_shift-5)*tail
              if (.not.strip_leading(imodel)) then
                beta_model=gamma(real(5-k_shift,8)/2.0d0)*gamma(real(im,8)+real(k_shift,8)/2.0d0+1.0d0) &
                  /gamma(real(im,8)+3.5d0)
                moment=moment+sign_h**h_power*q0**(real(h_power,8)/2.0d0) &
                  *mc%epsilon(j)**(h_power+k_shift-5)*beta_model/2.0d0
              endif
            endif
            M_model(ell+1,imodel)=moment
          enddo
        enddo
        C0_c1=dot_product(C0j,c1j); n_c1=dot_product(n_r0,c1j)
        e0_ray=0.0d0; e1_ray=0.0d0
        do ell=0,4
          e0_ray=e0_ray+JP_w0(ell+1,:)*M_model(ell+1,1)+JP_w1(ell+1,:)*M_model(ell+1,2) &
            -7.0d0/q0*JP_w0(ell+1,:)*(C0_c1*M_model(ell+1,3)-n_c1*M_model(ell+1,4))
          e1_ray=e1_ray+JP_w0(ell+1,:)*M_model(ell+1,2)
        enddo
        e0_ray=gp%g0(j)*e0_ray; e1_ray=gp%g0(j)*e1_ray
        ! Gauss sum of F-G at the scaled chord unit0+Ec1r=(lambda/rho)*c; the quartic numerator coefficients
        ! use s_lambda^(4-ell)*t_lambda^ell*(JP_w0,lambda*JP_w1,lambda^2*JP_wr).
        do i=1,nr
          rho(i)=hypot(lambda(i)*sqrt_q0,h)
        enddo
        s_lambda=lambda/rho; t_lambda=h/rho
        do ic=1,3
          unit0(:,ic)=s_lambda*C0j(ic)-t_lambda*n_r0(ic)
          Ec1(:,ic)=(lambda*s_lambda)*c1j(ic)
          Ecr(:,ic)=(lambda*lambda*s_lambda)*crj(:,ic)
        enddo
        Ec1r=Ec1+Ecr
        a_E1=2.0d0*sum(unit0*Ec1,dim=2)
        a_Er=2.0d0*sum(unit0*Ecr,dim=2)+sum(Ec1r**2,dim=2); a_E=a_E1+a_Er
        y_E=sqrt(1.0d0+a_E); g_E=1.0d0/y_E**7
        g_Ediff=-a_E*((((((y_E+1.0d0)*y_E+1.0d0)*y_E+1.0d0)*y_E+1.0d0)*y_E+1.0d0)*y_E+1.0d0) &
          /(y_E**7*(y_E+1.0d0))
        r_y=(((((((7.0d0*y_E+14.0d0)*y_E+12.0d0)*y_E+10.0d0)*y_E+8.0d0)*y_E+6.0d0)*y_E+4.0d0)*y_E+2.0d0) &
          /(2.0d0*y_E**7*(y_E+1.0d0)**2)
        g_Erem=-3.5d0*a_Er+a_E**2*r_y; w_rho=R%w*s_lambda/rho**2
        JP_E0=0.0d0; JP_E1=0.0d0; JP_Er=0.0d0
        do ell=0,4
          power_E=s_lambda**(4-ell)*t_lambda**ell
          do k=1,9
            JP_E0(:,k)=JP_E0(:,k)+power_E*JP_w0(ell+1,k)
            JP_E1(:,k)=JP_E1(:,k)+(power_E*lambda)*JP_w1(ell+1,k)
            JP_Er(:,k)=JP_Er(:,k)+(power_E*lambda**2)*JP_wr(:,k,ell+1)
          enddo
        enddo
        do k=1,9
          e0_ray(k)=e0_ray(k)+sum(w_rho*(g_E*JP_Er(:,k)+g_Ediff*JP_E1(:,k)+g_Erem*JP_E0(:,k)))
          e1_ray(k)=e1_ray(k)+sum((w_rho*lambda)*(g_E*(JP_E1(:,k)+JP_Er(:,k))+g_Ediff*JP_E0(:,k)))
        enddo
        ! Moment correction of F-G and the whole mur coefficient: d_curv=lambda*beta1r-h*gamma_h.
        ! h^ib*lambda^(im-ib) in d_curv^im and w^ell in JP_Dn take the shift 4-ell+q+im-ib+k_extra and h^(ell+ib).
        beta1r=(2.0d0*matmul(c1rj,C0j)+lambda*sum(c1rj**2,dim=2))/q0
        betar=(2.0d0*matmul(crj,C0j)+sum(c1rj**2,dim=2))/q0
        gamma_h=2.0d0*matmul(c1rj,n_r0)/q0; gamma_h1r=2.0d0*matmul(crj,n_r0)/q0
        T_mur=0.0d0
        do im=0,mc%K_m
          scale=R%c_y(im+1,4)*gp%g0(j)
          do ell=0,4
            JP_wr_ell=JP_wr(:,:,ell+1)
            do k=1,9
              JP_w1r_ell(:,k)=JP_w1(ell+1,k)+lambda*JP_wr_ell(:,k)
              JP_w_ell(:,k)=JP_w0(ell+1,k)+lambda*JP_w1r_ell(:,k)
            enddo
            do ib=0,im
              d_curv_term=R%binom(ib+1,im+1)*beta1r**(im-ib)*(-gamma_h)**ib
              do q=0,2
                k_extra=0_8
                if (im==0 .and. q==0) then   ! mu0: constant and linear terms are in the model
                  k_extra=2_8; smooth=JP_wr_ell
                elseif (im==0 .and. q==1) then   ! mu1: the constant term is in the model
                  k_extra=1_8; smooth=JP_w1r_ell
                elseif (im==1 .and. q==0) then   ! mu0: constant and linear terms of d_curv are in the model
                  k_extra=1_8
                  do k=1,9
                    if (ib==0) then
                      smooth(:,k)=JP_w1r_ell(:,k)*beta1r+JP_w0(ell+1,k)*betar
                    else
                      smooth(:,k)=-JP_w1r_ell(:,k)*gamma_h-JP_w0(ell+1,k)*gamma_h1r
                    endif
                  enddo
                else
                  do k=1,9
                    smooth(:,k)=JP_w_ell(:,k)*d_curv_term
                  enddo
                endif
                k_shift=4-ell+q+im-ib+k_extra
                do k=1,9
                  de_ray(:,k)=(scale*h**(ell+ib))*mc%omega_shift(:,im+1,k_shift+1,j)*smooth(:,k)
                enddo
                if (q==0) then
                  e0_ray=e0_ray+sum(de_ray,dim=1)
                elseif (q==1) then
                  e1_ray=e1_ray+sum(de_ray,dim=1)
                else
                  T_mur=T_mur+de_ray
                endif
              enddo
            enddo
          enddo
        enddo
        e0_ray=factor*e0_ray; e1_ray=factor*e1_ray
        do k=1,9
          wk%weight_r(:,j,k)=wk%weight_r(:,j,k)+factor*T_mur(:,k)
        enddo
      endif
      weight_0=weight_0+e0_ray; weight_u=weight_u+d(1,j)*e1_ray; weight_v=weight_v+d(2,j)*e1_ray
    enddo
    end associate

    ! Boundary completion: the missing fan's 1/|h| and log|h| terms at an exact edge or corner projection.
    block
      integer(8) :: i_first,n_zero,iboundary,corner,it,jm
      real(8) :: tangents(3,4),b1(3),b2(3),b_last(3),sides(3,2),corners(2,4),span,sine,cosine,one_minus_cosine
      real(8) :: int_b(3),int_ee(3,3),int_ccc,int_scc,int_css,int_sss,tnx_b1,tnx_b2,int_tee(3,3),int_ee_tnx(3),int_tnx_e
      real(8) :: leading(3,3),theta,theta_w(12),dual(3,2),N_u0(3),N_v0(3),model(3,3,3),e_t(3),d_t(2),A_t(3),N1_t(3)
      real(8) :: tnx_e_t,ea,an,tnx_A,tnx_N1,log1_t(3,3),log0_t(3,3)
      i_first=0_8; n_zero=0_8
      do k=4,1,-1
        if (geo%x_E(k)/=0.0d0) i_first=k
        if (geo%x_E(k)==0.0d0) then
          n_zero=n_zero+1_8; iboundary=k
        endif
      enddo
      if (subtract_common .and. iedge==i_first .and. n_zero>0_8 .and. all(abs(geo%t0)<=1.0d0)) then
        if (n_zero==1_8) then
          tangents(:,1)=geo%r_u0; tangents(:,2)=geo%r_v0
          tangents(:,3)=-geo%r_u0; tangents(:,4)=-geo%r_v0
          b1=tangents(:,iboundary)/norm2(tangents(:,iboundary))
          span=pi; sine=0.0d0; cosine=-1.0d0; one_minus_cosine=2.0d0
        else
          corners=reshape([-1.0d0,-1.0d0,1.0d0,-1.0d0,1.0d0,1.0d0,-1.0d0,1.0d0],[2,4])
          corner=0_8
          do k=4,1,-1
            if (all(corners(:,k)==geo%t0)) corner=k
          enddo
          sides(:,1)=geo%r_u0*(corners(1,mod(corner,4_8)+1)-geo%t0(1))+geo%r_v0*(corners(2,mod(corner,4_8)+1)-geo%t0(2))
          sides(:,2)=geo%r_u0*(corners(1,mod(corner-2_8+4_8,4_8)+1)-geo%t0(1)) &
            +geo%r_v0*(corners(2,mod(corner-2_8+4_8,4_8)+1)-geo%t0(2))
          b1=sides(:,1)/norm2(sides(:,1)); b_last=sides(:,2)/norm2(sides(:,2))
          sine=norm2([b1(2)*b_last(3)-b1(3)*b_last(2),b1(3)*b_last(1)-b1(1)*b_last(3),b1(1)*b_last(2)-b1(2)*b_last(1)])
          cosine=dot_product(b1,b_last); span=atan2(sine,cosine); one_minus_cosine=2.0d0*sin(span/2.0d0)**2
        endif
        b2=[n_r0(2)*b1(3)-n_r0(3)*b1(2),n_r0(3)*b1(1)-n_r0(1)*b1(3),n_r0(1)*b1(2)-n_r0(2)*b1(1)]
        ! Integrate leading_C over e=cos(theta)*b1+sin(theta)*b2.
        int_b=sine*b1+one_minus_cosine*b2
        do jc=1,3
          do ic=1,3
            int_ee(ic,jc)=-sine*cosine*(b1(ic)*b1(jc)-b2(ic)*b2(jc))-sine**2*(b1(ic)*b2(jc)+b2(ic)*b1(jc))
          enddo
        enddo
        int_ccc=sine-sine**3/3.0d0; int_scc=one_minus_cosine*(1.0d0+cosine+cosine**2)/3.0d0
        int_css=sine**3/3.0d0; int_sss=one_minus_cosine**2*(cosine+2.0d0)/3.0d0
        tnx_b1=dot_product(b1,tnx_t); tnx_b2=dot_product(b2,tnx_t)
        do jc=1,3
          do ic=1,3
            int_tee(ic,jc)=(tnx_b1*int_ccc+tnx_b2*int_scc)*(b1(ic)*b1(jc)) &
              +(tnx_b1*int_css+tnx_b2*int_sss)*(b2(ic)*b2(jc)) &
              +(tnx_b1*int_scc+tnx_b2*int_css)*(b1(ic)*b2(jc)+b2(ic)*b1(jc))
          enddo
        enddo
        int_ee_tnx=matmul(int_ee,tnx_t); int_tnx_e=dot_product(tnx_t,int_b)
        do jc=1,3
          do ic=1,3
            leading(ic,jc)=tnx_n*int_ee(ic,jc)+2.0d0*sign_h*tnx_n*(int_b(ic)*n_r0(jc)+n_r0(ic)*int_b(jc)) &
              +2.0d0*int_ee_tnx(ic)*n_r0(jc)+n_r0(ic)*int_ee_tnx(jc) &
              +sign_h*(-int_b(ic)*tnx_t(jc)+6.0d0*int_tee(ic,jc)-int_tnx_e*proj_t(ic,jc)+2.0d0*int_tnx_e*nn(ic,jc))
          enddo
        enddo
        do it=1,12   ! exact for the degree-5 angular polynomials of the log terms
          theta=2.0d0*pi*real(it-1,8)/12.0d0; theta_w(it)=span
          do jm=1,5
            theta_w(it)=theta_w(it)+2.0d0*(sin(real(jm,8)*span)*cos(real(jm,8)*theta) &
              +2.0d0*sin(real(jm,8)*span/2.0d0)**2*sin(real(jm,8)*theta))/real(jm,8)
          enddo
          theta_w(it)=theta_w(it)/12.0d0
        enddo
        dual(:,1)=[geo%r_v0(2)*n_r0(3)-geo%r_v0(3)*n_r0(2),geo%r_v0(3)*n_r0(1)-geo%r_v0(1)*n_r0(3), &
          geo%r_v0(1)*n_r0(2)-geo%r_v0(2)*n_r0(1)]/geo%J0
        dual(:,2)=[n_r0(2)*geo%r_u0(3)-n_r0(3)*geo%r_u0(2),n_r0(3)*geo%r_u0(1)-n_r0(1)*geo%r_u0(3), &
          n_r0(1)*geo%r_u0(2)-n_r0(2)*geo%r_u0(1)]/geo%J0
        N_u0=[geo%r_uu0(2)*geo%r_v0(3)-geo%r_uu0(3)*geo%r_v0(2), &
          geo%r_uu0(3)*geo%r_v0(1)-geo%r_uu0(1)*geo%r_v0(3),geo%r_uu0(1)*geo%r_v0(2)-geo%r_uu0(2)*geo%r_v0(1)] &
          +[geo%r_u0(2)*geo%r_uv0(3)-geo%r_u0(3)*geo%r_uv0(2),geo%r_u0(3)*geo%r_uv0(1)-geo%r_u0(1)*geo%r_uv0(3), &
          geo%r_u0(1)*geo%r_uv0(2)-geo%r_u0(2)*geo%r_uv0(1)]
        N_v0=[geo%r_uv0(2)*geo%r_v0(3)-geo%r_uv0(3)*geo%r_v0(2), &
          geo%r_uv0(3)*geo%r_v0(1)-geo%r_uv0(1)*geo%r_v0(3),geo%r_uv0(1)*geo%r_v0(2)-geo%r_uv0(2)*geo%r_v0(1)] &
          +[geo%r_u0(2)*geo%r_vv0(3)-geo%r_u0(3)*geo%r_vv0(2),geo%r_u0(3)*geo%r_vv0(1)-geo%r_u0(1)*geo%r_vv0(3), &
          geo%r_u0(1)*geo%r_vv0(2)-geo%r_u0(2)*geo%r_vv0(1)]
        model=0.0d0
        do it=1,12
          theta=2.0d0*pi*real(it-1,8)/12.0d0; e_t=cos(theta)*b1+sin(theta)*b2
          d_t=matmul(transpose(dual),e_t)
          A_t=(geo%r_uu0*d_t(1)**2+geo%r_uv0*(2.0d0*d_t(1)*d_t(2))+geo%r_vv0*d_t(2)**2)/2.0d0
          N1_t=(N_u0*d_t(1)+N_v0*d_t(2))/geo%J0
          tnx_e_t=dot_product(e_t,tnx); ea=dot_product(e_t,A_t); an=dot_product(n_r0,A_t)
          tnx_A=dot_product(tnx,A_t); tnx_N1=dot_product(tnx,N1_t)
          do jc=1,3
            do ic=1,3
              log1_t(ic,jc)=2.0d0*tnx(ic)*n_r0(jc)+3.0d0*tnx_e_t*n_r0(ic)*e_t(jc)+3.0d0*tnx_n*(e_t(ic)*e_t(jc))
              log0_t(ic,jc)=2.0d0*tnx(ic)*N1_t(jc)-6.0d0*ea*tnx(ic)*n_r0(jc)+3.0d0*tnx_e_t*N1_t(ic)*e_t(jc) &
                +3.0d0*(tnx_A-5.0d0*ea*tnx_e_t)*n_r0(ic)*e_t(jc)+3.0d0*tnx_e_t*n_r0(ic)*A_t(jc) &
                -3.0d0*an*e_t(ic)*tnx(jc)+3.0d0*tnx_n*(A_t(ic)*e_t(jc)+e_t(ic)*A_t(jc)) &
                +3.0d0*(tnx_N1-5.0d0*ea*tnx_n+10.0d0*tnx_e_t*an)*(e_t(ic)*e_t(jc))-3.0d0*tnx_e_t*an*eye3(ic,jc)
            enddo
          enddo
          model(:,:,1)=model(:,:,1)+theta_w(it)*log0_t
          model(:,:,2)=model(:,:,2)+theta_w(it)*d_t(1)*log1_t
          model(:,:,3)=model(:,:,3)+theta_w(it)*d_t(2)*log1_t
        enddo
        model=(-log(abs_h)/(4.0d0*pi))*model
        model(:,:,1)=model(:,:,1)+leading/(4.0d0*pi*abs_h)
        weight_0=weight_0+reshape(model(:,:,1),[9])
        weight_u=weight_u+reshape(model(:,:,2),[9])
        weight_v=weight_v+reshape(model(:,:,3),[9])
      endif
    end block
  end subroutine compute_singular_correction_Dnn_off

  subroutine compute_singular_correction_helm_off(zk,p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
    ! Step 8, Helmholtz off surface (h/=0): compute_singular_correction_off on the Laplace part, its weights times the even
    ! Helmholtz factor of step 7, a function of dist^2 (smooth in lambda, zero gradient at t'): SLP (p=0) A = cos(z);
    ! DLP, SLPn (p=1) A = c_e = cos(z)+z*sin(z); z = zk*dist, dist = |lambda*(C0+lambda*c1r)-h*n_r0|; A at the radial
    ! nodes, A0 = A(t') (dist = |h|) on the endpoint weight. compute_singular_correction_off starts from zero weights (wk_L),
    ! so only this edge's correction is scaled. The odd part of step 7 is smooth and needs no correction.
    ! DLPn (p=2): compute_singular_correction_Dn_helm_off.
    ! Components 1:ncomp_L real parts, ncomp_L+1:2*ncomp_L imaginary parts.
    complex(8), intent(in) :: zk
    integer(8), intent(in) :: p, order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    type(kernel_weights) :: wk_L
    integer(8) :: i, j, ic, ncomp_L
    real(8) :: n_r0(3), chord(nr,3,n_E), weight_0_L(3), weight_u_L(3), weight_v_L(3)
    complex(8) :: z(nr,n_E), z0, A(nr,n_E), A0

    n_r0=geo%N0/geo%J0
    do j=1,n_E
      do i=1,nr
        chord(i,:,j)=R%lambda(i)*(X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j))-geo%h*n_r0
        z(i,j)=zk*sqrt((chord(i,1,j)*chord(i,1,j)+chord(i,2,j)*chord(i,2,j))+chord(i,3,j)*chord(i,3,j))
      enddo
    enddo
    z0=zk*abs(geo%h)
    if (p<=1_8) then
      if (p==0_8) then
        A=cos(z); A0=cos(z0)
      else
        A=cos(z)+z*sin(z); A0=cos(z0)+z0*sin(z0)
      endif
      ncomp_L=size(wk%weight_lambda,3,kind=8)/2_8
      allocate(wk_L%weight_lambda(nr,n_E,ncomp_L))
      wk_L%weight_lambda=0.0d0
      weight_0_L=0.0d0; weight_u_L=0.0d0; weight_v_L=0.0d0
      call compute_singular_correction_off(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk_L,weight_0_L(1:ncomp_L), &
        weight_u_L(1:ncomp_L),weight_v_L(1:ncomp_L))
      do ic=1,ncomp_L
        do j=1,n_E
          do i=1,nr
            wk%weight_lambda(i,j,ic)=wk%weight_lambda(i,j,ic)+real(A(i,j),8)*wk_L%weight_lambda(i,j,ic)
            wk%weight_lambda(i,j,ncomp_L+ic)=wk%weight_lambda(i,j,ncomp_L+ic)+aimag(A(i,j))*wk_L%weight_lambda(i,j,ic)
          enddo
        enddo
        weight_0(ic)=weight_0(ic)+real(A0,8)*weight_0_L(ic)
        weight_0(ncomp_L+ic)=weight_0(ncomp_L+ic)+aimag(A0)*weight_0_L(ic)
      enddo
    else
      call compute_singular_correction_Dn_helm_off(zk,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v, &
        chord,z,z0)
    endif
  end subroutine compute_singular_correction_helm_off

  subroutine compute_singular_correction_Dn_helm_off(zk,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v, &
    chord,z,z0)
    ! DLPn (p=2). With dist = |chord| and z = zk*dist, the target Hessian of exp(1i*z)/dist is
    !   A*H + A_I*I/dist + 1i*(s*chord*chord'-q*I),  A_I = (-zk^2/3)*cos(z),
    ! H the Laplace DLPn kernel, A = cos(z)+z*sin(z)-(z^2/3)*cos(z), q = zk^3*j1(z)/z, s = zk^5*j2(z)/z^2 (functions of dist^2).
    ! 1) A*H: step 7 (Laplace) and compute_singular_correction_off integrate H*psi, psi = A*mu, giving weights on psi_r,
    !    psi0, psi1. As in compute_singular_correction_stok_off, with A_u = A_v = 0 and A1r = lambda*A_r here:
    !      psi0 = A0*mu0, psi1 = A0*mu1, psi_r = A*mu_r + A1r*mu1 + A_r*mu0,
    !      A_r = (A-A0)/lambda^2 = zk^2*Q_lambda*F[z^2,z0^2], A0 = A(t') (dist = |h|).
    !    Q_lambda = (dist^2-h^2)/lambda^2 = |C0+lambda*c1r|^2-2*h*n_r0.c1r, F[z^2,z0^2] the divided difference of A as a
    !    function F of z^2, int_0^1 F'(z0^2+tau*(z^2-z0^2)) dtau with F'(x) = (cos(sqrt(x))+sqrt(x)*sin(sqrt(x)))/6, by
    !    16-point Gauss in tau.
    ! 2) A_I*I/dist: the SLP Gauss sum f*w_R*J*(lambda/dist)*A_I and A_I times its p=0 correction,
    !    W_I, and 3) 1i*(s*chord*chord'-q*I) by Gauss, f*w_R*lambda*J*(...): weights W on mu, put on mur, mu0, mu1 through
    !    mu = mu0+lambda*mu1+lambda^2*mur. Components (a,k) at a+3*(k-1), as compute_kernelradial_Dn_off; 1:9 real parts,
    !    10:18 imaginary parts.
    complex(8), intent(in) :: zk, z(:,:), z0
    integer(8), intent(in) :: order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:), chord(:,:,:)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    type(distance_factor) :: gp_I
    type(kernel_weights) :: wk_I
    integer(8) :: i, j, k, ia, ik, q_tau, n_term
    real(8) :: factor, n_r0(3), c(3), Q_lambda, tau(16), w_tau(16), dist, dij, delta
    real(8) :: weight_0_I(1), weight_u_I(1), weight_v_I(1)
    real(8) :: weight_0_in(9), weight_u_in(9), weight_v_in(9)
    complex(8) :: A(nr), A0, A_r(nr), F_div, s_tau, conv_0, conv_1, A_I, W_I(nr,n_E)
    complex(8) :: q, s, term, w_smooth, W, sum_0(9), sum_u(9), sum_v(9)

    n_r0=geo%N0/geo%J0
    ! 1) Laplace on psi, then psi -> mu (endpoint terms of this edge, A_u = A_v = 0)
    weight_0_in=weight_0(1:9); weight_u_in=weight_u(1:9); weight_v_in=weight_v(1:9)
    call compute_singular_correction_off(2_8,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0(1:9),weight_u(1:9),weight_v(1:9))
    A0=cos(z0)+z0*sin(z0)-z0**2/3.0d0*cos(z0)
    do k=1,9
      delta=weight_0(k)-weight_0_in(k)
      weight_0(9+k)=weight_0(9+k)+aimag(A0)*delta
      weight_0(k)=weight_0(k)+real(A0-1.0d0,8)*delta
      delta=weight_u(k)-weight_u_in(k)
      weight_u(9+k)=weight_u(9+k)+aimag(A0)*delta
      weight_u(k)=weight_u(k)+real(A0-1.0d0,8)*delta
      delta=weight_v(k)-weight_v_in(k)
      weight_v(9+k)=weight_v(9+k)+aimag(A0)*delta
      weight_v(k)=weight_v(k)+real(A0-1.0d0,8)*delta
    enddo
    call compute_gauss_legendre(16_8,tau,w_tau)
    tau=(tau+1.0d0)/2.0d0; w_tau=w_tau/2.0d0
    do j=1,n_E
      do i=1,nr
        A(i)=cos(z(i,j))+z(i,j)*sin(z(i,j))-z(i,j)**2/3.0d0*cos(z(i,j))
        c=X%C0(:,j)+R%lambda(i)*X%c1r(i,:,j)
        Q_lambda=((c(1)*c(1)+c(2)*c(2))+c(3)*c(3)) &
          -2.0d0*geo%h*((n_r0(1)*X%c1r(i,1,j)+n_r0(2)*X%c1r(i,2,j))+n_r0(3)*X%c1r(i,3,j))
        F_div=(0.0d0,0.0d0)
        do q_tau=1,16
          s_tau=sqrt(z0**2+tau(q_tau)*(z(i,j)**2-z0**2))
          F_div=F_div+w_tau(q_tau)*(cos(s_tau)+s_tau*sin(s_tau))/6.0d0
        enddo
        A_r(i)=zk**2*Q_lambda*F_div
      enddo
      do k=1,9
        conv_0=(0.0d0,0.0d0)
        conv_1=(0.0d0,0.0d0)
        do i=1,nr
          conv_0=conv_0+wk%weight_r(i,j,k)*A_r(i)
          conv_1=conv_1+wk%weight_r(i,j,k)*(R%lambda(i)*A_r(i))
          wk%weight_r(i,j,9+k)=wk%weight_r(i,j,k)*aimag(A(i))
          wk%weight_r(i,j,k)=wk%weight_r(i,j,k)*real(A(i),8)
        enddo
        weight_0(k)=weight_0(k)+real(conv_0,8); weight_0(9+k)=weight_0(9+k)+aimag(conv_0)
        weight_u(k)=weight_u(k)+real(conv_1,8)*d(1,j); weight_u(9+k)=weight_u(9+k)+aimag(conv_1)*d(1,j)
        weight_v(k)=weight_v(k)+real(conv_1,8)*d(2,j); weight_v(9+k)=weight_v(9+k)+aimag(conv_1)*d(2,j)
      enddo
    enddo
    ! 2) the I/dist term: the correction of its SLP Gauss sum, compute_singular_correction_off (p=0, g0 = q0^(-1/2)) from
    !    zero weights (wk_I), times a = (-zk^2/3)*cos(z) at the nodes: W_I
    allocate(gp_I%g0(n_E),gp_I%d_curv(nr,n_E))
    gp_I%g0=X%q0**(-0.5d0); gp_I%d_curv=gp%d_curv
    allocate(wk_I%weight_lambda(nr,n_E,1))
    wk_I%weight_lambda=0.0d0
    weight_0_I=0.0d0; weight_u_I=0.0d0; weight_v_I=0.0d0
    call compute_singular_correction_off(0_8,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp_I,wk_I,weight_0_I,weight_u_I,weight_v_I)
    do j=1,n_E
      do i=1,nr
        A_I=(-zk**2/3.0d0)*cos(z(i,j))
        W_I(i,j)=A_I*wk_I%weight_lambda(i,j,1)
      enddo
    enddo
    ! 3) the I/dist Gauss sum and the smooth term; q = zk^3*j1(z)/z, s = zk^5*j2(z)/z^2
    !    = zk^5*sum_m (-z^2/2)^m/(m!*(2m+5)!!): the series for |z|<1
    sum_0=(0.0d0,0.0d0); sum_u=(0.0d0,0.0d0); sum_v=(0.0d0,0.0d0)
    do j=1,n_E
      factor=geo%x_E(iedge)*w_s(j)/(2.0d0*pi)
      do i=1,nr
        if (abs(z(i,j))<1.0d0) then
          term=(1.0d0,0.0d0)/3.0d0; q=term
          do n_term=1,20
            term=term*(-z(i,j)**2/2.0d0)/real(n_term*(2*n_term+3),8)
            q=q+term
          enddo
          term=(1.0d0,0.0d0)/15.0d0; s=term
          do n_term=1,20
            term=term*(-z(i,j)**2/2.0d0)/real(n_term*(2*n_term+5),8)
            s=s+term
          enddo
        else
          q=(sin(z(i,j))-z(i,j)*cos(z(i,j)))/z(i,j)**3
          s=((3.0d0-z(i,j)**2)*sin(z(i,j))-3.0d0*z(i,j)*cos(z(i,j)))/z(i,j)**5
        endif
        q=zk**3*q; s=zk**5*s
        dist=sqrt((chord(i,1,j)*chord(i,1,j)+chord(i,2,j)*chord(i,2,j))+chord(i,3,j)*chord(i,3,j))
        A_I=(-zk**2/3.0d0)*cos(z(i,j))
        W_I(i,j)=W_I(i,j)+factor*R%w(i)*(Jx%J(i,j)*(R%lambda(i)/dist)*A_I)
        w_smooth=factor*R%w(i)*((0.0d0,1.0d0)*R%lambda(i)*Jx%J(i,j))
        do ik=1,3
          do ia=1,3
            k=ia+3*(ik-1)
            dij=merge(1.0d0,0.0d0,ia==ik)
            W=w_smooth*(s*(chord(i,ia,j)*chord(i,ik,j))-dij*q)+dij*W_I(i,j)
            wk%weight_r(i,j,k)=wk%weight_r(i,j,k)+R%lambda(i)**2*real(W,8)
            wk%weight_r(i,j,9+k)=wk%weight_r(i,j,9+k)+R%lambda(i)**2*aimag(W)
            sum_0(k)=sum_0(k)+W
            sum_u(k)=sum_u(k)+(R%lambda(i)*d(1,j))*W
            sum_v(k)=sum_v(k)+(R%lambda(i)*d(2,j))*W
          enddo
        enddo
      enddo
    enddo
    weight_0(1:9)=weight_0(1:9)+real(sum_0,8); weight_0(10:18)=weight_0(10:18)+aimag(sum_0)
    weight_u(1:9)=weight_u(1:9)+real(sum_u,8); weight_u(10:18)=weight_u(10:18)+aimag(sum_u)
    weight_v(1:9)=weight_v(1:9)+real(sum_v,8); weight_v(10:18)=weight_v(10:18)+aimag(sum_v)
  end subroutine compute_singular_correction_Dn_helm_off

  subroutine compute_singular_correction_stok_off(kind,p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v, &
    tnx)
    ! Step 8 of the Stokes off-surface drivers (h/=0): the Stokes parts that are a Laplace kernel times a smooth factor A(y).
    ! compute_singular_correction_off gives this edge's weights on psi = A*mu along the rays (radial nodes, t'); A is then
    ! applied there, before steps 9-10, so the weights integrate kernel*A*I16[mu] (multiplying finished GL16 weights by A at
    ! the nodes would integrate kernel*I16[A*mu]). chord = y-x = lambda*(C0+lambda*c1r)-h*n_r0, n_r0 = N0/J0.
    !   kind 'S' (SLPn pipeline, order 0, p=1): psi_b = A_b*mu_b, A_b = chord_b; weights (a,b) at a+3*(b-1), a the Sn
    !     component: weight_lambda(:,:,a+3*(b-1)) = weight_lambda_Sn(:,:,a)*A_b, weight_0(a+3*(b-1)) += weight_0_Sn(a)*A_b(t'),
    !     A_b(t') = -h*n_r0(b).
    !   kind 'D', 'Sn' (DLPn pipeline, order 2, p=2): A = (x-y).n_s, n_s = -N/J (the normal of WD, snx), or A = (x-y).tnx.
    !     With A = A0+lambda*(d_u*A_u+d_v*A_v)+lambda^2*A_r along a ray, A1r = (A-A0)/lambda, mu = mu0+lambda*mu1+lambda^2*mur,
    !     mu1 = d_u*mu_u+d_v*mu_v:
    !       (A*mu)(t') = A0*mu0,  grad(A*mu)(t') = A0*grad mu + mu0*grad A,  (A*mu)_r = A*mur + A1r*mu1 + A_r*mu0.
    !     'D':  A0 = -h, grad A(t') = 0, A_r = -H1.n_c+(h/2)*|B|^2 = (C0.N1r+c1r.N)/J+(h/2)*|B|^2, H1 = V1r-c1r,
    !           n_c = N/J, N = N0+lambda*N1r, N1r = N1+lambda*Nr, B = (n_c-n_r0)/lambda = (N1r-n_r0*J1r)/J
    !           (V.N = 0 gives C0.N1r = -V1r.N; stable: no difference of O(h) numbers);
    !     'Sn': A0 = h*n_r0.tnx, grad A(t') = -(r_u0.tnx, r_v0.tnx), A_r = -c1r.tnx.
    !     Endpoint weights are transformed for this edge only (weight - weight before step 8).
    !   kind 'Dn' (p=3): the shared correction handles the traction kernel at tnx directly, without conversion.
    character(len=*), intent(in) :: kind
    integer(8), intent(in) :: p, order, iedge, n_E
    type(target_geometry), intent(in) :: geo
    type(radial_tables), intent(in) :: R
    real(8), intent(in) :: d(:,:), w_s(:), tnx(3)
    type(displacement_expansion), intent(in) :: X
    type(jacobian_expansion), intent(in) :: Jx
    type(distance_factor), intent(in) :: gp
    type(kernel_weights), intent(inout) :: wk
    real(8), intent(inout) :: weight_0(:), weight_u(:), weight_v(:)
    integer(8) :: i, j, k, ia, ib
    real(8) :: n_r0(3), chord_b, weight_0_Sn(3), weight_u_Sn(3), weight_v_Sn(3), W_Sn(nr,n_E,3)
    real(8) :: weight_0_in(9), weight_u_in(9), weight_v_in(9), dw0, dwu, dwv, conv_0, conv_1, sum_0, sum_u, sum_v
    real(8) :: A0, A_u, A_v, A_r(nr,n_E), A1r(nr,n_E), A(nr,n_E), N1r(3), N(3), B(3), w

    n_r0=geo%N0/geo%J0
    if (kind=='S') then
      weight_0_Sn=0.0d0; weight_u_Sn=0.0d0; weight_v_Sn=0.0d0
      call compute_singular_correction_off(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0_Sn,weight_u_Sn,weight_v_Sn)
      W_Sn=wk%weight_lambda(:,1:n_E,1:3)
      do ib=1,3
        do ia=1,3
          k=ia+3*(ib-1)
          do j=1,n_E
            do i=1,nr
              chord_b=R%lambda(i)*(X%C0(ib,j)+R%lambda(i)*X%c1r(i,ib,j))-geo%h*n_r0(ib)
              wk%weight_lambda(i,j,k)=W_Sn(i,j,ia)*chord_b
            enddo
          enddo
          weight_0(k)=weight_0(k)+weight_0_Sn(ia)*(-geo%h*n_r0(ib))
        enddo
      enddo
      return
    endif
    if (kind=='D' .or. kind=='Sn') then
      weight_0_in=weight_0(1:9); weight_u_in=weight_u(1:9); weight_v_in=weight_v(1:9)
      call compute_singular_correction_off(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v)
      if (kind=='D') then
        A0=-geo%h; A_u=0.0d0; A_v=0.0d0
        do j=1,n_E
          do i=1,nr
            N1r=Jx%N1(1,:,j)+R%lambda(i)*Jx%Nr(i,:,j)
            N=Jx%N0(1,:)+R%lambda(i)*N1r
            B=(N1r-n_r0*Jx%J1r(i,j))/Jx%J(i,j)
            A_r(i,j)=(((X%C0(1,j)*N1r(1)+X%C0(2,j)*N1r(2))+X%C0(3,j)*N1r(3)) &
              +((X%c1r(i,1,j)*N(1)+X%c1r(i,2,j)*N(2))+X%c1r(i,3,j)*N(3)))/Jx%J(i,j) &
              +0.5d0*geo%h*((B(1)*B(1)+B(2)*B(2))+B(3)*B(3))
            A1r(i,j)=R%lambda(i)*A_r(i,j)
          enddo
        enddo
      else
        A0=geo%h*((n_r0(1)*tnx(1)+n_r0(2)*tnx(2))+n_r0(3)*tnx(3))
        A_u=-((geo%r_u0(1)*tnx(1)+geo%r_u0(2)*tnx(2))+geo%r_u0(3)*tnx(3))
        A_v=-((geo%r_v0(1)*tnx(1)+geo%r_v0(2)*tnx(2))+geo%r_v0(3)*tnx(3))
        do j=1,n_E
          do i=1,nr
            A_r(i,j)=-((X%c1r(i,1,j)*tnx(1)+X%c1r(i,2,j)*tnx(2))+X%c1r(i,3,j)*tnx(3))
            A1r(i,j)=(d(1,j)*A_u+d(2,j)*A_v)+R%lambda(i)*A_r(i,j)
          enddo
        enddo
      endif
      do j=1,n_E
        do i=1,nr
          A(i,j)=A0+R%lambda(i)*A1r(i,j)
        enddo
      enddo
      do k=1,9
        ! endpoint terms of this edge
        dw0=weight_0(k)-weight_0_in(k); dwu=weight_u(k)-weight_u_in(k); dwv=weight_v(k)-weight_v_in(k)
        weight_0(k)=((weight_0_in(k)+A0*dw0)+A_u*dwu)+A_v*dwv
        weight_u(k)=weight_u_in(k)+A0*dwu
        weight_v(k)=weight_v_in(k)+A0*dwv
        ! radial weights on mur: (A*mu)_r = A*mur + A1r*mu1 + A_r*mu0
        sum_0=0.0d0; sum_u=0.0d0; sum_v=0.0d0
        do j=1,n_E
          conv_0=0.0d0; conv_1=0.0d0
          do i=1,nr
            w=wk%weight_r(i,j,k)
            conv_0=conv_0+w*A_r(i,j)
            conv_1=conv_1+w*A1r(i,j)
            wk%weight_r(i,j,k)=w*A(i,j)
          enddo
          sum_0=sum_0+conv_0
          sum_u=sum_u+d(1,j)*conv_1
          sum_v=sum_v+d(2,j)*conv_1
        enddo
        weight_0(k)=weight_0(k)+sum_0
        weight_u(k)=weight_u(k)+sum_u
        weight_v(k)=weight_v(k)+sum_v
      enddo
      return
    endif
    ! DLPn (p=3): direct traction correction
    call compute_singular_correction_off(p,order,iedge,geo,R,n_E,d,w_s,X,Jx,gp,wk,weight_0,weight_u,weight_v,tnx)
  end subroutine compute_singular_correction_stok_off

  subroutine compute_muhat_weights(order,iedge,R,T_radial,n_E,d,wk,W_leg)
    ! Step 9: the step-7 weights on the density samples of one edge, turned into weights on the Legendre
    ! coefficients muhat_ab of mu = sum_ab muhat_ab*P_a(u)*P_b(v); W_leg(a,b,k), a u-mode, b v-mode.
    ! Fixed axis f, moving axis m: sum_ij w(i,j)*mu(t'+lambda_i*d_j) = sum_ab muhat_ab*[Pfixed'*C'(w*Pedge)]_ab,
    !   (C'M)(i,b) = sum_k T(i,k,b)*M(i,k) (T, T_dd: compute_legendre_translation); the remainders use T_dd the same way:
    !   mu1r = Pf1r x Pm0 + Pfixed x (dm*Pm1 + lambda*dm^2*T_dd part)   (order 1, SLPn),
    !   mur  = Pfr  x Pm0 + Pf1r  x dm*Pm1 + Pfixed x dm^2*T_dd part      (order 2, DLPn).
    ! wk%weight_lambda (SLP, DLP, SLPn) is on mu, wk%weight_r on mu1r or mur. Endpoint terms: step 10.
    ! With v fixed (edges 1, 3) the edge term is added transposed.
    integer(8), intent(in) :: order, iedge, n_E
    type(radial_tables), intent(in) :: R
    type(radial_expansion_tables), intent(in) :: T_radial
    real(8), intent(in) :: d(:,:)
    type(kernel_weights), intent(in) :: wk
    real(8), intent(inout) :: W_leg(:,:,:)
    integer(8) :: fixed_axis, moving_axis, k, i, j, ell, ell2, ib
    real(8) :: M(nr,16), M_c(nr,16), W_edge(16,16), W_r(16,16), sum_w(nr), w_dm(nr)
    real(8) :: Pw_m0(16), Pw_m1(16)
    real(8), allocatable :: w(:,:)
    external :: dgemm, dgemv

    fixed_axis=1_8+mod(iedge,2_8)
    moving_axis=3_8-fixed_axis
    if (order>=1_8) allocate(w(nr,n_E))
    do k=1,size(W_leg,3,kind=8)
      W_edge=0.0d0
      if (order<=1_8) then
        call dgemm('N','N',nr,16_8,n_E,1.0d0,wk%weight_lambda(1,1,k),nr,T_radial%Pedge,size(T_radial%Pedge,1,kind=8),0.0d0,M,nr)
        do ell2=1,16
          !GCC$ ivdep
          do i=1,nr
            M_c(i,ell2)=T_radial%T(i,1,ell2)*M(i,1)
          enddo
          do ell=2,ell2
            !GCC$ ivdep
            do i=1,nr
              M_c(i,ell2)=M_c(i,ell2)+T_radial%T(i,ell,ell2)*M(i,ell)
            enddo
          enddo
        enddo
        call dgemm('T','N',16_8,16_8,nr,1.0d0,T_radial%Pfixed,nr,M_c,nr,0.0d0,W_edge,16_8)
      endif
      if (order>=1_8) then
        w(:,:)=wk%weight_r(:,1:n_E,k)
        sum_w=0.0d0
        w_dm=0.0d0
        do j=1,n_E
          sum_w=sum_w+w(:,j)
          w_dm=w_dm+w(:,j)*d(moving_axis,j)
        enddo
        if (order==1_8) then
          call dgemv('T',nr,16_8,1.0d0,T_radial%Pf1r,nr,sum_w,1_8,0.0d0,Pw_m0,1_8)
          call dgemv('T',nr,16_8,1.0d0,T_radial%Pfixed,nr,w_dm,1_8,0.0d0,Pw_m1,1_8)
          do j=1,n_E
            w(:,j)=R%lambda*w(:,j)*d(moving_axis,j)**2
          enddo
        else
          call dgemv('T',nr,16_8,1.0d0,T_radial%Pfr,nr,sum_w,1_8,0.0d0,Pw_m0,1_8)
          call dgemv('T',nr,16_8,1.0d0,T_radial%Pf1r,nr,w_dm,1_8,0.0d0,Pw_m1,1_8)
          do j=1,n_E
            w(:,j)=w(:,j)*d(moving_axis,j)**2
          enddo
        endif
        call dgemm('N','N',nr,16_8,n_E,1.0d0,w,nr,T_radial%Pedge,size(T_radial%Pedge,1,kind=8),0.0d0,M,nr)
        ! T_dd(:,ell,ell2) = 0 for ell > ell2-2, so only those terms are summed.
        M_c(:,1:2)=0.0d0
        do ell2=3,16
          !GCC$ ivdep
          do i=1,nr
            M_c(i,ell2)=T_radial%T_dd(i,1,ell2)*M(i,1)
          enddo
          do ell=2,ell2-2
            !GCC$ ivdep
            do i=1,nr
              M_c(i,ell2)=M_c(i,ell2)+T_radial%T_dd(i,ell,ell2)*M(i,ell)
            enddo
          enddo
        enddo
        call dgemm('T','N',16_8,16_8,nr,1.0d0,T_radial%Pfixed,nr,M_c,nr,0.0d0,W_r,16_8)
        do ib=1,16
          W_r(:,ib)=W_r(:,ib)+Pw_m0*T_radial%Pm0(ib)+Pw_m1*T_radial%Pm1(ib)
        enddo
        W_edge=W_edge+W_r
      endif
      if (fixed_axis==1_8) then
        W_leg(:,:,k)=W_leg(:,:,k)+W_edge
      else
        W_leg(:,:,k)=W_leg(:,:,k)+transpose(W_edge)
      endif
    enddo
  end subroutine compute_muhat_weights

  subroutine compute_final_weights(density_map,geo,weight_0,weight_u,weight_v,W_leg,W)
    ! Step 10: add the endpoint terms once and map Legendre-coefficient weights to the GL16 nodes.
    ! W_leg(:,:,k) + weight_0*P(u') x P(v') + weight_u*P'(u') x P(v') + weight_v*P(u') x P'(v'),
    ! then density_map*W_leg*density_map'; W is row itarg of the driver's W(m,256,ncomp), here (256,ncomp).
    real(8), intent(in) :: density_map(16,16), weight_0(:), weight_u(:), weight_v(:), W_leg(:,:,:)
    type(target_geometry), intent(in) :: geo
    real(8), intent(inout) :: W(:,:)
    integer(8) :: k, ia, ib
    real(8) :: W_k(16,16), W_tmp(16,16), W_nodes(16,16)
    external :: dgemm

    do k=1,size(W_leg,3,kind=8)
      do ib=1,16
        do ia=1,16
          W_k(ia,ib)=((W_leg(ia,ib,k)+weight_0(k)*(geo%P0(ia,1)*geo%P0(ib,2)))+weight_u(k)*(geo%dP0(ia,1)*geo%P0(ib,2))) &
            +weight_v(k)*(geo%P0(ia,1)*geo%dP0(ib,2))
        enddo
      enddo
      call dgemm('N','N',16_8,16_8,16_8,1.0d0,density_map,16_8,W_k,16_8,0.0d0,W_tmp,16_8)
      call dgemm('N','T',16_8,16_8,16_8,1.0d0,W_tmp,16_8,density_map,16_8,0.0d0,W_nodes,16_8)
      W(:,k)=reshape(W_nodes,[256])
    enddo
  end subroutine compute_final_weights

end module surfquadr
