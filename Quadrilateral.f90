!======================================================
!======================================================
!=====           FINITE ELEMENT METHOD            =====
!=====    QUADRILATERAL ELEMENTS WITH 4-9 NODES   =====
!=====             MOHAMMAD SAJEDNIA              =====
!======================================================
!======================================================
program main
  implicit none
  integer :: nnd, nel, nbc, nload, ntrac, nvol, i, ios
  integer, allocatable :: elements(:,:)
  integer, parameter :: maxx=10000
  real(8), dimension(:,:), allocatable :: nodes, nodal_loads, tractions, volume_loads, bc
  real :: E, nu, plane
  real(8) :: thickness
  real(8), allocatable :: D(:,:), K(:,:), F(:), U(:), reaction_forces(:)
  character(len=20) :: filename

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!             !!!!!!!!!!!!!
  !!!!!!!!!!!!! READ INPUTS !!!!!!!!!!!!!
  !!!!!!!!!!!!!             !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  filename = 'nodes.txt'
  open(unit=10, file=filename, status='old', action='read')
  allocate(nodes(maxx,3))
  nodes = 0.0
  nnd=0
  do i=1,maxx
    read(10,*,iostat=ios) nodes(i,:)
    if (ios/=0) exit
    nnd=nnd + 1
  end do
  close(10)

  do i=1, nnd
    nodes(i,1)=nodes(i,2)
    nodes(i,2)=nodes(i,3)
  end do

  filename = 'elements.txt'
  open(unit=11, file=filename, status='old', action='read')
  allocate(elements(maxx,10))
  elements = 0
  nel=0
  do i=1,maxx
    read(11,*,iostat=ios) elements(i,:)
    if (ios/=0) exit
    nel=nel+1
  end do
  close(11)

  do i=1, nel
    elements(i,1)=elements(i,2)
    elements(i,2)=elements(i,3)
    elements(i,3)=elements(i,4)
    elements(i,4)=elements(i,5)
    elements(i,5)=elements(i,6)
    elements(i,6)=elements(i,7)
    elements(i,7)=elements(i,8)
    elements(i,8)=elements(i,9)
    elements(i,9)=elements(i,10)
  end do

  filename = 'boundaries.txt'
  open(unit=12, file=filename, status='old', action='read')
  allocate(bc(maxx,3))
  bc = 0.0
  nbc=0
  do i=1,maxx
    read(12,*,iostat=ios) bc(i,:)
    if (ios/=0) exit
    nbc=nbc+1
  end do
  close(12)

  filename = 'nloads.txt'
  open(unit=13, file=filename, status='old', action='read')
  allocate(nodal_loads(maxx,3))
  nodal_loads = 0.0
  nload=0
  do i=1,maxx
    read(13,*,iostat=ios) nodal_loads(i,:)
    if (ios/=0) exit
    nload=nload+1
  end do
  close(13)

  filename = 'tractions.txt'
  open(unit=14, file=filename, status='old', action='read')
  allocate(tractions(maxx,8))
  tractions = 0.0
  ntrac=0
  do i=1,maxx
    read(14,*,iostat=ios) tractions(i,:)
    if (ios/=0) exit
    ntrac=ntrac+1
  end do
  close(14)

  filename = 'vloads.txt'
  open(unit=15, file=filename, status='old', action='read')
  allocate(volume_loads(maxx,3))
  volume_loads = 0.0
  nvol=0
  do i=1,maxx
    read(15,*,iostat=ios) volume_loads(i,:)
    if (ios/=0) exit
    nvol=nvol+1
  end do
  close(15)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!! MATERIAL'S PROPERTIES !!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  E = 210e6
  nu = 0.3
  thickness = 0.025
  plane = 0 ! 0 if plane stress , 1 if plane strain

  if (plane==0) then
    allocate(D(3,3))
    D = (E/(1.0-nu**2.0)) * reshape([1.0, nu, 0.0, nu, 1.0, 0.0, 0.0, 0.0, (1.0-nu)/2.0], [3,3])
  elseif (plane==1) then
    allocate(D(3,3))
    D = (E/(1.0+nu)/(1.0-2.0*nu)) * reshape([1.0-nu, -nu, 0.0, -nu, 1.0-nu, 0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0], [3,3])
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                     !!!!!!!!
  !!!!!!!! CALLING SUBROUTINES !!!!!!!!
  !!!!!!!!                     !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(F(2*nnd))
  F = 0.0
  call assemble_forces(thickness, nodes, nodal_loads, nload, tractions, ntrac, volume_loads, nvol, elements, F)

  allocate(K(2*nnd,2*nnd))
  K = 0.0
  call assemble_stiffness_matrix(nodes, nnd, elements, nel, thickness, D, K)

  allocate(U(2*nnd))
  U = 0.0
  call bcs_displacements(K, F, U, nbc, bc, nnd)

  allocate(reaction_forces(2*nnd))
  reaction_forces = 0.0
  call r_forces(K, U, bc, nbc, reaction_forces)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!    WRITING RESULTS    !!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  filename = 'Outputs.txt'
  open(unit=16, file=filename, status='replace', action='write')
  write(16,*) 'Nodal Displacements :'
  do i = 1, nnd
    write(16,*) 'Node', i, ': X =', U(2*i-1), ', Y =', U(2*i)
  end do
  write(16,*) 'Reaction Forces :'
  do i = 1, nnd
    write(16,*) 'Node', i, ': X =', reaction_forces(2*i-1), ', Y =', reaction_forces(2*i)
  end do
  close(16)

  print*, ''
  print*, 'Check Outputs.txt !'

  deallocate(nodes, elements, bc, D, K, F, U, reaction_forces, nodal_loads, tractions,volume_loads)

  stop

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!! GAUSS ELIMINATION METHOD !!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gauss(k,f,u)
  implicit none
  integer :: n
  real(8), dimension(:), intent(inout) :: u, f
  real(8), dimension(:,:), intent(inout) :: k
  real :: temp
  integer :: i, j, p, q
  real :: factor
  n = size(u)

  do i = 1, n-1
    p = i
    do j = i+1, n
      if (abs(k(j,i)) > abs(k(p,i))) then
        p = j
      end if
    end do

    if (p /= i) then
      do q = 1, n
        temp = k(i,q)
        k(i,q) = k(p,q)
        k(p,q) = temp
      end do
      temp = f(i)
      f(i) = f(p)
      f(p) = temp
    end if

    do j = i+1, n
      factor = k(j,i) / k(i,i)
      k(j,i) = 0.0
      do q = i+1, n
        k(j,q) = k(j,q) - factor*k(i,q)
      end do
      f(j) = f(j) - factor*f(i)
    end do
  end do

  u(n) = f(n) / k(n,n)
  do i = n-1, 1, -1
    u(i) = f(i)
    do j = i+1, n
      u(i) = u(i) - k(i,j)*u(j)
    end do
    u(i) = u(i) / k(i,i)
  end do

end subroutine gauss

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                           !!!!!!!!!!!!!
  !!!!!!!!!!!!!  MATRICES MULTIPLICATION  !!!!!!!!!!!!!
  !!!!!!!!!!!!!                           !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sub_matmul(A, B, C, m, n, p)
  implicit none
  integer :: m, n, p
  real(8) :: A(m,n), B(n,p), C(m,p)
  integer :: i, j, k
  real :: sum0
  do i = 1, m
    do j = 1, p
      sum0 = 0.0
      do k = 1, n
        sum0 = sum0 + A(i,k) * B(k,j)
      C(i,j) = sum0
      end do
    end do
  end do
end subroutine sub_matmul

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!    MATRICES TRANSPOSE    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sub_transpose(a, at)
  implicit none
  real(8), dimension(:,:), intent(in) :: a
  real(8), dimension(size(a, 2), size(a, 1)), intent(out) :: at
  integer :: i, j
  do i = 1, size(a, 1)
    do j = 1, size(a, 2)
      at(j, i) = a(i, j)
    end do
  end do
end subroutine sub_transpose

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                     !!!!!!!!!!!!!
  !!!!!!!!!!!!!    FORCES VECTOR    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                     !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assemble_forces(thickness, nodes, nodal_loads, nload, tractions, ntrac, volume_loads, nvol, elements, F)
    implicit none

    real(8), intent(in) :: thickness
    real(8), dimension(:, :), intent(in) :: nodes
    real(8), dimension(:, :), intent(in) :: nodal_loads, tractions, volume_loads
    integer, dimension(:, :), intent(in) :: elements
    integer, intent(in) :: nload, ntrac, nvol
    real(8), dimension(2*nnd, 1), intent(out) :: F

    integer :: node_counter(nel,1)
    integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, nnod
    integer :: i, j, ki, Fnode, num_nodes, start_node, middle_node, end_node, vol_element
    integer :: mm0, nn0, pp0, ndel
    real(8) :: x1, x2, x3, x4, x5, x6, x7, x8, x9, y1, y2, y3, y4, y5, y6, y7, y8, y9, r, s
    real(8) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, detj
    real(8) :: f1r, f2r, f3r, f4r, f5r, f6r, f7r, f8r, f9r
    real(8) :: f1s, f2s, f3s, f4s, f5s, f6s, f7s, f8s, f9s
    real(8) :: fx, fy, Dx, Dy, L, A, bx, by, pi
    real(8) :: Fmag, Fangle, start_mag, middle_mag, end_mag, traction_angle, vol_mag, vol_angle
    real(8), dimension(2,2) :: dj
    real(8), allocatable, dimension(:,:) :: t, t1, n, nnt, nt, txi, tyi, tx, ty, drs, xy, hxy, hm

    pi = 3.14159265359

do i = 1, nload
    Fnode = nodal_loads(i, 1)
    Fmag = nodal_loads(i, 2)
    Fangle = nodal_loads(i, 3) / 180 * pi

    fx = Fmag * cos(Fangle)
    fy = Fmag * sin(Fangle)

    F(2 * Fnode - 1, 1) = F(2 * Fnode - 1, 1) + fx
    F(2 * Fnode, 1) = F(2 * Fnode, 1) + fy
end do

    node_counter=0
    do i=1, nel
        do j=1,9
            if (elements(i,j)/=0) then
            node_counter(i,1)=node_counter(i,1)+1
            end if
        end do
    end do

do j = 1, ntrac
    num_nodes = tractions(j, 1)
    start_node = tractions(j, 2)
    end_node = tractions(j, 4)
    start_mag = tractions(j, 5)
    end_mag = tractions(j, 7)
    allocate(t(num_nodes,1))
    allocate(t1(num_nodes,1))
    allocate(n(num_nodes,1))
    allocate(nt(1,num_nodes))
    allocate(txi(num_nodes,1))
    allocate(tyi(num_nodes,1))
    allocate(tx(num_nodes,1))
    allocate(ty(num_nodes,1))
    allocate(nnt(num_nodes,num_nodes))

    t(1,1) = start_mag
    t(2,1) = end_mag
    if (num_nodes==3) then
        middle_node = tractions(j, 3)
        middle_mag = tractions(j, 6)
        t(2,1) = middle_mag
        t(3,1) = end_mag
    end if
    traction_angle = tractions(j, 8) / 180 * pi

    x1 = nodes(start_node, 1)
    x2 = nodes(end_node, 1)
    Dx = x2 - x1
    y1 = nodes(start_node, 2)
    y2 = nodes(end_node, 2)
    Dy = y2 - y1
    L = sqrt(Dx**2 + Dy**2)
    detj = L/2

    do i=1,2
        txi = 0.0
        tyi = 0.0
        n = 0.0
        nt = 0.0
        t1 = 0.0
        if (i==1) then
            r = -0.5773502692
            else
            r = 0.5773502692
        end if

    !!!!! SHAPE FUNCTIONS !!!!!

        if (num_nodes==2) then
            f1 = 0.5*(1-r)
            f2 = 0.5*(1+r)
            n(1,1) = f1
            n(2,1) = f2
        end if
        if (num_nodes==3) then
            f1 = -0.5*r*(1-r)
            f3 = 0.5*r*(1+r)
            f2 = (1-r)*(1+r)
            n(1,1) = f1
            n(2,1) = f2
            n(3,1) = f3
        end if

        call sub_transpose(n, nt)
        call sub_matmul(n, nt, nnt, num_nodes, 1, num_nodes)
        call sub_matmul(nnt, t, t1, num_nodes, num_nodes, 1)
        Txi = t1 * thickness * cos(traction_angle) * detj
        Tyi = t1 * thickness * sin(traction_angle) * detj
        Tx = Tx + Txi
        Ty = Ty + Tyi
    end do

    if (num_nodes==2) then
    F(2 * start_node - 1, 1) = F(2 * start_node - 1, 1) + Tx(1,1)
    F(2 * start_node, 1) = F(2 * start_node, 1) + Ty(1,1)
    F(2 * end_node - 1, 1) = F(2 * end_node - 1, 1) + Tx(2,1)
    F(2 * end_node, 1) = F(2 * end_node, 1) + Ty(2,1)
        else
        F(2 * start_node - 1, 1) = F(2 * start_node - 1, 1) + Tx(1,1)
        F(2 * start_node, 1) = F(2 * start_node, 1) + Ty(1,1)
        F(2 * end_node - 1, 1) = F(2 * end_node - 1, 1) + Tx(3,1)
        F(2 * end_node, 1) = F(2 * end_node, 1) + Ty(3,1)
        F(2 * middle_node - 1, 1) = F(2 * middle_node - 1, 1) + Tx(2,1)
        F(2 * middle_node, 1) = F(2 * middle_node, 1) + Ty(2,1)
    end if
    deallocate(n, nt, nnt, t, t1, txi, tyi, tx, ty)
end do

do ki = 1, nvol
    vol_element = volume_loads(ki, 1)
    vol_mag = volume_loads(ki, 2)
    vol_angle = volume_loads(ki, 3) / 180 * pi
    bx = vol_mag * cos(vol_angle)
    by = vol_mag * sin(vol_angle)
    nnod = node_counter(vol_element,1)
    ndel = nnod * 2

    allocate(drs(nnod,2))
    allocate(hm(nnod,2))
    allocate(xy(2,nnod))
    allocate(hxy(nnod,2))

    hxy = 0.0

    n1 = elements(vol_element, 1)
    x1 = nodes(n1, 1)
    y1 = nodes(n1, 2)
    n2 = elements(vol_element, 2)
    x2 = nodes(n2, 1)
    y2 = nodes(n2, 2)
    n3 = elements(vol_element, 3)
    x3 = nodes(n3, 1)
    y3 = nodes(n3, 2)
    n4 = elements(vol_element, 4)
    x4 = nodes(n4, 1)
    y4 = nodes(n4, 2)
    if (5 <= nnod) then
        n5 = elements(vol_element, 5)
        x5 = nodes(n5, 1)
        y5 = nodes(n5, 2)
        if (6 <= nnod) then
            n6 = elements(vol_element, 6)
            x6 = nodes(n6, 1)
            y6 = nodes(n6, 2)
            if (7 <= nnod) then
                n7 = elements(vol_element, 7)
                x7 = nodes(n7, 1)
                y7 = nodes(n7, 2)
                if (8 <= nnod) then
                    n8 = elements(vol_element, 8)
                    x8 = nodes(n8, 1)
                    y8 = nodes(n8, 2)
                    if (9 <= nnod) then
                        n9 = elements(vol_element, 9)
                        x9 = nodes(n9, 1)
                        y9 = nodes(n9, 2)
                    end if
                end if
            end if
        end if
    end if

    do j=1,4
        if (j==1) then
            r=-0.5773502692
            s=-0.5773502692
        end if
        if (j==2) then
            r=0.5773502692
            s=-0.5773502692
        end if
        if (j==3) then
            r=-0.5773502692
            s=0.5773502692
        end if
        if (j==4) then
            r=0.5773502692
            s=0.5773502692
        end if

    !!!!! SHAPE FUNCTIONS AND DERIVATIONS !!!!!

    f1=0.25*(1+r)*(1-s)
    f1r=0.25*(1-s)
    f1s=-0.25*(1+r)
    f2=0.25*(1+r)*(1+s)
    f2r=0.25*(1+s)
    f2s=0.25*(1+r)
    f3=0.25*(1-r)*(1+s)
    f3r=-0.25*(1+s)
    f3s=0.25*(1-r)
    f4=0.25*(1-r)*(1-s)
    f4r=-0.25*(1-s)
    f4s=-0.25*(1-r)
    f5=0.5*(1+r)*(1-s**2)
    f5r=0.5*(1-s**2)
    f5s=-(1+r)*s
    f6=0.5*(1+s)*(1-r**2)
    f6r=-(1+s)*r
    f6s=0.5*(1-r**2)
    f7=0.5*(1-r)*(1-s**2)
    f7r=-0.5*(1-s**2)
    f7s=-(1-r)*s
    f8=0.5*(1-s)*(1-r**2)
    f8r=-(1-s)*r
    f8s=-0.5*(1-r**2)
    f9=(1-r**2)*(1-s**2)
    f9r=-2*r*(1-s**2)
    f9s=-2*s*(1-r**2)

    hm(1,1)=f1
    hm(1,2)=f1
    hm(2,1)=f2
    hm(2,2)=f2
    hm(3,1)=f3
    hm(3,2)=f3
    hm(4,1)=f4
    hm(4,2)=f4

    if (5<=nnod) then
        f1=f1-0.5*f5
        f1r=f1r-0.5*f5r
        f1s=f1s-0.5*f5s
        f2=f2-0.5*f5
        f2r=f2r-0.5*f5r
        f2s=f2s-0.5*f5s
        hm(5,1)=f5
        hm(5,2)=f5
    end if

    if (6<=nnod) then
        f2=f2-0.5*f6
        f2r=f2r-0.5*f6r
        f2s=f2s-0.5*f6s
        f3=f3-0.5*f6
        f3r=f3r-0.5*f6r
        f3s=f3s-0.5*f6s
        hm(6,1)=f6
        hm(6,2)=f6
    end if

    if (7<=nnod) then
        f3=f3-0.5*f7
        f3r=f3r-0.5*f7r
        f3s=f3s-0.5*f7s
        f4=f4-0.5*f7
        f4r=f4r-0.5*f7r
        f4s=f4s-0.5*f7s
        hm(7,1)=f7
        hm(7,2)=f7
    end if

    if (8<=nnod) then
        f1=f1-0.5*f8
        f1r=f1r-0.5*f8r
        f1s=f1s-0.5*f8s
        f4=f4-0.5*f8
        f4r=f4r-0.5*f8r
        f4s=f4s-0.5*f8s
        hm(8,1)=f8
        hm(8,2)=f8
    end if

    if (9<=nnod) then
        f1=f1+0.25*f9
        f1r=f1r+0.25*f9r
        f1s=f1s+0.25*f9s
        f2=f2+0.25*f9
        f2r=f2r+0.25*f9r
        f2s=f2s+0.25*f9s
        f3=f3+0.25*f9
        f3r=f3r+0.25*f9r
        f3s=f3s+0.25*f9s
        f4=f4+0.25*f9
        f4r=f4r+0.25*f9r
        f4s=f4s+0.25*f9s
        f5=f5-0.5*f9
        f5r=f5r-0.5*f9r
        f5s=f5s-0.5*f9s
        f6=f6-0.5*f9
        f6r=f6r-0.5*f9r
        f6s=f6s-0.5*f9s
        f7=f7-0.5*f9
        f7r=f7r-0.5*f9r
        f7s=f7s-0.5*f9s
        f8=f8-0.5*f9
        f8r=f8r-0.5*f9r
        f8s=f8s-0.5*f9s
        hm(9,1)=f9
        hm(9,2)=f9
    end if

    !!!!! JACOBIAN !!!!!

    xy(1,1)=x1
    xy(2,1)=y1
    xy(1,2)=x2
    xy(2,2)=y2
    xy(1,3)=x3
    xy(2,3)=y3
    xy(1,4)=x4
    xy(2,4)=y4

    if (5<=nnod) then
        xy(1,5)=x5
        xy(2,5)=y5
    end if

    if (6<=nnod) then
        xy(1,6)=x6
        xy(2,6)=y6
    end if

    if (7<=nnod) then
        xy(1,7)=x7
        xy(2,7)=y7
    end if

    if (8<=nnod) then
        xy(1,8)=x8
        xy(2,8)=y8
    end if

    if (9<=nnod) then
        xy(1,9)=x9
        xy(2,9)=y9
    end if

    drs(1,1)=f1r
    drs(1,2)=f1s
    drs(2,1)=f2r
    drs(2,2)=f2s
    drs(3,1)=f3r
    drs(3,2)=f3s
    drs(4,1)=f4r
    drs(4,2)=f4s

    if (5<=nnod) then
        drs(5,1)=f5r
        drs(5,2)=f5s
    end if

    if (6<=nnod) then
        drs(6,1)=f6r
        drs(6,2)=f6s
    end if

    if (7<=nnod) then
        drs(7,1)=f7r
        drs(7,2)=f7s
    end if

    if (8<=nnod) then
        drs(8,1)=f8r
        drs(8,2)=f8s
    end if

    if (9<=nnod) then
        drs(9,1)=f9r
        drs(9,2)=f9s
    end if

    mm0 = size(xy,1)
    nn0 = size(xy,2)
    pp0 = size(drs,2)

    call sub_matmul(xy, drs, dj, mm0, nn0, pp0)

    detj=abs(dj(1,1)*dj(2,2)-dj(1,2)*dj(2,1))

    do i=1,nnod
        hm(i,1) = hm(i,1) * detj * thickness * bx
        hm(i,2) = hm(i,2) * detj * thickness * by
    end do

    hxy = hxy + hm

end do

    !!!!! ASSEMBLING !!!!!

    F(2 * n1 - 1, 1) = F(2 * n1 - 1, 1) + hxy(1,1)
    F(2 * n1, 1) = F(2 * n1, 1) + hxy(1,2)
    F(2 * n2 - 1, 1) = F(2 * n2 - 1, 1) + hxy(2,1)
    F(2 * n2, 1) = F(2 * n2, 1) + hxy(2,2)
    F(2 * n3 - 1, 1) = F(2 * n3 - 1, 1) + hxy(3,1)
    F(2 * n3, 1) = F(2 * n3, 1) + hxy(3,2)
    F(2 * n4 - 1, 1) = F(2 * n4 - 1, 1) + hxy(4,1)
    F(2 * n4, 1) = F(2 * n4, 1) + hxy(4,2)
    if (5<=nnod) then
    F(2 * n5 - 1, 1) = F(2 * n5 - 1, 1) + hxy(5,1)
    F(2 * n5, 1) = F(2 * n5, 1) + hxy(5,2)
    end if
    if (6<=nnod) then
    F(2 * n6 - 1, 1) = F(2 * n6 - 1, 1) + hxy(6,1)
    F(2 * n6, 1) = F(2 * n6, 1) + hxy(6,2)
    end if
    if (7<=nnod) then
    F(2 * n7 - 1, 1) = F(2 * n7 - 1, 1) + hxy(7,1)
    F(2 * n7, 1) = F(2 * n7, 1) + hxy(7,2)
    end if
    if (8<=nnod) then
    F(2 * n8 - 1, 1) = F(2 * n8 - 1, 1) + hxy(8,1)
    F(2 * n8, 1) = F(2 * n8, 1) + hxy(8,2)
    end if
    if (9<=nnod) then
    F(2 * n9 - 1, 1) = F(2 * n9 - 1, 1) + hxy(9,1)
    F(2 * n9, 1) = F(2 * n9, 1) + hxy(9,2)
    end if
    deallocate(drs,hxy,xy,hm)
end do

end subroutine assemble_forces

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                               !!!!!!!!!!!!!
  !!!!!!!!!!!!!    GLOBAL STIFFNESS MATRIX    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                               !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assemble_stiffness_matrix(nodes, nnd, elements, nel, thickness, D, K)

  implicit none

  integer, intent(in) :: nnd, nel
  real(8), dimension(:, :), intent(in) :: nodes
  integer, dimension(:, :), intent(in) :: elements
  real(8), intent(in) :: thickness
  real(8), dimension(3, 3), intent(in) :: D
  real(8), dimension(2*nnd, 2*nnd), intent(out) :: K

  integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9
  integer :: i, j, mm1, nn1, pp1, mm2, nn2, pp2, ndel, nnod
  integer :: node_counter(nel,1)
  real(8) :: x1, x2, x3, x4, x5, x6, x7, x8, x9, y1, y2, y3, y4, y5, y6, y7, y8, y9
  real(8) :: A, m11, m21, m31, m12, m22, m32, m13, m23, m33
  real(8), allocatable, dimension(:,:) :: k1

  node_counter=0
  do i=1, nel
    do j=1,9
        if (elements(i,j)/=0) then
        node_counter(i,1)=node_counter(i,1)+1
        end if
    end do
  end do

  K = 0.0

  do i = 1, nel
    nnod = node_counter(i,1)
    ndel = nnod*2

    n1 = elements(i, 1)
    n2 = elements(i, 2)
    n3 = elements(i, 3)
    n4 = elements(i, 4)
    if (5 <= nnod) then
        n5 = elements(i, 5)
        if (6 <= nnod) then
            n6 = elements(i, 6)
            if (7 <= nnod) then
                n7 = elements(i, 7)
                if (8 <= nnod) then
                    n8 = elements(i, 8)
                    if (9 <= nnod) then
                        n9 = elements(i, 9)
                    end if
                end if
            end if
        end if
    end if

    allocate(k1(ndel,ndel))
    k1 = 0.0
    call local_stiffness(i, nodes, nnd, elements, nel, thickness, ndel, D, k1)

    K(2*n1-1,2*n1-1) = K(2*n1-1,2*n1-1) + k1(1,1);
    K(2*n1-1,2*n1) = K(2*n1-1,2*n1) + k1(1,2);
    K(2*n1-1,2*n2-1) = K(2*n1-1,2*n2-1) + k1(1,3);
    K(2*n1-1,2*n2) = K(2*n1-1,2*n2) + k1(1,4);
    K(2*n1-1,2*n3-1) = K(2*n1-1,2*n3-1) + k1(1,5);
    K(2*n1-1,2*n3) = K(2*n1-1,2*n3) + k1(1,6);
    K(2*n1-1,2*n4-1) = K(2*n1-1,2*n4-1) + k1(1,7);
    K(2*n1-1,2*n4) = K(2*n1-1,2*n4) + k1(1,8);

    K(2*n1,2*n1-1) = K(2*n1,2*n1-1) + k1(2,1);
    K(2*n1,2*n1) = K(2*n1,2*n1) + k1(2,2);
    K(2*n1,2*n2-1) = K(2*n1,2*n2-1) + k1(2,3);
    K(2*n1,2*n2) = K(2*n1,2*n2) + k1(2,4);
    K(2*n1,2*n3-1) = K(2*n1,2*n3-1) + k1(2,5);
    K(2*n1,2*n3) = K(2*n1,2*n3) + k1(2,6);
    K(2*n1,2*n4-1) = K(2*n1,2*n4-1) + k1(2,7);
    K(2*n1,2*n4) = K(2*n1,2*n4) + k1(2,8);

    K(2*n2-1,2*n1-1) = K(2*n2-1,2*n1-1) + k1(3,1);
    K(2*n2-1,2*n1) = K(2*n2-1,2*n1) + k1(3,2);
    K(2*n2-1,2*n2-1) = K(2*n2-1,2*n2-1) + k1(3,3);
    K(2*n2-1,2*n2) = K(2*n2-1,2*n2) + k1(3,4);
    K(2*n2-1,2*n3-1) = K(2*n2-1,2*n3-1) + k1(3,5);
    K(2*n2-1,2*n3) = K(2*n2-1,2*n3) + k1(3,6);
    K(2*n2-1,2*n4-1) = K(2*n2-1,2*n4-1) + k1(3,7);
    K(2*n2-1,2*n4) = K(2*n2-1,2*n4) + k1(3,8);

    K(2*n2,2*n1-1) = K(2*n2,2*n1-1) + k1(4,1);
    K(2*n2,2*n1) = K(2*n2,2*n1) + k1(4,2);
    K(2*n2,2*n2-1) = K(2*n2,2*n2-1) + k1(4,3);
    K(2*n2,2*n2) = K(2*n2,2*n2) + k1(4,4);
    K(2*n2,2*n3-1) = K(2*n2,2*n3-1) + k1(4,5);
    K(2*n2,2*n3) = K(2*n2,2*n3) + k1(4,6);
    K(2*n2,2*n4-1) = K(2*n2,2*n4-1) + k1(4,7);
    K(2*n2,2*n4) = K(2*n2,2*n4) + k1(4,8);

    K(2*n3-1,2*n1-1) = K(2*n3-1,2*n1-1) + k1(5,1);
    K(2*n3-1,2*n1) = K(2*n3-1,2*n1) + k1(5,2);
    K(2*n3-1,2*n2-1) = K(2*n3-1,2*n2-1) + k1(5,3);
    K(2*n3-1,2*n2) = K(2*n3-1,2*n2) + k1(5,4);
    K(2*n3-1,2*n3-1) = K(2*n3-1,2*n3-1) + k1(5,5);
    K(2*n3-1,2*n3) = K(2*n3-1,2*n3) + k1(5,6);
    K(2*n3-1,2*n4-1) = K(2*n3-1,2*n4-1) + k1(5,7);
    K(2*n3-1,2*n4) = K(2*n3-1,2*n4) + k1(5,8);

    K(2*n3,2*n1-1) = K(2*n3,2*n1-1) + k1(6,1);
    K(2*n3,2*n1) = K(2*n3,2*n1) + k1(6,2);
    K(2*n3,2*n2-1) = K(2*n3,2*n2-1) + k1(6,3);
    K(2*n3,2*n2) = K(2*n3,2*n2) + k1(6,4);
    K(2*n3,2*n3-1) = K(2*n3,2*n3-1) + k1(6,5);
    K(2*n3,2*n3) = K(2*n3,2*n3) + k1(6,6);
    K(2*n3,2*n4-1) = K(2*n3,2*n4-1) + k1(6,7);
    K(2*n3,2*n4) = K(2*n3,2*n4) + k1(6,8);

    K(2*n4-1,2*n1-1) = K(2*n4-1,2*n1-1) + k1(7,1);
    K(2*n4-1,2*n1) = K(2*n4-1,2*n1) + k1(7,2);
    K(2*n4-1,2*n2-1) = K(2*n4-1,2*n2-1) + k1(7,3);
    K(2*n4-1,2*n2) = K(2*n4-1,2*n2) + k1(7,4);
    K(2*n4-1,2*n3-1) = K(2*n4-1,2*n3-1) + k1(7,5);
    K(2*n4-1,2*n3) = K(2*n4-1,2*n3) + k1(7,6);
    K(2*n4-1,2*n4-1) = K(2*n4-1,2*n4-1) + k1(7,7);
    K(2*n4-1,2*n4) = K(2*n4-1,2*n4) + k1(7,8);

    K(2*n4,2*n1-1) = K(2*n4,2*n1-1) + k1(8,1);
    K(2*n4,2*n1) = K(2*n4,2*n1) + k1(8,2);
    K(2*n4,2*n2-1) = K(2*n4,2*n2-1) + k1(8,3);
    K(2*n4,2*n2) = K(2*n4,2*n2) + k1(8,4);
    K(2*n4,2*n3-1) = K(2*n4,2*n3-1) + k1(8,5);
    K(2*n4,2*n3) = K(2*n4,2*n3) + k1(8,6);
    K(2*n4,2*n4-1) = K(2*n4,2*n4-1) + k1(8,7);
    K(2*n4,2*n4) = K(2*n4,2*n4) + k1(8,8);


    if (5 <= nnod) then

        K(2*n1-1,2*n5-1) = K(2*n1-1,2*n5-1) + k1(1,9);
        K(2*n1-1,2*n5) = K(2*n1-1,2*n5) + k1(1,10);

        K(2*n1,2*n5-1) = K(2*n1,2*n5-1) + k1(2,9);
        K(2*n1,2*n5) = K(2*n1,2*n5) + k1(2,10);

        K(2*n2-1,2*n5-1) = K(2*n2-1,2*n5-1) + k1(3,9);
        K(2*n2-1,2*n5) = K(2*n2-1,2*n5) + k1(3,10);

        K(2*n2,2*n5-1) = K(2*n2,2*n5-1) + k1(4,9);
        K(2*n2,2*n5) = K(2*n2,2*n5) + k1(4,10);

        K(2*n3-1,2*n5-1) = K(2*n3-1,2*n5-1) + k1(5,9);
        K(2*n3-1,2*n5) = K(2*n3-1,2*n5) + k1(5,10);

        K(2*n3,2*n5-1) = K(2*n3,2*n5-1) + k1(6,9);
        K(2*n3,2*n5) = K(2*n3,2*n5) + k1(6,10);

        K(2*n4-1,2*n5-1) = K(2*n4-1,2*n5-1) + k1(7,9);
        K(2*n4-1,2*n5) = K(2*n4-1,2*n5) + k1(7,10);

        K(2*n4,2*n5-1) = K(2*n4,2*n5-1) + k1(8,9);
        K(2*n4,2*n5) = K(2*n4,2*n5) + k1(8,10);

        K(2*n5-1,2*n1-1) = K(2*n5-1,2*n1-1) + k1(9,1);
        K(2*n5-1,2*n1) = K(2*n5-1,2*n1) + k1(9,2);
        K(2*n5-1,2*n2-1) = K(2*n5-1,2*n2-1) + k1(9,3);
        K(2*n5-1,2*n2) = K(2*n5-1,2*n2) + k1(9,4);
        K(2*n5-1,2*n3-1) = K(2*n5-1,2*n3-1) + k1(9,5);
        K(2*n5-1,2*n3) = K(2*n5-1,2*n3) + k1(9,6);
        K(2*n5-1,2*n4-1) = K(2*n5-1,2*n4-1) + k1(9,7);
        K(2*n5-1,2*n4) = K(2*n5-1,2*n4) + k1(9,8);
        K(2*n5-1,2*n5-1) = K(2*n5-1,2*n5-1) + k1(9,9);
        K(2*n5-1,2*n5) = K(2*n5-1,2*n5) + k1(9,10);

        K(2*n5,2*n1-1) = K(2*n5,2*n1-1) + k1(10,1);
        K(2*n5,2*n1) = K(2*n5,2*n1) + k1(10,2);
        K(2*n5,2*n2-1) = K(2*n5,2*n2-1) + k1(10,3);
        K(2*n5,2*n2) = K(2*n5,2*n2) + k1(10,4);
        K(2*n5,2*n3-1) = K(2*n5,2*n3-1) + k1(10,5);
        K(2*n5,2*n3) = K(2*n5,2*n3) + k1(10,6);
        K(2*n5,2*n4-1) = K(2*n5,2*n4-1) + k1(10,7);
        K(2*n5,2*n4) = K(2*n5,2*n4) + k1(10,8);
        K(2*n5,2*n5-1) = K(2*n5,2*n5-1) + k1(10,9);
        K(2*n5,2*n5) = K(2*n5,2*n5) + k1(10,10);

        if (6 <= nnod) then

            K(2*n1-1,2*n6-1) = K(2*n1-1,2*n6-1) + k1(1,11);
            K(2*n1-1,2*n6) = K(2*n1-1,2*n6) + k1(1,12);

            K(2*n1,2*n6-1) = K(2*n1,2*n6-1) + k1(2,11);
            K(2*n1,2*n6) = K(2*n1,2*n6) + k1(2,12);

            K(2*n2-1,2*n6-1) = K(2*n2-1,2*n6-1) + k1(3,11);
            K(2*n2-1,2*n6) = K(2*n2-1,2*n6) + k1(3,12);

            K(2*n2,2*n6-1) = K(2*n2,2*n6-1) + k1(4,11);
            K(2*n2,2*n6) = K(2*n2,2*n6) + k1(4,12);

            K(2*n3-1,2*n6-1) = K(2*n3-1,2*n6-1) + k1(5,11);
            K(2*n3-1,2*n6) = K(2*n3-1,2*n6) + k1(5,12);

            K(2*n3,2*n6-1) = K(2*n3,2*n6-1) + k1(6,11);
            K(2*n3,2*n6) = K(2*n3,2*n6) + k1(6,12);

            K(2*n4-1,2*n6-1) = K(2*n4-1,2*n6-1) + k1(7,11);
            K(2*n4-1,2*n6) = K(2*n4-1,2*n6) + k1(7,12);

            K(2*n4,2*n6-1) = K(2*n4,2*n6-1) + k1(8,11);
            K(2*n4,2*n6) = K(2*n4,2*n6) + k1(8,12);

            K(2*n5-1,2*n6-1) = K(2*n5-1,2*n6-1) + k1(9,11);
            K(2*n5-1,2*n6) = K(2*n5-1,2*n6) + k1(9,12);

            K(2*n5,2*n6-1) = K(2*n5,2*n6-1) + k1(10,11);
            K(2*n5,2*n6) = K(2*n5,2*n6) + k1(10,12);


            K(2*n6-1,2*n1-1) = K(2*n6-1,2*n1-1) + k1(11,1);
            K(2*n6-1,2*n1) = K(2*n6-1,2*n1) + k1(11,2);
            K(2*n6-1,2*n2-1) = K(2*n6-1,2*n2-1) + k1(11,3);
            K(2*n6-1,2*n2) = K(2*n6-1,2*n2) + k1(11,4);
            K(2*n6-1,2*n3-1) = K(2*n6-1,2*n3-1) + k1(11,5);
            K(2*n6-1,2*n3) = K(2*n6-1,2*n3) + k1(11,6);
            K(2*n6-1,2*n4-1) = K(2*n6-1,2*n4-1) + k1(11,7);
            K(2*n6-1,2*n4) = K(2*n6-1,2*n4) + k1(11,8);
            K(2*n6-1,2*n5-1) = K(2*n6-1,2*n5-1) + k1(11,9);
            K(2*n6-1,2*n5) = K(2*n6-1,2*n5) + k1(11,10);
            K(2*n6-1,2*n6-1) = K(2*n6-1,2*n6-1) + k1(11,11);
            K(2*n6-1,2*n6) = K(2*n6-1,2*n6) + k1(11,12);

            K(2*n6,2*n1-1) = K(2*n6,2*n1-1) + k1(12,1);
            K(2*n6,2*n1) = K(2*n6,2*n1) + k1(12,2);
            K(2*n6,2*n2-1) = K(2*n6,2*n2-1) + k1(12,3);
            K(2*n6,2*n2) = K(2*n6,2*n2) + k1(12,4);
            K(2*n6,2*n3-1) = K(2*n6,2*n3-1) + k1(12,5);
            K(2*n6,2*n3) = K(2*n6,2*n3) + k1(12,6);
            K(2*n6,2*n4-1) = K(2*n6,2*n4-1) + k1(12,7);
            K(2*n6,2*n4) = K(2*n6,2*n4) + k1(12,8);
            K(2*n6,2*n5-1) = K(2*n6,2*n5-1) + k1(12,9);
            K(2*n6,2*n5) = K(2*n6,2*n5) + k1(12,10);
            K(2*n6,2*n6-1) = K(2*n6,2*n6-1) + k1(12,11);
            K(2*n6,2*n6) = K(2*n6,2*n6) + k1(12,12);

            if (7 <= nnod) then

                K(2*n1-1,2*n7-1) = K(2*n1-1,2*n7-1) + k1(1,13);
                K(2*n1-1,2*n7) = K(2*n1-1,2*n7) + k1(1,14);

                K(2*n1,2*n7-1) = K(2*n1,2*n7-1) + k1(2,13);
                K(2*n1,2*n7) = K(2*n1,2*n7) + k1(2,14);

                K(2*n2-1,2*n7-1) = K(2*n2-1,2*n7-1) + k1(3,13);
                K(2*n2-1,2*n7) = K(2*n2-1,2*n7) + k1(3,14);

                K(2*n2,2*n7-1) = K(2*n2,2*n7-1) + k1(4,13);
                K(2*n2,2*n7) = K(2*n2,2*n7) + k1(4,14);

                K(2*n3-1,2*n7-1) = K(2*n3-1,2*n7-1) + k1(5,13);
                K(2*n3-1,2*n7) = K(2*n3-1,2*n7) + k1(5,14);

                K(2*n3,2*n7-1) = K(2*n3,2*n7-1) + k1(6,13);
                K(2*n3,2*n7) = K(2*n3,2*n7) + k1(6,14);

                K(2*n4-1,2*n7-1) = K(2*n4-1,2*n7-1) + k1(7,13);
                K(2*n4-1,2*n7) = K(2*n4-1,2*n7) + k1(7,14);

                K(2*n4,2*n7-1) = K(2*n4,2*n7-1) + k1(8,13);
                K(2*n4,2*n7) = K(2*n4,2*n7) + k1(8,14);

                K(2*n5-1,2*n7-1) = K(2*n5-1,2*n7-1) + k1(9,13);
                K(2*n5-1,2*n7) = K(2*n5-1,2*n7) + k1(9,14);

                K(2*n5,2*n7-1) = K(2*n5,2*n7-1) + k1(10,13);
                K(2*n5,2*n7) = K(2*n5,2*n7) + k1(10,14);

                K(2*n6-1,2*n7-1) = K(2*n6-1,2*n7-1) + k1(11,13);
                K(2*n6-1,2*n7) = K(2*n6-1,2*n7) + k1(11,14);

                K(2*n6,2*n7-1) = K(2*n6,2*n7-1) + k1(12,13);
                K(2*n6,2*n7) = K(2*n6,2*n7) + k1(12,14);

                K(2*n7-1,2*n1-1) = K(2*n7-1,2*n1-1) + k1(13,1);
                K(2*n7-1,2*n1) = K(2*n7-1,2*n1) + k1(13,2);
                K(2*n7-1,2*n2-1) = K(2*n7-1,2*n2-1) + k1(13,3);
                K(2*n7-1,2*n2) = K(2*n7-1,2*n2) + k1(13,4);
                K(2*n7-1,2*n3-1) = K(2*n7-1,2*n3-1) + k1(13,5);
                K(2*n7-1,2*n3) = K(2*n7-1,2*n3) + k1(13,6);
                K(2*n7-1,2*n4-1) = K(2*n7-1,2*n4-1) + k1(13,7);
                K(2*n7-1,2*n4) = K(2*n7-1,2*n4) + k1(13,8);
                K(2*n7-1,2*n5-1) = K(2*n7-1,2*n5-1) + k1(13,9);
                K(2*n7-1,2*n5) = K(2*n7-1,2*n5) + k1(13,10);
                K(2*n7-1,2*n6-1) = K(2*n7-1,2*n6-1) + k1(13,11);
                K(2*n7-1,2*n6) = K(2*n7-1,2*n6) + k1(13,12);
                K(2*n7-1,2*n7-1) = K(2*n7-1,2*n7-1) + k1(13,13);
                K(2*n7-1,2*n7) = K(2*n7-1,2*n7) + k1(13,14);

                K(2*n7,2*n1-1) = K(2*n7,2*n1-1) + k1(14,1);
                K(2*n7,2*n1) = K(2*n7,2*n1) + k1(14,2);
                K(2*n7,2*n2-1) = K(2*n7,2*n2-1) + k1(14,3);
                K(2*n7,2*n2) = K(2*n7,2*n2) + k1(14,4);
                K(2*n7,2*n3-1) = K(2*n7,2*n3-1) + k1(14,5);
                K(2*n7,2*n3) = K(2*n7,2*n3) + k1(14,6);
                K(2*n7,2*n4-1) = K(2*n7,2*n4-1) + k1(14,7);
                K(2*n7,2*n4) = K(2*n7,2*n4) + k1(14,8);
                K(2*n7,2*n5-1) = K(2*n7,2*n5-1) + k1(14,9);
                K(2*n7,2*n5) = K(2*n7,2*n5) + k1(14,10);
                K(2*n7,2*n6-1) = K(2*n7,2*n6-1) + k1(14,11);
                K(2*n7,2*n6) = K(2*n7,2*n6) + k1(14,12);
                K(2*n7,2*n7-1) = K(2*n7,2*n7-1) + k1(14,13);
                K(2*n7,2*n7) = K(2*n7,2*n7) + k1(14,14);

                if (8 <= nnod) then

                    K(2*n1-1,2*n8-1) = K(2*n1-1,2*n8-1) + k1(1,15);
                    K(2*n1-1,2*n8) = K(2*n1-1,2*n8) + k1(1,16);

                    K(2*n1,2*n8-1) = K(2*n1,2*n8-1) + k1(2,15);
                    K(2*n1,2*n8) = K(2*n1,2*n8) + k1(2,16);

                    K(2*n2-1,2*n8-1) = K(2*n2-1,2*n8-1) + k1(3,15);
                    K(2*n2-1,2*n8) = K(2*n2-1,2*n8) + k1(3,16);

                    K(2*n2,2*n8-1) = K(2*n2,2*n8-1) + k1(4,15);
                    K(2*n2,2*n8) = K(2*n2,2*n8) + k1(4,16);

                    K(2*n3-1,2*n8-1) = K(2*n3-1,2*n8-1) + k1(5,15);
                    K(2*n3-1,2*n8) = K(2*n3-1,2*n8) + k1(5,16);

                    K(2*n3,2*n8-1) = K(2*n3,2*n8-1) + k1(6,15);
                    K(2*n3,2*n8) = K(2*n3,2*n8) + k1(6,16);

                    K(2*n4-1,2*n8-1) = K(2*n4-1,2*n8-1) + k1(7,15);
                    K(2*n4-1,2*n8) = K(2*n4-1,2*n8) + k1(7,16);

                    K(2*n4,2*n8-1) = K(2*n4,2*n8-1) + k1(8,15);
                    K(2*n4,2*n8) = K(2*n4,2*n8) + k1(8,16);

                    K(2*n5-1,2*n8-1) = K(2*n5-1,2*n8-1) + k1(9,15);
                    K(2*n5-1,2*n8) = K(2*n5-1,2*n8) + k1(9,16);

                    K(2*n5,2*n8-1) = K(2*n5,2*n8-1) + k1(10,15);
                    K(2*n5,2*n8) = K(2*n5,2*n8) + k1(10,16);

                    K(2*n6-1,2*n8-1) = K(2*n6-1,2*n8-1) + k1(11,15);
                    K(2*n6-1,2*n8) = K(2*n6-1,2*n8) + k1(11,16);

                    K(2*n6,2*n8-1) = K(2*n6,2*n8-1) + k1(12,15);
                    K(2*n6,2*n8) = K(2*n6,2*n8) + k1(12,16);

                    K(2*n7-1,2*n8-1) = K(2*n7-1,2*n8-1) + k1(13,15);
                    K(2*n7-1,2*n8) = K(2*n7-1,2*n8) + k1(13,16);

                    K(2*n7,2*n8-1) = K(2*n7,2*n8-1) + k1(14,15);
                    K(2*n7,2*n8) = K(2*n7,2*n8) + k1(14,16);

                    K(2*n8-1,2*n1-1) = K(2*n8-1,2*n1-1) + k1(15,1);
                    K(2*n8-1,2*n1) = K(2*n8-1,2*n1) + k1(15,2);
                    K(2*n8-1,2*n2-1) = K(2*n8-1,2*n2-1) + k1(15,3);
                    K(2*n8-1,2*n2) = K(2*n8-1,2*n2) + k1(15,4);
                    K(2*n8-1,2*n3-1) = K(2*n8-1,2*n3-1) + k1(15,5);
                    K(2*n8-1,2*n3) = K(2*n8-1,2*n3) + k1(15,6);
                    K(2*n8-1,2*n4-1) = K(2*n8-1,2*n4-1) + k1(15,7);
                    K(2*n8-1,2*n4) = K(2*n8-1,2*n4) + k1(15,8);
                    K(2*n8-1,2*n5-1) = K(2*n8-1,2*n5-1) + k1(15,9);
                    K(2*n8-1,2*n5) = K(2*n8-1,2*n5) + k1(15,10);
                    K(2*n8-1,2*n6-1) = K(2*n8-1,2*n6-1) + k1(15,11);
                    K(2*n8-1,2*n6) = K(2*n8-1,2*n6) + k1(15,12);
                    K(2*n8-1,2*n7-1) = K(2*n8-1,2*n7-1) + k1(15,13);
                    K(2*n8-1,2*n7) = K(2*n8-1,2*n7) + k1(15,14);
                    K(2*n8-1,2*n8-1) = K(2*n8-1,2*n8-1) + k1(15,15);
                    K(2*n8-1,2*n8) = K(2*n8-1,2*n8) + k1(15,16);

                    K(2*n8,2*n1-1) = K(2*n8,2*n1-1) + k1(16,1);
                    K(2*n8,2*n1) = K(2*n8,2*n1) + k1(16,2);
                    K(2*n8,2*n2-1) = K(2*n8,2*n2-1) + k1(16,3);
                    K(2*n8,2*n2) = K(2*n8,2*n2) + k1(16,4);
                    K(2*n8,2*n3-1) = K(2*n8,2*n3-1) + k1(16,5);
                    K(2*n8,2*n3) = K(2*n8,2*n3) + k1(16,6);
                    K(2*n8,2*n4-1) = K(2*n8,2*n4-1) + k1(16,7);
                    K(2*n8,2*n4) = K(2*n8,2*n4) + k1(16,8);
                    K(2*n8,2*n5-1) = K(2*n8,2*n5-1) + k1(16,9);
                    K(2*n8,2*n5) = K(2*n8,2*n5) + k1(16,10);
                    K(2*n8,2*n6-1) = K(2*n8,2*n6-1) + k1(16,11);
                    K(2*n8,2*n6) = K(2*n8,2*n6) + k1(16,12);
                    K(2*n8,2*n7-1) = K(2*n8,2*n7-1) + k1(16,13);
                    K(2*n8,2*n7) = K(2*n8,2*n7) + k1(16,14);
                    K(2*n8,2*n8-1) = K(2*n8,2*n8-1) + k1(16,15);
                    K(2*n8,2*n8) = K(2*n8,2*n8) + k1(16,16);

                    if (9 <= nnod) then

                        K(2*n1-1,2*n9-1) = K(2*n1-1,2*n9-1) + k1(1,17);
                        K(2*n1-1,2*n9) = K(2*n1-1,2*n9) + k1(1,18);

                        K(2*n1,2*n9-1) = K(2*n1,2*n9-1) + k1(2,17);
                        K(2*n1,2*n9) = K(2*n1,2*n9) + k1(2,18);

                        K(2*n2-1,2*n9-1) = K(2*n2-1,2*n9-1) + k1(3,17);
                        K(2*n2-1,2*n9) = K(2*n2-1,2*n9) + k1(3,18);

                        K(2*n2,2*n9-1) = K(2*n2,2*n9-1) + k1(4,17);
                        K(2*n2,2*n9) = K(2*n2,2*n9) + k1(4,18);

                        K(2*n3-1,2*n9-1) = K(2*n3-1,2*n9-1) + k1(5,17);
                        K(2*n3-1,2*n9) = K(2*n3-1,2*n9) + k1(5,18);

                        K(2*n3,2*n9-1) = K(2*n3,2*n9-1) + k1(6,17);
                        K(2*n3,2*n9) = K(2*n3,2*n9) + k1(6,18);

                        K(2*n4-1,2*n9-1) = K(2*n4-1,2*n9-1) + k1(7,17);
                        K(2*n4-1,2*n9) = K(2*n4-1,2*n9) + k1(7,18);

                        K(2*n4,2*n9-1) = K(2*n4,2*n9-1) + k1(8,17);
                        K(2*n4,2*n9) = K(2*n4,2*n9) + k1(8,18);

                        K(2*n5-1,2*n9-1) = K(2*n5-1,2*n9-1) + k1(9,17);
                        K(2*n5-1,2*n9) = K(2*n5-1,2*n9) + k1(9,18);

                        K(2*n5,2*n9-1) = K(2*n5,2*n9-1) + k1(10,17);
                        K(2*n5,2*n9) = K(2*n5,2*n9) + k1(10,18);

                        K(2*n6-1,2*n9-1) = K(2*n6-1,2*n9-1) + k1(11,17);
                        K(2*n6-1,2*n9) = K(2*n6-1,2*n9) + k1(11,18);

                        K(2*n6,2*n9-1) = K(2*n6,2*n9-1) + k1(12,17);
                        K(2*n6,2*n9) = K(2*n6,2*n9) + k1(12,18);

                        K(2*n7-1,2*n9-1) = K(2*n7-1,2*n9-1) + k1(13,17);
                        K(2*n7-1,2*n9) = K(2*n7-1,2*n9) + k1(13,18);

                        K(2*n7,2*n9-1) = K(2*n7,2*n9-1) + k1(14,17);
                        K(2*n7,2*n9) = K(2*n7,2*n9) + k1(14,18);

                        K(2*n8-1,2*n9-1) = K(2*n8-1,2*n9-1) + k1(15,17);
                        K(2*n8-1,2*n9) = K(2*n8-1,2*n9) + k1(15,18);

                        K(2*n8,2*n9-1) = K(2*n8,2*n9-1) + k1(16,17);
                        K(2*n8,2*n9) = K(2*n8,2*n9) + k1(16,18);

                        K(2*n9-1,2*n1-1) = K(2*n9-1,2*n1-1) + k1(17,1);
                        K(2*n9-1,2*n1) = K(2*n9-1,2*n1) + k1(17,2);
                        K(2*n9-1,2*n2-1) = K(2*n9-1,2*n2-1) + k1(17,3);
                        K(2*n9-1,2*n2) = K(2*n9-1,2*n2) + k1(17,4);
                        K(2*n9-1,2*n3-1) = K(2*n9-1,2*n3-1) + k1(17,5);
                        K(2*n9-1,2*n3) = K(2*n9-1,2*n3) + k1(17,6);
                        K(2*n9-1,2*n4-1) = K(2*n9-1,2*n4-1) + k1(17,7);
                        K(2*n9-1,2*n4) = K(2*n9-1,2*n4) + k1(17,8);
                        K(2*n9-1,2*n5-1) = K(2*n9-1,2*n5-1) + k1(17,9);
                        K(2*n9-1,2*n5) = K(2*n9-1,2*n5) + k1(17,10);
                        K(2*n9-1,2*n6-1) = K(2*n9-1,2*n6-1) + k1(17,11);
                        K(2*n9-1,2*n6) = K(2*n9-1,2*n6) + k1(17,12);
                        K(2*n9-1,2*n7-1) = K(2*n9-1,2*n7-1) + k1(17,13);
                        K(2*n9-1,2*n7) = K(2*n9-1,2*n7) + k1(17,14);
                        K(2*n9-1,2*n8-1) = K(2*n9-1,2*n8-1) + k1(17,15);
                        K(2*n9-1,2*n8) = K(2*n9-1,2*n8) + k1(17,16);
                        K(2*n9-1,2*n9-1) = K(2*n9-1,2*n9-1) + k1(17,17);
                        K(2*n9-1,2*n9) = K(2*n9-1,2*n9) + k1(17,18);

                        K(2*n9,2*n1-1) = K(2*n9,2*n1-1) + k1(18,1);
                        K(2*n9,2*n1) = K(2*n9,2*n1) + k1(18,2);
                        K(2*n9,2*n2-1) = K(2*n9,2*n2-1) + k1(18,3);
                        K(2*n9,2*n2) = K(2*n9,2*n2) + k1(18,4);
                        K(2*n9,2*n3-1) = K(2*n9,2*n3-1) + k1(18,5);
                        K(2*n9,2*n3) = K(2*n9,2*n3) + k1(18,6);
                        K(2*n9,2*n4-1) = K(2*n9,2*n4-1) + k1(18,7);
                        K(2*n9,2*n4) = K(2*n9,2*n4) + k1(18,8);
                        K(2*n9,2*n5-1) = K(2*n9,2*n5-1) + k1(18,9);
                        K(2*n9,2*n5) = K(2*n9,2*n5) + k1(18,10);
                        K(2*n9,2*n6-1) = K(2*n9,2*n6-1) + k1(18,11);
                        K(2*n9,2*n6) = K(2*n9,2*n6) + k1(18,12);
                        K(2*n9,2*n7-1) = K(2*n9,2*n7-1) + k1(18,13);
                        K(2*n9,2*n7) = K(2*n9,2*n7) + k1(18,14);
                        K(2*n9,2*n8-1) = K(2*n9,2*n8-1) + k1(18,15);
                        K(2*n9,2*n8) = K(2*n9,2*n8) + k1(18,16);
                        K(2*n9,2*n9-1) = K(2*n9,2*n9-1) + k1(18,17);
                        K(2*n9,2*n9) = K(2*n9,2*n9) + k1(18,18);
                    end if
                end if
            end if
        end if
    end if

    deallocate(k1)

  end do

end subroutine assemble_stiffness_matrix


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                      !!!!!!!!!!!!!
  !!!!!!!!!!!!!   LOCAL STIFFNESS    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                      !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine local_stiffness(noel, nodes, nnd, elements, nel, thickness, ndel, D, k1)

  implicit none

  integer, intent(in) :: noel, nnd, nel, ndel
  real(8), dimension(:, :), intent(in) :: nodes
  integer, dimension(:, :), intent(in) :: elements
  real(8), intent(in) :: thickness
  real(8), dimension(3, 3), intent(in) :: D
  real(8), dimension(ndel, ndel), intent(out) :: k1

  integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, nnod, nnum
  integer :: j, m, mm0, nn0, pp0, mm1, nn1, pp1, mm2, nn2, pp2, mm3, nn3, pp3
  integer :: node_counter(nel,1)
  real(8) :: x1, x2, x3, x4, x5, x6, x7, x8, x9, y1, y2, y3, y4, y5, y6, y7, y8, y9, r, s
  real(8) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, detj
  real(8) :: f1r, f2r, f3r, f4r, f5r, f6r, f7r, f8r, f9r
  real(8) :: f1s, f2s, f3s, f4s, f5s, f6s, f7s, f8s, f9s
  real(8), allocatable, dimension(:,:) :: B, Bt, BTD, dxy, drs, xy, kga
  real(8), dimension(2,2) :: dj, dji

  nnod=ndel/2

  allocate(B(3,ndel))
  allocate(Bt(ndel,3))
  allocate(BTD(ndel,3))
  allocate(drs(nnod,2))
  allocate(dxy(nnod,2))
  allocate(xy(2,nnod))
  allocate(kga(ndel,ndel))

    n1 = elements(noel, 1)
    x1 = nodes(n1, 1)
    y1 = nodes(n1, 2)
    xy(1,1)=x1
    xy(2,1)=y1
    n2 = elements(noel, 2)
    x2 = nodes(n2, 1)
    y2 = nodes(n2, 2)
    xy(1,2)=x2
    xy(2,2)=y2
    n3 = elements(noel, 3)
    x3 = nodes(n3, 1)
    y3 = nodes(n3, 2)
    xy(1,3)=x3
    xy(2,3)=y3
    n4 = elements(noel, 4)
    x4 = nodes(n4, 1)
    y4 = nodes(n4, 2)
    xy(1,4)=x4
    xy(2,4)=y4

    if (5 <= nnod) then
        n5 = elements(noel, 5)
        x5 = nodes(n5, 1)
        y5 = nodes(n5, 2)
        xy(1,5)=x5
        xy(2,5)=y5
        if (6 <= nnod) then
            n6 = elements(noel, 6)
            x6 = nodes(n6, 1)
            y6 = nodes(n6, 2)
            xy(1,6)=x6
            xy(2,6)=y6
            if (7 <= nnod) then
                n7 = elements(noel, 7)
                x7 = nodes(n7, 1)
                y7 = nodes(n7, 2)
                xy(1,7)=x7
                xy(2,7)=y7
                if (8 <= nnod) then
                    n8 = elements(noel, 8)
                    x8 = nodes(n8, 1)
                    y8 = nodes(n8, 2)
                    xy(1,8)=x8
                    xy(2,8)=y8
                    if (9 <= nnod) then
                        n9 = elements(noel, 9)
                        x9 = nodes(n9, 1)
                        y9 = nodes(n9, 2)
                        xy(1,9)=x9
                        xy(2,9)=y9
                    end if
                end if
            end if
        end if
    end if

do j=1,4
    if (j==1) then
        r=-0.5773502692
        s=-0.5773502692
    end if
    if (j==2) then
        r=0.5773502692
        s=-0.5773502692
    end if
    if (j==3) then
        r=-0.5773502692
        s=0.5773502692
    end if
    if (j==4) then
        r=0.5773502692
        s=0.5773502692
    end if

    kga = 0.0
    B = 0.0
    Bt = 0.0
    BTD = 0.0

    !!!!! SHAPE FUNCTIONS AND DERIVATIONS !!!!!

    f1=0.25*(1+r)*(1-s)
    f1r=0.25*(1-s)
    f1s=-0.25*(1+r)
    f2=0.25*(1+r)*(1+s)
    f2r=0.25*(1+s)
    f2s=0.25*(1+r)
    f3=0.25*(1-r)*(1+s)
    f3r=-0.25*(1+s)
    f3s=0.25*(1-r)
    f4=0.25*(1-r)*(1-s)
    f4r=-0.25*(1-s)
    f4s=-0.25*(1-r)
    f5=0.5*(1+r)*(1-s**2)
    f5r=0.5*(1-s**2)
    f5s=-(1+r)*s
    f6=0.5*(1+s)*(1-r**2)
    f6r=-(1+s)*r
    f6s=0.5*(1-r**2)
    f7=0.5*(1-r)*(1-s**2)
    f7r=-0.5*(1-s**2)
    f7s=-(1-r)*s
    f8=0.5*(1-s)*(1-r**2)
    f8r=-(1-s)*r
    f8s=-0.5*(1-r**2)
    f9=(1-r**2)*(1-s**2)
    f9r=-2*r*(1-s**2)
    f9s=-2*s*(1-r**2)

    if (5<=nnod) then
        f1=f1-0.5*f5
        f1r=f1r-0.5*f5r
        f1s=f1s-0.5*f5s
        f2=f2-0.5*f5
        f2r=f2r-0.5*f5r
        f2s=f2s-0.5*f5s
    end if

    if (6<=nnod) then
        f2=f2-0.5*f6
        f2r=f2r-0.5*f6r
        f2s=f2s-0.5*f6s
        f3=f3-0.5*f6
        f3r=f3r-0.5*f6r
        f3s=f3s-0.5*f6s
    end if

    if (7<=nnod) then
        f3=f3-0.5*f7
        f3r=f3r-0.5*f7r
        f3s=f3s-0.5*f7s
        f4=f4-0.5*f7
        f4r=f4r-0.5*f7r
        f4s=f4s-0.5*f7s
    end if

    if (8<=nnod) then
        f1=f1-0.5*f8
        f1r=f1r-0.5*f8r
        f1s=f1s-0.5*f8s
        f4=f4-0.5*f8
        f4r=f4r-0.5*f8r
        f4s=f4s-0.5*f8s
    end if

    if (9<=nnod) then
        f1=f1+0.25*f9
        f1r=f1r+0.25*f9r
        f1s=f1s+0.25*f9s
        f2=f2+0.25*f9
        f2r=f2r+0.25*f9r
        f2s=f2s+0.25*f9s
        f3=f3+0.25*f9
        f3r=f3r+0.25*f9r
        f3s=f3s+0.25*f9s
        f4=f4+0.25*f9
        f4r=f4r+0.25*f9r
        f4s=f4s+0.25*f9s
        f5=f5-0.5*f9
        f5r=f5r-0.5*f9r
        f5s=f5s-0.5*f9s
        f6=f6-0.5*f9
        f6r=f6r-0.5*f9r
        f6s=f6s-0.5*f9s
        f7=f7-0.5*f9
        f7r=f7r-0.5*f9r
        f7s=f7s-0.5*f9s
        f8=f8-0.5*f9
        f8r=f8r-0.5*f9r
        f8s=f8s-0.5*f9s
    end if

    !!!!! JACOBIAN AND DXY !!!!!

    drs(1,1)=f1r
    drs(1,2)=f1s
    drs(2,1)=f2r
    drs(2,2)=f2s
    drs(3,1)=f3r
    drs(3,2)=f3s
    drs(4,1)=f4r
    drs(4,2)=f4s

    if (5<=nnod) then
        drs(5,1)=f5r
        drs(5,2)=f5s
    end if

    if (6<=nnod) then
        drs(6,1)=f6r
        drs(6,2)=f6s
    end if

    if (7<=nnod) then
        drs(7,1)=f7r
        drs(7,2)=f7s
    end if

    if (8<=nnod) then
        drs(8,1)=f8r
        drs(8,2)=f8s
    end if

    if (9<=nnod) then
        drs(9,1)=f9r
        drs(9,2)=f9s
    end if

    mm0 = size(xy,1)
    nn0 = size(xy,2)
    pp0 = size(drs,2)

    call sub_matmul(xy, drs, dj, mm0, nn0, pp0)

    detj=abs(dj(1,1)*dj(2,2)-dj(1,2)*dj(2,1))

    dji(1,1)=dj(2,2)/detj
    dji(2,2)=dj(1,1)/detj
    dji(1,2)=-dj(1,2)/detj
    dji(2,1)=-dj(2,1)/detj

    mm3 = size(drs,1)
    nn3 = size(drs,2)
    pp3 = size(dji,2)

    call sub_matmul(drs, dji, dxy, mm3, nn3, pp3)

    !!!!! B MATRIX !!!!!

    B = 0.0

    do m=1,ndel
        if (mod(m, 2) /= 0) then
            nnum=(m+1)/2
            B(1,m)=dxy(nnum,1)
            B(2,m)=0
            B(3,m)=dxy(nnum,2)
        else
            nnum=m/2
            B(1,m)=0
            B(2,m)=dxy(nnum,2)
            B(3,m)=dxy(nnum,1)
        end if
    end do

    !!!!! K FOR GAUSS POINT !!!!!

    call sub_transpose(B, Bt)

    mm1 = size(Bt,1)
    nn1 = size(Bt,2)
    pp1 = size(D,2)

    call sub_matmul(Bt, D, BTD, mm1, nn1, pp1)

    mm2 = size(BTD,1)
    nn2 = size(BTD,2)
    pp2 = size(B,2)

    call sub_matmul(BTD , B , kga, mm2, nn2, pp2)

    kga = kga * thickness * detj

    k1=k1+kga

end do

    deallocate(B, Bt, BTD, dxy, drs, xy, kga)

end subroutine local_stiffness

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                                             !!!!!!!!!!!!!
  !!!!!!!!!!!!!    BOUNDARY CONDITIONS AND DISPLACEMENTS    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                                             !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bcs_displacements(K, F, U, nbc, bc, nnd)

  implicit none

  real(8), dimension(:,:), intent(in) :: K
  real(8), dimension(:), intent(inout) :: F
  real(8), dimension(:, :), intent(in) :: bc
  integer, intent(in) :: nnd, nbc
  real(8), dimension(:), intent(out) :: U

  integer :: i, j, node, direction, cnt
  real(8) :: value0
  integer, dimension(:), allocatable :: uu_zero
  real(8), dimension(:, :), allocatable :: K_F
  real(8), dimension(:), allocatable :: F_F, uf
  integer :: counter, cnti, cntj, cntii, adof

  adof = 2 * nnd - nbc ! Active degrees of freedom

  allocate(K_F(adof, adof))
  allocate(F_F(adof))
  allocate(uf(adof))
  allocate(uu_zero(nbc))

  U = 0.0d0
  cnt = 0

  do i = 1, nbc
    node = bc(i, 1)
    direction = bc(i, 2)
    value0 = bc(i, 3)

    U(2 * node - (2 - direction)) = value0

    F = F - K(:,2*node-(2-direction)) * value0

    cnt = cnt + 1
    uu_zero(cnt) = 2 * node - (2 - direction)

  end do

  K_F = 0.0d0
  cnti = 1

  do i = 1, 2*nnd
    if (all(uu_zero /= i)) then
      cntj = 1
      do j = 1, 2*nnd
        if (all(uu_zero /= j)) then
          K_F(cnti,cntj) = k(i, j)
          cntj = cntj + 1
        end if
      end do
      cnti = cnti + 1
    end if
  end do

  F_F = 0.0d0
  cntii = 1

  do i = 1, 2*nnd
    if (all(uu_zero /= i)) then
      F_F(cntii) = F(i)
      cntii = cntii + 1
    end if
  end do

  call gauss(K_F, F_F, uf)

  counter = 1
  do i = 1, 2 * nnd
    if (all(uu_zero /= i)) then
      U(i) = uf(counter)
      counter = counter + 1
    end if
  end do

  deallocate(K_F, F_F, uf, uu_zero)

end subroutine bcs_displacements

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!    REACTION FORCES    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine r_forces(K, U, bc, nbc, reaction_forces)
  implicit none
  real(8), intent(in) :: K(:,:)
  real(8), intent(in) :: U(:)
  real(8), intent(in) :: bc(:, :)
  integer, intent(in) :: nbc
  real(8), intent(out) :: reaction_forces(:)
  integer :: i, node0, direction0, mm, nn, pp

  reaction_forces = 0.0d0

  do i = 1, nbc
    node0 = bc(i, 1)
    direction0 = bc(i, 2)

  mm = 1
  nn = size(K,2)
  pp = 1

  call sub_matmul(K(2 * node0 - (2 - direction0), :), U, reaction_forces(2 * node0 - (2 - direction0):1), mm, nn, pp)

  end do

end subroutine r_forces

end program main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!        THE END        !!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
