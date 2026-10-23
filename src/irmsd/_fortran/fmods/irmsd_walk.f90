module irmsd_walk
  !*****************************************************************
  !* Bond-preserving atom assignment by a graph walk from the centre
  !* of each fragment outwards.
  !*
  !* Heavy atoms of unique rank are matched first. The others are
  !* visited in BFS order from the graph centre; each is matched to a
  !* free neighbour (same rank) of its matched neighbour's partner,
  !* such that all of its matched neighbours map onto bonded atoms
  !* (the VF2 feasibility test). Dead ends backtrack. Once the matched
  !* atoms define a frame, a choice between candidates goes to the one
  !* closest after a Kabsch fit on the matched atoms; before that, all
  !* candidates are enumerated (up to maxleaf complete assignments).
  !* Every complete assignment is then polished: H atoms follow their
  !* parent, a Kabsch fit on all atoms gives a frame, and the walk is
  !* repeated in that fixed frame until the assignment is stable. A
  !* stable assignment is then tested against all walks that differ
  !* from it in a single bond-constrained choice; an improvement
  !* restarts the polish.
  !*
  !* Clusters (largest fragment below half of the heavy atoms) are
  !* left to min_rmsd: matching whole molecules is an assignment
  !* problem the LSAP handles, and the walk would only enumerate.
  !*****************************************************************
  use crest_parameters
  use strucrd
  use adjacency,only:adjacency_to_neighbours,bfs_visit
  use ls_rmsd,only:quat_rmsd => rmsd
  implicit none
  private

  public :: walk_assignment

  !> complete assignments enumerated before the frame is defined
  integer,parameter :: maxleaf = 64
  !> walk/Kabsch alternations and accepted single-choice changes
  !> per assignment
  integer,parameter :: maxpolish = 20
  !> matched atoms define a frame once one is this far off the line
  !> through two others (bohr)
  real(wp),parameter :: framethr = 0.3_wp*aatoau

  type :: walk_data
    integer :: n = 0
    integer,allocatable :: rank(:)
    logical,allocatable :: heavy(:)
    !> full neighbour lists of ref (nnr/nbr) and mol (nnm/nbm)
    integer,allocatable :: nnr(:),nbr(:,:),nnm(:),nbm(:,:)
    !> atoms of each rank: rmem(rstart(r):rstart(r+1)-1)
    integer,allocatable :: rstart(:),rmem(:)
    !> heavy atoms to walk, centre outwards
    integer :: nord = 0
    integer,allocatable :: order(:)
    !> heavy atoms of unique rank, matched up front
    logical,allocatable :: fixed(:)
    real(wp),allocatable :: xr(:,:),xm(:,:)
  end type walk_data

contains

  subroutine walk_assignment(ref,mol,rank,aref,amol,perm,rmsdval,io)
    !*******************************************************************
    !* Bond-preserving assignment of mol onto ref.
    !*
    !* ref     : reference structure
    !* mol     : structure in matched order, i.e. mol atom i carries
    !*           rank(i) (as min_rmsd leaves it); coordinates are only
    !*           read
    !* rank    : canonical rank of each atom
    !* aref    : adjacency matrix of ref
    !* amol    : adjacency matrix of mol
    !* perm    : (out) mol atom perm(i) is assigned to ref atom i
    !* rmsdval : (out) RMSD of that assignment (bohr), H atoms assigned
    !*           greedily; the caller refines them
    !* io      : (out) 0 on success, 1 if no assignment was found
    !*******************************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(in) :: rank(:)
    integer,intent(in) :: aref(:,:),amol(:,:)
    integer,intent(out) :: perm(:)
    real(wp),intent(out) :: rmsdval
    integer,intent(out) :: io

    type(walk_data) :: wd
    integer,allocatable :: map(:),inv(:),leaves(:,:),full(:)
    integer :: nleaf,l,iol
    real(wp) :: U(3,3),cm(3),cr(3),r

    io = 1
    rmsdval = huge(rmsdval)
    if (ref%nat < 3.or.mol%nat /= ref%nat) return
    if (any(rank(1:ref%nat) < 1)) return
    if (all(ref%at(1:ref%nat) == 1)) return
    call setup_walk(wd,ref,mol,rank,aref,amol,iol)
    if (iol /= 0) return

    allocate (map(wd%n),inv(wd%n),full(wd%n))
    allocate (leaves(wd%n,maxleaf))
    call start_map(wd,map,inv)
    U = 0.0_wp
    cm = 0.0_wp
    cr = 0.0_wp
    call dfs_walk(wd,map,inv,.false.,U,cm,cr,leaves,nleaf)

!>--- polish every start by alternation only, then the best one fully
    do l = 1,nleaf
      call polish(wd,leaves(:,l),.false.,full,r,iol)
      if (iol /= 0) cycle
      if (r < rmsdval) then
        rmsdval = r
        perm(1:wd%n) = full(:)
        io = 0
      end if
    end do
    if (io /= 0) return
    map(:) = perm(1:wd%n)
    where (.not.wd%heavy) map = 0
    call polish(wd,map,.true.,full,r,iol)
    if (iol == 0.and.r < rmsdval) then
      rmsdval = r
      perm(1:wd%n) = full(:)
    end if
  end subroutine walk_assignment

!========================================================================================!

  subroutine setup_walk(wd,ref,mol,rank,aref,amol,io)
    !****************************************************************
    !* Neighbour lists, rank members and the walk order: per heavy
    !* fragment (largest first) a BFS from its centre, the atom of
    !* smallest eccentricity (ties: smaller rank class, lower rank).
    !****************************************************************
    implicit none
    type(walk_data),intent(inout) :: wd
    type(coord),intent(in) :: ref,mol
    integer,intent(in) :: rank(:)
    integer,intent(in) :: aref(:,:),amol(:,:)
    integer,intent(out) :: io
    integer,allocatable :: nnh(:),nbh(:,:),dist(:),queue(:),ecc(:),frag(:)
    integer,allocatable :: csize(:),centre(:),fsize(:),cnt(:),forder(:)
    integer :: n,i,j,k,l,f,nfrag,nvis,maxr,c,best

    io = 1
    n = ref%nat
    wd%n = n
    wd%rank = rank(1:n)
    wd%heavy = ref%at(1:n) /= 1
    wd%xr = ref%xyz(:,1:n)
    wd%xm = mol%xyz(:,1:n)
    call adjacency_to_neighbours(n,aref,wd%nnr,wd%nbr)
    call adjacency_to_neighbours(n,amol,wd%nnm,wd%nbm)
    if (any(wd%nnr /= wd%nnm)) return

!>--- atoms of each rank
    maxr = maxval(wd%rank)
    allocate (cnt(maxr+1),source=0)
    do i = 1,n
      cnt(wd%rank(i)) = cnt(wd%rank(i))+1
    end do
    allocate (wd%rstart(maxr+1))
    wd%rstart(1) = 1
    do k = 1,maxr
      wd%rstart(k+1) = wd%rstart(k)+cnt(k)
    end do
    allocate (wd%rmem(n))
    cnt(:) = 0
    do i = 1,n
      k = wd%rank(i)
      wd%rmem(wd%rstart(k)+cnt(k)) = i
      cnt(k) = cnt(k)+1
    end do
    allocate (csize(n))
    do i = 1,n
      csize(i) = wd%rstart(wd%rank(i)+1)-wd%rstart(wd%rank(i))
    end do
    wd%fixed = wd%heavy.and.csize == 1

!>--- heavy-atom graph of ref, its fragments and eccentricities
    allocate (nnh(n),source=0)
    allocate (nbh(max(size(wd%nbr,1),1),n),source=0)
    do i = 1,n
      if (.not.wd%heavy(i)) cycle
      do l = 1,wd%nnr(i)
        j = wd%nbr(l,i)
        if (.not.wd%heavy(j)) cycle
        nnh(i) = nnh(i)+1
        nbh(nnh(i),i) = j
      end do
    end do
    allocate (dist(n),queue(n),ecc(n),frag(n))
    dist(:) = -1
    frag(:) = 0
    ecc(:) = 0
    nfrag = 0
    do i = 1,n
      if (.not.wd%heavy(i)) cycle
      call bfs_visit(i,nnh,nbh,dist,queue,nvis)
      ecc(i) = dist(queue(nvis))
      if (frag(i) == 0) then
        nfrag = nfrag+1
        frag(queue(1:nvis)) = nfrag
      end if
      dist(queue(1:nvis)) = -1
    end do

!>--- centre and size of each fragment; clusters are not walked
    allocate (centre(nfrag),source=0)
    allocate (fsize(nfrag),source=0)
    do i = 1,n
      if (.not.wd%heavy(i)) cycle
      f = frag(i)
      fsize(f) = fsize(f)+1
      c = centre(f)
      if (c == 0) then
        centre(f) = i
      else if (ecc(i) < ecc(c)) then
        centre(f) = i
      else if (ecc(i) == ecc(c)) then
        if (csize(i) < csize(c).or. &
        &  (csize(i) == csize(c).and.wd%rank(i) < wd%rank(c))) centre(f) = i
      end if
    end do

    if (2*maxval(fsize) < count(wd%heavy)) return

!>--- fragments by size (largest first), each walked from its centre
    allocate (forder(nfrag))
    do f = 1,nfrag
      forder(f) = f
    end do
    do f = 2,nfrag
      c = forder(f)
      k = f-1
      do while (k >= 1)
        best = forder(k)
        if (fsize(best) >= fsize(c)) exit
        forder(k+1) = best
        k = k-1
      end do
      forder(k+1) = c
    end do
    allocate (wd%order(n))
    wd%nord = 0
    do k = 1,nfrag
      f = forder(k)
      call bfs_visit(centre(f),nnh,nbh,dist,queue,nvis)
      do l = 1,nvis
        if (wd%fixed(queue(l))) cycle
        wd%nord = wd%nord+1
        wd%order(wd%nord) = queue(l)
      end do
      dist(queue(1:nvis)) = -1
    end do
    io = 0
  end subroutine setup_walk

!========================================================================================!

  subroutine start_map(wd,map,inv)
    !*********************************************
    !* Empty assignment plus the unique-rank atoms
    !*********************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(out) :: map(:),inv(:)
    integer :: i
    map(:) = 0
    inv(:) = 0
    do i = 1,wd%n
      if (wd%fixed(i)) then
        map(i) = i
        inv(i) = i
      end if
    end do
  end subroutine start_map

!========================================================================================!

  subroutine dfs_walk(wd,map,inv,fixframe,U,cm,cr,leaves,nleaf, &
  &                   forcepos,forceidx,ndec,decpos,decnum)
    !*****************************************************************
    !* Depth-first walk over wd%order, extending the assignment map
    !* (ref -> mol; inv is mol -> ref) to all heavy atoms.
    !*
    !* fixframe = .true. : every choice by distance in the given frame
    !*                     (U,cm,cr); stops at the first complete
    !*                     assignment
    !* fixframe = .false.: choices by distance once the matched atoms
    !*                     define a frame, otherwise enumerated; collects
    !*                     up to maxleaf complete assignments
    !*
    !* leaves(:,1:nleaf) : (out) the complete heavy-atom assignments
    !*
    !* Fixed frame only:
    !* forcepos, forceidx : take the forceidx-th nearest candidate at
    !*                      walk position forcepos instead of the nearest
    !* ndec, decpos(1:ndec), decnum(1:ndec) : (out) positions of the
    !*                      assignment with more than one candidate drawn
    !*                      from a matched neighbour, and their number of
    !*                      candidates
    !*****************************************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(inout) :: map(:),inv(:)
    logical,intent(in) :: fixframe
    real(wp),intent(inout) :: U(3,3),cm(3),cr(3)
    integer,intent(out) :: leaves(:,:)
    integer,intent(out) :: nleaf
    integer,intent(in),optional :: forcepos,forceidx
    integer,intent(out),optional :: ndec
    integer,intent(out),optional :: decpos(:),decnum(:)
    integer,allocatable :: cbeg(:),cnum(:),cidx(:),stack(:)
    logical,allocatable :: branch(:),anchored(:)
    integer :: pos,top,x,k,nmapped,nalign,steps,maxsteps,p
    real(wp) :: Uw(3,3),cmw(3),crw(3)
    logical :: framed

    nleaf = 0
    if (present(ndec)) ndec = 0
    allocate (cbeg(wd%nord+1),cnum(wd%nord+1),cidx(wd%nord+1))
    allocate (branch(wd%nord+1),source=.false.)
    allocate (anchored(wd%nord+1),source=.false.)
    allocate (stack(max(4*wd%n,16)))
    top = 0
    nmapped = count(map > 0)
    nalign = -1
    Uw = U
    cmw = cm
    crw = cr
    maxsteps = 50*wd%nord+10000
    steps = 0
    pos = 1

    do
      steps = steps+1
      if (steps > maxsteps) exit

      if (pos > wd%nord) then
!>--- complete assignment
        nleaf = nleaf+1
        leaves(:,nleaf) = map(:)
        if (present(ndec)) then
          do p = 1,wd%nord
            if (cnum(p) < 2.or..not.anchored(p)) cycle
            ndec = ndec+1
            decpos(ndec) = p
            decnum(ndec) = cnum(p)
          end do
        end if
        if (fixframe.or.nleaf >= size(leaves,2)) exit
        call backtrack(.true.)
        if (pos < 1) exit
        cycle
      end if

      x = wd%order(pos)
      cbeg(pos) = top+1
      call candidates(x,k,anchored(pos))
      cnum(pos) = k
      if (k == 0) then
        call backtrack(.false.)
        if (pos < 1) exit
        cycle
      end if

!>--- several candidates: nearest first in a defined frame, else branch
      branch(pos) = .false.
      if (k > 1) then
        if (fixframe) then
          framed = .true.
        else
          framed = defines_frame(wd,map)
          if (framed.and.nmapped /= nalign) then
            call fit_mapped(wd,map,Uw,cmw,crw)
            nalign = nmapped
          end if
        end if
        if (framed) then
          call sort_by_distance(x,stack(cbeg(pos):top))
        else
          branch(pos) = .true.
        end if
      end if
      cidx(pos) = 1
      if (present(forcepos)) then
        if (pos == forcepos) cidx(pos) = min(forceidx,k)
      end if
      call assign(x,stack(cbeg(pos)+cidx(pos)-1))
      pos = pos+1
    end do

  contains

    subroutine candidates(xx,kk,anch)
      !> feasible mol partners of ref atom xx, pushed onto stack; anch
      !> tells whether they come from a matched neighbour
      integer,intent(in) :: xx
      integer,intent(out) :: kk
      logical,intent(out) :: anch
      integer :: anchor,ll,y,r1,r2,jj
      kk = 0
      anchor = 0
      do ll = 1,wd%nnr(xx)
        if (map(wd%nbr(ll,xx)) > 0) then
          anchor = wd%nbr(ll,xx)
          exit
        end if
      end do
      anch = anchor > 0
      if (anch) then
        do ll = 1,wd%nnm(map(anchor))
          y = wd%nbm(ll,map(anchor))
          if (feasible(xx,y)) call push(y,kk)
        end do
      else
        r1 = wd%rstart(wd%rank(xx))
        r2 = wd%rstart(wd%rank(xx)+1)-1
        do jj = r1,r2
          y = wd%rmem(jj)
          if (feasible(xx,y)) call push(y,kk)
        end do
      end if
    end subroutine candidates

    subroutine push(yy,kk)
      !> put candidate yy on the stack (grown as needed)
      integer,intent(in) :: yy
      integer,intent(inout) :: kk
      integer,allocatable :: tmp(:)
      if (top == size(stack)) then
        allocate (tmp(2*size(stack)))
        tmp(1:top) = stack(1:top)
        call move_alloc(tmp,stack)
      end if
      top = top+1
      stack(top) = yy
      kk = kk+1
    end subroutine push

    logical function feasible(xx,y)
      !> same rank and degree, free, and the matched neighbours of xx
      !> map onto exactly the matched neighbours of y
      integer,intent(in) :: xx,y
      integer :: ll,z,nx,ny
      feasible = .false.
      if (inv(y) /= 0) return
      if (wd%rank(y) /= wd%rank(xx)) return
      if (wd%nnm(y) /= wd%nnr(xx)) return
      nx = 0
      do ll = 1,wd%nnr(xx)
        z = map(wd%nbr(ll,xx))
        if (z == 0) cycle
        nx = nx+1
        if (.not.any(wd%nbm(1:wd%nnm(y),y) == z)) return
      end do
      ny = 0
      do ll = 1,wd%nnm(y)
        if (inv(wd%nbm(ll,y)) /= 0) ny = ny+1
      end do
      feasible = (nx == ny)
    end function feasible

    subroutine sort_by_distance(xx,c)
      !> insertion sort of the candidates c by distance to ref atom xx
      integer,intent(in) :: xx
      integer,intent(inout) :: c(:)
      real(wp) :: d(size(c)),dt
      integer :: a,b,t
      do a = 1,size(c)
        d(a) = sum((matmul(Uw,wd%xm(:,c(a))-cmw)+crw-wd%xr(:,xx))**2)
      end do
      do a = 2,size(c)
        t = c(a)
        dt = d(a)
        b = a-1
        do while (b >= 1)
          if (d(b) <= dt) exit
          c(b+1) = c(b)
          d(b+1) = d(b)
          b = b-1
        end do
        c(b+1) = t
        d(b+1) = dt
      end do
    end subroutine sort_by_distance

    subroutine assign(xx,y)
      integer,intent(in) :: xx,y
      map(xx) = y
      inv(y) = xx
      nmapped = nmapped+1
    end subroutine assign

    subroutine backtrack(branch_only)
      !> undo positions until one has an untried candidate (for
      !> branch_only, only enumerated choices count); pos < 1 if none
      logical,intent(in) :: branch_only
      integer :: xx
      do
        pos = pos-1
        if (pos < 1) return
        xx = wd%order(pos)
        inv(map(xx)) = 0
        map(xx) = 0
        nmapped = nmapped-1
        if (cidx(pos) < cnum(pos).and.(branch(pos).or..not.branch_only)) then
          cidx(pos) = cidx(pos)+1
          call assign(xx,stack(cbeg(pos)+cidx(pos)-1))
          pos = pos+1
          return
        end if
        top = cbeg(pos)-1
      end do
    end subroutine backtrack

  end subroutine dfs_walk

!========================================================================================!

  logical function defines_frame(wd,map)
    !****************************************************************
    !* .true. if the matched atoms fix an orientation: at least three
    !* of them, not on one line (one is > framethr off the line
    !* through the first and the one farthest from it).
    !****************************************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(in) :: map(:)
    integer :: i,i0,i1
    real(wp) :: dmax,d,e(3),v(3)
    defines_frame = .false.
    if (count(map > 0) < 3) return
    i0 = 0
    i1 = 0
    dmax = -1.0_wp
    do i = 1,wd%n
      if (map(i) == 0) cycle
      if (i0 == 0) then
        i0 = i
        cycle
      end if
      d = sum((wd%xr(:,i)-wd%xr(:,i0))**2)
      if (d > dmax) then
        dmax = d
        i1 = i
      end if
    end do
    if (dmax <= 0.0_wp) return
    e = (wd%xr(:,i1)-wd%xr(:,i0))/sqrt(dmax)
    do i = 1,wd%n
      if (map(i) == 0) cycle
      v = wd%xr(:,i)-wd%xr(:,i0)
      v = v-dot_product(v,e)*e
      if (sum(v**2) > framethr**2) then
        defines_frame = .true.
        return
      end if
    end do
  end function defines_frame

!========================================================================================!

  subroutine fit_mapped(wd,map,U,cm,cr,r)
    !****************************************************************
    !* Kabsch fit of the matched mol atoms onto their ref partners:
    !* matmul(U,xm-cm)+cr is the fitted mol position. r is the RMSD.
    !****************************************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(in) :: map(:)
    real(wp),intent(out) :: U(3,3),cm(3),cr(3)
    real(wp),intent(out),optional :: r
    real(wp),allocatable :: X(:,:),Y(:,:)
    real(wp) :: err,gdum(3,1)
    integer :: i,k
    k = count(map > 0)
    if (k < 1) then
      U = reshape([1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp,0.0_wp,0.0_wp,0.0_wp,1.0_wp],[3,3])
      cm = 0.0_wp
      cr = 0.0_wp
      if (present(r)) r = huge(1.0_wp)
      return
    end if
    allocate (X(3,k),Y(3,k))
    k = 0
    do i = 1,wd%n
      if (map(i) == 0) cycle
      k = k+1
      X(:,k) = wd%xm(:,map(i))
      Y(:,k) = wd%xr(:,i)
    end do
    call quat_rmsd(k,X,Y,1,U,cm,cr,err,.false.,gdum)
    if (present(r)) r = err
  end subroutine fit_mapped

!========================================================================================!

  subroutine assign_light(wd,map,U,cm,cr,io)
    !****************************************************************
    !* Extend a heavy-atom assignment to the H atoms: each goes to the
    !* nearest free H of the same rank on its parent's partner (or, if
    !* not bound to a heavy atom, among all free H of that rank).
    !****************************************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(inout) :: map(:)
    real(wp),intent(in) :: U(3,3),cm(3),cr(3)
    integer,intent(out) :: io
    logical,allocatable :: used(:)
    integer :: x,l,p,y,best,jj
    real(wp) :: d,dbest
    io = 1
    allocate (used(wd%n),source=.false.)
    do x = 1,wd%n
      if (map(x) > 0) used(map(x)) = .true.
    end do
    do x = 1,wd%n
      if (map(x) > 0) cycle
      p = 0
      do l = 1,wd%nnr(x)
        if (wd%heavy(wd%nbr(l,x))) then
          p = wd%nbr(l,x)
          exit
        end if
      end do
      best = 0
      dbest = huge(dbest)
      if (p > 0) then
        do l = 1,wd%nnm(map(p))
          y = wd%nbm(l,map(p))
          call consider(y)
        end do
      else
        do jj = wd%rstart(wd%rank(x)),wd%rstart(wd%rank(x)+1)-1
          call consider(wd%rmem(jj))
        end do
      end if
      if (best == 0) return
      map(x) = best
      used(best) = .true.
    end do
    io = 0
  contains
    subroutine consider(yy)
      integer,intent(in) :: yy
      if (used(yy).or.wd%rank(yy) /= wd%rank(x)) return
      d = sum((matmul(U,wd%xm(:,yy)-cm)+cr-wd%xr(:,x))**2)
      if (d < dbest) then
        dbest = d
        best = yy
      end if
    end subroutine consider
  end subroutine assign_light

!========================================================================================!

  subroutine polish(wd,heavymap,search,full,r,io)
    !****************************************************************
    !* Alternate a Kabsch fit on all atoms with a re-walk of the heavy
    !* atoms in that fixed frame, from the heavy-atom assignment
    !* heavymap, until the assignment no longer changes. With search,
    !* then try the walks that differ in one choice; an improvement
    !* restarts this.
    !* full : (out) best complete assignment found (ref -> mol)
    !* r    : (out) its RMSD
    !****************************************************************
    implicit none
    type(walk_data),intent(in) :: wd
    integer,intent(in) :: heavymap(:)
    logical,intent(in) :: search
    integer,intent(out) :: full(:)
    real(wp),intent(out) :: r
    integer,intent(out) :: io
    integer,allocatable :: map(:),inv(:),cur(:),leaf(:,:),trial(:)
    integer,allocatable :: decpos(:),decnum(:)
    real(wp) :: U(3,3),cm(3),cr(3),Ut(3,3),cmt(3),crt(3),rit,rt
    integer :: it,nleaf,iol,ndec,d,a
    logical :: improved

    io = 1
    r = huge(r)
    allocate (map(wd%n),inv(wd%n),cur(wd%n),trial(wd%n),leaf(wd%n,1))
    allocate (decpos(wd%nord),decnum(wd%nord))
    cur(:) = heavymap(:)
    call fit_mapped(wd,cur,U,cm,cr)
    do it = 1,maxpolish
      call complete(cur,U,cm,cr,map,rit,iol)
      if (iol /= 0) exit
      if (rit < r) then
        r = rit
        full(:) = map(:)
        io = 0
      end if

!>--- re-walk in the frame of the complete assignment
      call start_map(wd,map,inv)
      call dfs_walk(wd,map,inv,.true.,U,cm,cr,leaf,nleaf, &
      &             ndec=ndec,decpos=decpos,decnum=decnum)
      if (nleaf < 1) exit
      if (any(leaf(:,1) /= cur)) then
        cur(:) = leaf(:,1)
        cycle
      end if

!>--- stable: try the walks that differ in a single choice
      if (.not.search) exit
      improved = .false.
      DECISIONS: do d = 1,ndec
        do a = 2,decnum(d)
          call start_map(wd,map,inv)
          call dfs_walk(wd,map,inv,.true.,U,cm,cr,leaf,nleaf, &
          &             forcepos=decpos(d),forceidx=a)
          if (nleaf < 1) cycle
          Ut = U
          cmt = cm
          crt = cr
          call complete(leaf(:,1),Ut,cmt,crt,trial,rt,iol)
          if (iol /= 0) cycle
          if (rt < r-1.0e-10_wp) then
            cur(:) = leaf(:,1)
            U = Ut
            cm = cmt
            cr = crt
            improved = .true.
            exit DECISIONS
          end if
        end do
      end do DECISIONS
      if (.not.improved) exit
    end do

  contains

    subroutine complete(hmap,Uc,cmc,crc,fmap,rc,ioc)
      !> H atoms by the frame (Uc,cmc,crc), then a Kabsch fit on all
      !> atoms, which replaces the frame; rc is the RMSD of fmap
      integer,intent(in) :: hmap(:)
      real(wp),intent(inout) :: Uc(3,3),cmc(3),crc(3)
      integer,intent(out) :: fmap(:)
      real(wp),intent(out) :: rc
      integer,intent(out) :: ioc
      fmap(:) = hmap(:)
      call assign_light(wd,fmap,Uc,cmc,crc,ioc)
      if (ioc /= 0) return
      call fit_mapped(wd,fmap,Uc,cmc,crc,rc)
    end subroutine complete

  end subroutine polish

end module irmsd_walk
