
module canonical_mod
!*************************************************************************
!* Implementation of different algorithms for determining atom identities
!*
!* A) Implementation of the CANGEN algorithm by Weininger et al.
!* D.Weininger et al., J. Chem. Inf. Comput. Sci., 1989, 29, 97-101.
!* doi.org/10.1021/ci00062a008
!* The algorithm is used as the backend implementation to
!* determine the unqiue atom sequence in canonical SMILES.
!* It is pretty useful for all kinds of stuff, e.g.
!* for determining tautomers, or trying to restore atom order
!*
!* B) A custom implementation based on the all-pair-shortest-path
!* of the molecular graph, inspired by the Morgan algorithm
!*
!* C) a Helper algorithm to determine the (relative) stereo orientation
!* after a given canonical atom ranking has been generated
!*************************************************************************
  use crest_parameters
  use strucrd
  use adjacency
  use geo
  use utilities,only:first_primes,sorted_order,isort
  implicit none
  private

  public :: canonical_sorter
  type :: canonical_sorter
    !> system mapping
    integer :: nat = 0    !> total number of atoms
    integer :: hatms = 0  !> number of heavy atoms (>H)
    integer :: nfrag = 0  !> number of molecules/subgraphs in structure
    integer,allocatable :: nmap(:)  !> map atom to heavy-atom list (H->0)
    integer,allocatable :: hmap(:)  !> map heavy-atom to full atom order
    integer :: maxnei = 0
    integer,allocatable :: neigh(:,:)  !> neighbour list (neigh(j,i) = j-th neighbour of atom i)
    integer,allocatable :: nhnei(:)    !> number of heavy neighbours of each heavy atom
    integer,allocatable :: hneigh(:,:) !> heavy neighbours, in heavy-atom indices (never H)

    !> the important bit for the algorithm
    integer(int64),allocatable :: invariants(:)
    integer,allocatable :: invariants0(:)
    logical,allocatable :: ovf(:)      !> prime product overflowed int64
    integer,allocatable :: prime(:)
    integer,allocatable :: primetab(:) !> the first hatms primes
    integer,allocatable :: rank(:)        !> That's what we are after

    !> workspace helpers
    integer,allocatable :: newrank(:)

  contains
    procedure :: deallocate => deallocate_canonical_sorter
    procedure :: shrink => shrink_canonical_sorter
    procedure :: init => init_canonical_sorter
    procedure :: init_graph => init_canonical_graph
    procedure :: init_connect => init_canonical_sorter_connect
    procedure :: update_ranks
    procedure :: update_invariants
    procedure :: iterate
    procedure :: rankprint
    procedure :: stereo => analyze_stereo
    procedure :: hasstereo => has_stereo
    procedure :: compare => compare_canonical_sorter
    procedure :: add_h_ranks
    procedure :: split_pair_ranks
    procedure :: same_graph
  end type canonical_sorter

  logical,parameter :: debug = .false.

!========================================================================================!
!========================================================================================!
contains  !> MODULE PROCEDURES START HERE
!========================================================================================!
!========================================================================================!

  subroutine deallocate_canonical_sorter(self)
!********************************************
!* reset canonical_sorter entirely
!********************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    self%nat = 0
    self%hatms = 0
    self%nfrag = 0
    if (allocated(self%nmap)) deallocate (self%nmap)
    if (allocated(self%hmap)) deallocate (self%hmap)
    if (allocated(self%invariants)) deallocate (self%invariants)
    if (allocated(self%invariants0)) deallocate (self%invariants0)
    if (allocated(self%rank)) deallocate (self%rank)
    if (allocated(self%prime)) deallocate (self%prime)
    if (allocated(self%primetab)) deallocate (self%primetab)
    if (allocated(self%ovf)) deallocate (self%ovf)
    if (allocated(self%nhnei)) deallocate (self%nhnei)
    if (allocated(self%hneigh)) deallocate (self%hneigh)
    if (allocated(self%newrank)) deallocate (self%newrank)
    if (allocated(self%neigh)) deallocate (self%neigh)
  end subroutine deallocate_canonical_sorter

!========================================================================================!

  subroutine shrink_canonical_sorter(self)
!***************************************************************
!* Reduce memory consumption of canconical_sorter
!* by deallocating everything except the (original) invariants
!* mappings and the determined ranks
!***************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    if (allocated(self%invariants)) deallocate (self%invariants)
    if (allocated(self%prime)) deallocate (self%prime)
    if (allocated(self%primetab)) deallocate (self%primetab)
    if (allocated(self%ovf)) deallocate (self%ovf)
    if (allocated(self%nhnei)) deallocate (self%nhnei)
    if (allocated(self%hneigh)) deallocate (self%hneigh)
    if (allocated(self%newrank)) deallocate (self%newrank)
    if (allocated(self%neigh)) deallocate (self%neigh)
  end subroutine shrink_canonical_sorter

!========================================================================================!

  subroutine init_canonical_sorter(self,mol,wbo,invtype,heavy)
!*****************************************************************
!* Initializes the canonical_sorter and runs the CANGEN algorithm
!*****************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    real(wp),intent(in),optional :: wbo(mol%nat,mol%nat)
    character(len=*),intent(in),optional :: invtype
    logical,intent(in),optional :: heavy
    integer :: nodes
    integer,allocatable :: Amat(:,:) !> adjacency matrix for FULL molecule
    integer :: counth,countb,countbo
    real(wp) :: countbo2
    integer :: i,j,k,ii,ati
    integer,allocatable :: ichrgs(:)
    logical,allocatable :: hnode(:)
    character(len=:),allocatable :: myinvtype
    logical :: use_icharges,include_H,anyH,nmr_split

!>--- optional argument handling
    if (present(invtype)) then
      myinvtype = invtype
    else
      myinvtype = 'cangen'
    end if
    !> 'apsp+nmr' is the apsp+ channel plus a magnetic-equivalency split (hack)
    nmr_split = (myinvtype .eq. 'apsp+nmr')
    if (present(heavy)) then
      include_H = .not.heavy
    else
      include_H = .false.
    end if
    anyH = any(mol%at(:) .eq. 1)

!>--- build the molecular graph (heavy-atom map, neighbour list, adjacency).
!>    Returns the full-molecule adjacency matrix Amat for the invariant setup.
    if (present(wbo)) then
      call self%init_graph(mol,include_H,Amat,wbo=wbo)
    else
      call self%init_graph(mol,include_H,Amat)
    end if
    nodes = self%nat
    k = self%hatms

!>--- get the first invatiants
    if (allocated(mol%qat)) then
      use_icharges = .true.
      allocate (ichrgs(mol%nat),source=0)
      ichrgs(:) = nint(mol%qat(:))
    else
      use_icharges = .false.
    end if
    select case (myinvtype)
    case ('apsp+','apsp+nmr') !> custom all-pair-shortest-path algo

      allocate (hnode(k))
      hnode(:) = mol%at(self%hmap(:)) .eq. 1
      call get_invariant0_apsp(k,self%nhnei,self%hneigh,hnode,self%invariants0)
      do i = 1,k
        ii = self%hmap(i)
        ati = mol%at(ii)
        counth = 0
        do j = 1,self%maxnei
          if (self%neigh(j,ii) .eq. 0) exit
          if (mol%at(self%neigh(j,ii)) .eq. 1) counth = counth+1 !> count H neighbours
        end do
        self%invariants0(i) = update_invariant0_apsp(self%invariants0(i),ati,counth)
      end do
      self%invariants(:) = int(self%invariants0(:))

    case default !> CANGEN

      if (.not.present(wbo)) then
        error stop 'CANGEN implementation requires wbo matrix as argument'
      end if
      do i = 1,k
        ii = self%hmap(i)
        ati = mol%at(ii)
        counth = 0
        countbo2 = 0.0_wp
        countb = 0
        do j = 1,nodes
          if (Amat(j,ii) .ne. 0) then
            if (mol%at(j) .eq. 1) then
              counth = counth+1 !> count H neighbours
              countbo2 = countbo2-wbo(j,ii) !> but NOT in total bond order
            end if
            countb = countb+1 !> count all neighbours
            !countbo2 = countbo2+wbo(j,ii)  !> sum the total bond order
          end if
          countbo2 = countbo2+wbo(j,ii)  !> sum the total bond order
        end do
        countb = countb-counth          !> only heavy atom connections
        !countbo = nint(countbo2)-counth !> same for number of bonds
        countbo = nint(countbo2)
        if (use_icharges) then
          self%invariants(i) = get_invariant0(ati,countb,countbo,ichrgs(ii),counth)
        else
          self%invariants(i) = get_invariant0(ati,countb,countbo,0,counth)
        end if
        self%invariants0(i) = int(self%invariants(i)) !> back up for later tasks
      end do
    end select

    deallocate (Amat)

    if (debug) then
      call debugprint(self,mol)
    end if
!>--- start assignment
    if (.not.allocated(self%newrank)) allocate (self%newrank(k),source=0) !> workspace
    call self%update_ranks()
    self%rank(:) = self%newrank(:)
    if (debug) then
      call debugprint(self,mol)
    end if
    call self%update_invariants()
    if (debug) then
      call debugprint(self,mol)
    end if
    call self%iterate(mol) !> iterate recursively until ranking doesn't change

!>--- NMR hack: split heavy-atom ranks shared by exactly two atoms 
    if (nmr_split) then
      call self%split_pair_ranks(mol)
    end if

!>--- finally, if required, add H atoms. H nodes normally hold rank 1 alone,
!>    so this renumbering is a shift by one; it stays correct if they do not.
    if (include_H.and.anyH) then
      call renumber_heavy_ranks(self,mol)
      call self%add_h_ranks(mol)
    end if

  end subroutine init_canonical_sorter

!========================================================================================!

  subroutine init_canonical_graph(self,mol,include_H,Amat,wbo)
!*****************************************************************
!* Builds only the molecular graph part of the canonical_sorter:
!* the heavy-atom maps (nmap/hmap), the neighbour list (neigh),
!* the heavy-atom neighbour list (hneigh) and the fragment count. The
!* CANGEN/APSP rank refinement is NOT run here.
!*
!* This is the common graph-setup shared by init_canonical_sorter
!* and has_stereo (which self-fills the graph when only externally
!* supplied ranks are available). The full-molecule adjacency matrix
!* is returned in Amat for any subsequent invariant construction.
!*
!* self      : canonical_sorter whose graph fields get populated
!* mol       : input structure
!* include_H : whether hydrogens are treated as graph nodes/ranks
!* Amat      : (out) full-molecule adjacency matrix
!* wbo       : (optional) bond-order matrix; if absent, connectivity
!*             is derived from coordination numbers
!*****************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    logical,intent(in) :: include_H
    integer,allocatable,intent(out) :: Amat(:,:)
    real(wp),intent(in),optional :: wbo(mol%nat,mol%nat)
    integer :: nodes,i,j,k,l,ii,maxnei
    integer,allocatable :: nnei(:),frag(:)

!>--- all atoms of the full mol. graph are nodes
    nodes = mol%nat

!>--- map to heavy atom-only representation
    k = 0
    do i = 1,mol%nat
      if (mol%at(i) .ne. 1.or.include_h) k = k+1
    end do
    self%nat = nodes
    self%hatms = k
    if (.not.allocated(self%nmap)) allocate (self%nmap(nodes))
    if (.not.allocated(self%hmap)) allocate (self%hmap(k))
    if (.not.allocated(self%invariants)) allocate (self%invariants(k),source=0_int64)
    if (.not.allocated(self%invariants0)) allocate (self%invariants0(k),source=0)
    if (.not.allocated(self%prime)) allocate (self%prime(k),source=2)
    if (.not.allocated(self%rank)) allocate (self%rank(k),source=1)
    if (.not.allocated(self%ovf)) allocate (self%ovf(k),source=.false.)
    self%primetab = first_primes(k)

    if (present(wbo)) then
!>--- get connectivity. Easiest is just via WBO (allocates Amat)
      call wbo2adjacency(nodes,wbo,Amat,0.02_wp)
    else
!>--- connectivity via CN
      call mol%cn_to_adjacency(Amat)
    end if

!>--- documment neighbour list
    maxnei = 0
    do i = 1,mol%nat
      k = count(Amat(:,i) > 0)
      if (k > maxnei) maxnei = k
    end do
    if (debug) write (stdout,*) 'maximum number of neighbours',maxnei
    self%maxnei = maxnei
    if (.not.allocated(self%neigh)) allocate (self%neigh(maxnei,mol%nat),source=0)

!>--- fill rest of self
    k = 0
    do i = 1,nodes
      l = 0
      if (mol%at(i) .ne. 1.or.include_h) then
        k = k+1
        self%nmap(i) = k
        self%hmap(k) = i
      else
        self%nmap(i) = 0
      end if
      do j = 1,nodes
        if (Amat(j,i) > 0) then
          l = l+1
          self%neigh(l,i) = j
        end if
      end do
    end do
    allocate (nnei(nodes),frag(nodes))
    nnei(:) = count(self%neigh(:,:) > 0,dim=1)
    call neighbour_fragments(nodes,nnei,self%neigh,frag,self%nfrag)

!>--- heavy-atom neighbour list; H's are excluded from it, always
    if (allocated(self%nhnei)) deallocate (self%nhnei)
    if (allocated(self%hneigh)) deallocate (self%hneigh)
    allocate (self%nhnei(k),source=0)
    allocate (self%hneigh(max(maxnei,1),k),source=0)
    do i = 1,k
      ii = self%hmap(i)
      if (mol%at(ii) .eq. 1) cycle
      do l = 1,maxnei
        j = self%neigh(l,ii)
        if (j .eq. 0) exit
        if (mol%at(j) .eq. 1) cycle
        self%nhnei(i) = self%nhnei(i)+1
        self%hneigh(self%nhnei(i),i) = self%nmap(j)
      end do
    end do
  end subroutine init_canonical_graph

!========================================================================================!

 subroutine init_canonical_sorter_connect(self,at,wbo,invtype,heavy)
!*****************************************************************
!* Initializes the canonical_sorter provided only atom types and
!* connectivity. No CN calculation etc.
!*****************************************************************
    implicit none
    !> IN-/OUTPUTS
    class(canonical_sorter),intent(inout) :: self
    integer,intent(in) :: at(:)
    real(wp),intent(in) :: wbo(size(at,1),size(at,1))
    character(len=*),intent(in),optional :: invtype
    logical,intent(in),optional :: heavy
    !> LOCAL
    integer :: nat
    type(coord) :: tmpmol

    nat = size(at,1)
    allocate(tmpmol%at(nat))
    tmpmol%at(:) = at(:)
    tmpmol%nat = nat

    call self%init(tmpmol,wbo=wbo,invtype=invtype,heavy=heavy)
end subroutine init_canonical_sorter_connect


!========================================================================================!

  subroutine get_invariant0_apsp(n,nnei,neigh,hnode,inv)
    !*******************************************************************
    !* Initial APSP invariant. Every atom j gets a distance-weighted sum
    !* over the graph, built up shell by shell,
    !*   rinv_d(j) = rinv_{d-1}(j) + sum of rinv_{d-1}(k), d(j,k) = d,
    !* with atoms of other fragments as the last ("infinitely distant")
    !* shell. inv(j) is the rank of rinv(j) in descending order.
    !*
    !* Distances come from one BFS per atom over the heavy-atom graph.
    !* H atoms (hnode) are isolated nodes of that graph and only add the
    !* same constant to every heavy-atom sum, so they are handled in
    !* closed form: they take rank 1 and the heavy atoms start at rank
    !* 2, unless a single-atom heavy fragment ties with them at the top.
    !* The invariants equal those of the former Floyd-Warshall version
    !* as long as the sums are exact (< 2^53). Beyond that, each shell
    !* is summed in ascending order of its values, so the result does
    !* not depend on the atom order and symmetry-equivalent atoms get
    !* bitwise identical sums. All values are rescaled by a power of
    !* two (exact) before they can overflow.
    !*
    !* n     : number of nodes (heavy-atom indices of canonical_sorter)
    !* nnei  : number of heavy neighbours of each node
    !* neigh : heavy neighbours of each node
    !* hnode : .true. for H atoms carried as (isolated) nodes
    !* inv   : (out) the invariant
    !*******************************************************************
    implicit none
    integer,intent(in)  :: n
    integer,intent(in)  :: nnei(n)
    integer,intent(in)  :: neigh(:,:)
    logical,intent(in)  :: hnode(n)
    integer,intent(out) :: inv(n)

    integer,allocatable :: comp(:),csize(:),queue(:),dist(:),ecc(:)
    integer,allocatable :: list(:),lsize(:),idx(:),cur(:),hidx(:)
    integer(int64),allocatable :: lpos(:),spos(:)
    real(wp),allocatable :: rinv(:),tmprinv(:),buf(:),tsum(:)
    integer :: i,l,c,d,m,ncomp,tail,maxecc,nb,shift
    integer(int64) :: nlist,nlay
    real(wp) :: s
    real(wp),parameter :: bigval = 2.0_wp**900

    inv(:) = 1
    if (n < 1) return

!>--- fragments of the heavy-atom graph (comp = 0 for H nodes)
    allocate (comp(n),queue(n),dist(n),ecc(n))
    call neighbour_fragments(n,nnei,neigh,comp,ncomp,mask=.not.hnode)
    if (ncomp == 0) return  !> only H nodes
    allocate (csize(ncomp),source=0)
    do i = 1,n
      if (comp(i) > 0) csize(comp(i)) = csize(comp(i))+1
    end do

!>--- first BFS sweep: eccentricities, to size the shell storage
    dist(:) = -1
    nlist = 0
    nlay = 0
    do i = 1,n
      if (comp(i) == 0) cycle
      call bfs_visit(i,nnei,neigh,dist,queue,tail)
      ecc(i) = dist(queue(tail))
      nlist = nlist+int(tail-1,int64)
      nlay = nlay+int(ecc(i),int64)
      dist(queue(1:tail)) = -1
    end do

!>--- second sweep: for every atom the atoms of each shell (list) and the
!>    shell sizes (lsize), in CSR layout starting at lpos/spos
    allocate (list(max(nlist,1_int64)),lsize(max(nlay,1_int64)))
    allocate (lpos(n),spos(n))
    nlist = 0
    nlay = 0
    do i = 1,n
      if (comp(i) == 0) cycle
      call bfs_visit(i,nnei,neigh,dist,queue,tail)
      lpos(i) = nlist
      spos(i) = nlay
      list(nlist+1:nlist+tail-1) = queue(2:tail)
      lsize(nlay+1:nlay+ecc(i)) = 0
      do l = 2,tail
        d = dist(queue(l))
        lsize(nlay+d) = lsize(nlay+d)+1
      end do
      nlist = nlist+int(tail-1,int64)
      nlay = nlay+int(ecc(i),int64)
      dist(queue(1:tail)) = -1
    end do
    deallocate (queue,dist)

!>--- accumulate the finite shells, all atoms in step
    allocate (rinv(n),tmprinv(n),source=1.0_wp)
    allocate (cur(n),source=0)
    allocate (buf(n),idx(n))
    maxecc = maxval(ecc,mask=comp > 0)
    do d = 1,maxecc
      tmprinv(:) = 0.0_wp
      do i = 1,n
        if (comp(i) == 0) cycle
        if (ecc(i) < d) cycle
        nb = lsize(spos(i)+d)
        do l = 1,nb
          buf(l) = rinv(list(lpos(i)+cur(i)+l))
        end do
        tmprinv(i) = ascending_sum(nb,buf,idx)
        cur(i) = cur(i)+nb
      end do
      rinv(:) = rinv(:)+tmprinv(:)
      if (maxval(rinv,mask=comp > 0) > bigval) rinv(:) = scale(rinv(:),-512)
    end do
    deallocate (list,lsize,lpos,spos,cur)

!>--- last shell: all heavy atoms of the other fragments
    m = count(comp > 0)
    allocate (hidx(m))
    hidx(:) = pack([(i,i=1,n)],comp > 0)
    if (ncomp > 1) then
      allocate (tsum(ncomp))
      buf(1:m) = rinv(hidx)
      call sorted_order(m,buf,idx)
      do c = 1,ncomp
        s = 0.0_wp
        do l = 1,m
          if (comp(hidx(idx(l))) == c) cycle
          s = s+buf(idx(l))
        end do
        tsum(c) = s
      end do
      do l = 1,m
        rinv(hidx(l)) = rinv(hidx(l))+tsum(comp(hidx(l)))
      end do
      deallocate (tsum)
    end if

!>--- dense rank in descending order; H atoms hold the top rank
    shift = 0
    if (any(hnode).and..not.any(csize(:) == 1)) shift = 1
    buf(1:m) = -rinv(hidx)
    call sorted_order(m,buf,idx)
    d = 1
    do l = 1,m
      if (l > 1) then
        if (buf(idx(l)) /= buf(idx(l-1))) d = d+1
      end if
      inv(hidx(idx(l))) = d+shift
    end do
  end subroutine get_invariant0_apsp

  function ascending_sum(n,v,idx) result(s)
    !***********************************************
    !* Sum of v(1:n), added in ascending order so
    !* the result depends only on the multiset.
    !***********************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: v(:)
    integer,intent(inout) :: idx(:)
    real(wp) :: s
    integer :: i
    s = 0.0_wp
    if (n < 1) return
    call sorted_order(n,v,idx)
    do i = 1,n
      s = s+v(idx(i))
    end do
  end function ascending_sum


!========================================================================================!

  function update_invariant0_apsp(inv_in,ati,hneigh) result(inv)
!************************************************************************
!* update invariant for atom by appending atom type and number of H atoms
!*************************************************************************
    implicit none
    integer :: inv
    integer,intent(in) :: inv_in !> previous invariant
    integer,intent(in) :: ati !> atomic number
    integer,intent(in) :: hneigh !> # H atoms bound to atom
    inv = inv_in*10**3
    inv = inv+ati*10
    inv = inv+hneigh
  end function update_invariant0_apsp

!========================================================================================!

  function get_invariant0(ati,nneigh,nbonds,chrg,hneigh) result(inv)
!**********************************************************************
!* Get invariant for atom according to the origin al CANGEN scheme
!**********************************************************************
    implicit none
    integer :: inv
    integer,intent(in) :: ati !> atomic number
    integer,intent(in) :: nneigh !> # neighbours
    integer,intent(in) :: nbonds !> # bonds (sum of wbo)
    integer,intent(in) :: chrg   !> charge on atom
    integer,intent(in) :: hneigh !> # H atoms bound to atom
    inv = 0
    inv = inv+nneigh*10**5
    inv = inv+nbonds*10**4
    inv = inv+ati*10**2
    inv = inv+chrg*10**1
    inv = inv+hneigh
  end function get_invariant0

!=========================================================================================!

  subroutine update_ranks(self)
    !*************************************************************
    !* Refine the ranks: the atoms of each rank are split by their
    !* invariant, and the new ranks are numbered densely in order
    !* of (old rank, invariant). A sort keeps this O(n log n).
    !* Also assigns the prime of each new rank.
    !*************************************************************
    implicit none
    class(canonical_sorter) :: self
    integer,allocatable :: idx(:),work(:)
    integer :: i,n,r
    n = self%hatms
    if (n < 1) return
    allocate (idx(n),work(n))
    do i = 1,n
      idx(i) = i
    end do
    call sort_atoms(self,idx,work,1,n)
    r = 1
    self%newrank(idx(1)) = r
    do i = 2,n
      if (compare_atoms(self,idx(i-1),idx(i)) /= 0) r = r+1
      self%newrank(idx(i)) = r
      if (debug) write (stdout,*) 'new rank',r,'assigned to atom',idx(i)
    end do
!>--- assign primes
    do i = 1,n
      self%prime(i) = self%primetab(self%newrank(i))
    end do
  end subroutine update_ranks

!=========================================================================================!

  recursive subroutine sort_atoms(self,idx,work,lo,hi)
    !*********************************************
    !* Merge sort of idx(lo:hi) by compare_atoms.
    !*********************************************
    implicit none
    class(canonical_sorter),intent(in) :: self
    integer,intent(inout) :: idx(:),work(:)
    integer,intent(in) :: lo,hi
    integer :: mid,i,j,k
    if (hi <= lo) return
    mid = (lo+hi)/2
    call sort_atoms(self,idx,work,lo,mid)
    call sort_atoms(self,idx,work,mid+1,hi)
    i = lo
    j = mid+1
    k = lo
    do while (i <= mid.and.j <= hi)
      if (compare_atoms(self,idx(j),idx(i)) < 0) then
        work(k) = idx(j)
        j = j+1
      else
        work(k) = idx(i)
        i = i+1
      end if
      k = k+1
    end do
    work(k:k+mid-i) = idx(i:mid)
    k = k+mid-i+1
    work(k:hi) = idx(j:hi)
    idx(lo:hi) = work(lo:hi)
  end subroutine sort_atoms

  integer function compare_atoms(self,a,b) result(c)
    !****************************************************************
    !* Order of atoms a and b by (rank, invariant): -1, 0 or 1.
    !* An overflowed prime product is larger than any int64 one; two
    !* overflowed ones are ordered by their sorted neighbour primes,
    !* which is equal exactly when the products are.
    !****************************************************************
    implicit none
    class(canonical_sorter),intent(in) :: self
    integer,intent(in) :: a,b
    integer :: pa(self%maxnei),pb(self%maxnei),na,nb,l
    c = 0
    if (self%rank(a) /= self%rank(b)) then
      c = merge(-1,1,self%rank(a) < self%rank(b))
    else if (self%ovf(a).neqv.self%ovf(b)) then
      c = merge(1,-1,self%ovf(a))
    else if (.not.self%ovf(a)) then
      if (self%invariants(a) /= self%invariants(b)) then
        c = merge(-1,1,self%invariants(a) < self%invariants(b))
      end if
    else
      na = self%nhnei(a)
      nb = self%nhnei(b)
      if (na /= nb) then
        c = merge(-1,1,na < nb)
        return
      end if
      pa(1:na) = self%prime(self%hneigh(1:na,a))
      pb(1:nb) = self%prime(self%hneigh(1:nb,b))
      call isort(pa(1:na))
      call isort(pb(1:nb))
      do l = 1,na
        if (pa(l) /= pb(l)) then
          c = merge(-1,1,pa(l) < pb(l))
          return
        end if
      end do
    end if
  end function compare_atoms


!=========================================================================================!

  subroutine update_invariants(self)
    !****************************************************************
    !* New invariant of each atom: the product of the primes of its
    !* heavy neighbours' ranks. A product that would overflow int64
    !* is flagged in ovf instead and compared via its factors.
    !****************************************************************
    implicit none
    class(canonical_sorter) :: self
    integer :: i,l
    integer(int64) :: invprod,p
    do i = 1,self%hatms
      invprod = 1
      self%ovf(i) = .false.
      do l = 1,self%nhnei(i)
        p = int(self%prime(self%hneigh(l,i)),int64)
        if (invprod > huge(invprod)/p) then
          self%ovf(i) = .true.
          invprod = 0
          exit
        end if
        invprod = invprod*p
      end do
      self%invariants(i) = invprod
    end do
    if (debug) write (stdout,*) 'update invariants'
  end subroutine update_invariants

!========================================================================================!

  subroutine iterate(self,mol)
    !*************************************************
    !* Refine ranks and invariants until the ranking
    !* no longer changes.
    !*************************************************
    implicit none
    class(canonical_sorter) :: self
    type(coord) :: mol
    do
      call self%update_ranks()
      if (debug) then
        call debugprint(self,mol)
      end if
      if (all(self%rank(:) .eq. self%newrank(:))) exit  !> termination condition
      self%rank(:) = self%newrank(:)
      call self%update_invariants()
      if (debug) then
        call debugprint(self,mol)
      end if
    end do
  end subroutine iterate

!========================================================================================!

  subroutine renumber_heavy_ranks(self,mol)
    !******************************************************************
    !* With H atoms carried as graph nodes (heavy=.false.), number the
    !* heavy-atom ranks densely from 1 while keeping their order, and
    !* set the H entries to 0 for add_h_ranks. The H nodes normally
    !* form rank 1 alone, which makes this a shift by one.
    !******************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    integer,allocatable :: newnum(:)
    integer :: i,r
    if (self%hatms < 1) return
    allocate (newnum(maxval(self%rank)),source=0)
    do i = 1,self%hatms
      if (mol%at(self%hmap(i)) .ne. 1) newnum(self%rank(i)) = 1
    end do
    r = 0
    do i = 1,size(newnum)
      if (newnum(i) > 0) then
        r = r+1
        newnum(i) = r
      end if
    end do
    do i = 1,self%hatms
      if (mol%at(self%hmap(i)) .eq. 1) then
        self%rank(i) = 0
      else
        self%rank(i) = newnum(self%rank(i))
      end if
    end do
  end subroutine renumber_heavy_ranks

!========================================================================================!

  subroutine analyze_stereo(self,mol)
    implicit none
    class(canonical_sorter) :: self
    type(coord),intent(in) :: mol
    integer :: i,ii,zero,nei,j,jj,maxrank
    integer :: k,rs
    integer,allocatable :: neiranks(:,:)
    real(wp) :: coords(3,4)
    logical,allocatable :: isstereo(:)
    allocate (isstereo(mol%nat),source=.false.)
    allocate (neiranks(4,mol%nat),source=0)
    maxrank = maxval(self%rank(:))
    do i = 1,self%hatms
      ii = self%hmap(i)
      zero = count(self%neigh(:,ii) == 0)
      nei = self%maxnei-zero
!>--- consider only atoms with 4 unique (in terms of CANGEN ranks) neighbours as stereocenter
      if (nei == 4) then
        do j = 1,4
          jj = self%neigh(j,ii)
          if (mol%at(jj) == 1) then !> one hydrogen allowed
            neiranks(j,ii) = maxrank+1
          else
            neiranks(j,ii) = self%rank(jj)
          end if
        end do
        !isstereo(ii) = unique_neighbours(nei,self%neigh(1:nei,ii))
        isstereo(ii) = unique_neighbours(4,neiranks(:,ii))
      end if
    end do
!>--- do some actual geometry checks to determine "R" and "S"
    do i = 1,mol%nat
      if (isstereo(i)) then
        !>--- transfer to dummy coords
        do j = 1,4
          k = minloc(neiranks(:,i),1)
          jj = self%neigh(k,i)
          neiranks(k,i) = huge(k)
          coords(1:3,j) = mol%xyz(1:3,jj)-mol%xyz(1:3,i) !> shift i to center of coordinates
        end do
        RS = determineRS(coords)
        ii = self%nmap(i)
        self%invariants0(ii) = (self%invariants0(ii)+abs(RS)*10**6)*RS
      end if
    end do

    if (debug) then
      do ii = 1,self%hatms
        self%invariants(ii) = self%invariants0(ii)
      end do
      call debugprint(self,mol)
    end if
    deallocate (neiranks,isstereo)
  end subroutine analyze_stereo

!===========================================================================================!

  function has_stereo(self,mol) result(yesno)
!*****************************************************************
!* Reports whether the molecule has a (CANGEN-rank) stereocenter.
!*
!* Only the molecular graph (neighbour list, heavy-atom map) and a
!* rank array are required. When the graph has not been initialized
!* yet - e.g. because the ranks were supplied externally rather than
!* computed via init - it is self-filled here from the connectivity,
!* so callers can run the stereo check on provided ranks without a
!* full canonical_sorter%init.
!*****************************************************************
    implicit none
    logical :: yesno
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    integer :: i,ii,zero,nei,j,jj,maxrank
    integer,allocatable :: neiranks(:,:)
    logical,allocatable :: isstereo(:)
    integer,allocatable :: Amat(:,:)

    !> self-fill the molecular graph if only ranks (but no graph) are present.
    !> include_H=.true. keeps hmap/hatms consistent with full-molecule ranks.
    if (.not.allocated(self%neigh)) then
      call self%init_graph(mol,.true.,Amat)
      if (allocated(Amat)) deallocate (Amat)
    end if

    allocate (isstereo(mol%nat),source=.false.)
    allocate (neiranks(4,mol%nat),source=0)
    maxrank = maxval(self%rank(:))
    do i = 1,self%hatms
      ii = self%hmap(i)
      zero = count(self%neigh(:,ii) == 0)
      nei = self%maxnei-zero
!>--- consider only atoms with 4 unique (in terms of ranks) neighbours as stereocenter
      if (nei == 4) then
        do j = 1,4
          jj = self%neigh(j,ii)
          if (mol%at(jj) == 1) then !> one hydrogen allowed
            neiranks(j,ii) = maxrank+1
          else
            neiranks(j,ii) = self%rank(jj)
          end if
        end do
        isstereo(ii) = unique_neighbours(4,neiranks(:,ii))
      end if
    end do
    yesno = any(isstereo(:))
    deallocate (neiranks,isstereo)
  end function has_stereo

!========================================================================================!

  function compare_canonical_sorter(self,other) result(yesno)
!*****************************************************
!* compare two canonical_sorter objects to determine
!* if both correspond to the same molecule
!*****************************************************
    implicit none
    class(canonical_sorter) :: self
    type(canonical_sorter) :: other
    logical :: yesno
    integer,allocatable :: sorted_invariants(:,:)
    integer :: maxrank_a,maxrank_b,hatms
    integer :: i,j,jj,k,kk

    yesno = .false.
!>--- the obvious cases first
    if (self%nat .ne. other%nat) return
    if (self%hatms .ne. other%hatms) return
    maxrank_a = maxval(self%rank(:),1)
    maxrank_b = maxval(other%rank(:),1)
    if (maxrank_a .ne. maxrank_b) return

!>--- if the easy checks passed, compare the actual invariants!
    hatms = self%hatms
    allocate (sorted_invariants(hatms,2),source=0)
    jj = 0
    kk = 0
    do i = 1,maxrank_a  !> maxrank_a == maxrank_b, see above
      do j = 1,hatms
        if (self%rank(j) == i) then
          jj = jj+1
          sorted_invariants(jj,1) = self%invariants0(j)
        end if
      end do
      do k = 1,hatms
        if (other%rank(k) == i) then
          kk = kk+1
          sorted_invariants(kk,2) = other%invariants0(k)
        end if
      end do
    end do

    if (debug) then
      do i = 1,hatms
        write (stdout,'(i10,i10,i10)') i,sorted_invariants(i,:)
      end do
    end if

!>--- assign identity
    yesno = all(sorted_invariants(:,1) .eq. sorted_invariants(:,2))

    deallocate (sorted_invariants)
    return
  end function compare_canonical_sorter

!========================================================================================!

  subroutine add_h_ranks(self,mol)
!******************************************************************
!* Mapps ranks of the heavy atoms back to the full molecule order
!* And continues ranks for H atoms, based on neighbor list
!******************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    integer,allocatable :: rankmap(:)
    integer :: i,ii,jj,maxrank,rr
!>--- self%rank must already have the correct dimension!
    if (size(self%rank,1) .ne. mol%nat) then
      stop 'wrong dimension for adding H to canonical ranks!'
    end if

!>--- otherwise, analyze and resize
    maxrank = maxval(self%rank(:),1)

!>--- cycle through atoms, assign ranks depending on neighbour list
    allocate (rankmap(maxrank),source=0)
    rr = 0
    do i = 1,self%hatms
      if (mol%at(i) .ne. 1) cycle
      ii = self%neigh(1,i)
      if (ii < 1) cycle  !> Edge-case: "unbound" hydrogens (H+, H2, etc.), skip those here
      if (mol%at(ii) .eq. 1) cycle
      jj = self%rank(ii)
      rankmap(jj) = 1
    end do
    do i = 1,maxrank
      if (rankmap(i) .eq. 1) then
        rr = rr+1
        rankmap(i) = maxrank+rr
      end if
    end do
    !> one more rank above all others for the unbound H
    maxrank = max(maxrank,maxval(rankmap(:),1))+1
    do i = 1,self%hatms
      if (mol%at(i) .ne. 1) cycle
      ii = self%neigh(1,i)
      unbound: block
        if (ii >= 1) then
          if (mol%at(ii) .ne. 1) exit unbound
        end if
        !> "unbound" H (no neighbour, or only H neighbours) share one extra rank
        self%rank(i) = maxrank
        cycle
      end block unbound
      jj = self%rank(ii)
      self%rank(i) = rankmap(jj)
    end do
    deallocate (rankmap)
  end subroutine add_h_ranks

!========================================================================================!

  function same_graph(self,ref,mol) result(yesno)
    !*****************************************************************
    !* .true. if mol has the atom types of ref and the connectivity
    !* stored in self, which must have been set up from ref by init
    !* (or init_graph) without a wbo matrix. The ranks depend on
    !* nothing else, so they can then be copied from ref to mol
    !* instead of being recomputed.
    !*****************************************************************
    implicit none
    class(canonical_sorter),intent(in) :: self
    type(coord),intent(in) :: ref,mol
    logical :: yesno
    integer,allocatable :: A(:,:)
    integer :: i,j,l
    yesno = .false.
    if (.not.allocated(self%neigh)) return
    if (ref%nat .ne. self%nat.or.mol%nat .ne. self%nat) return
    if (any(mol%at(:) .ne. ref%at(:))) return
    call mol%cn_to_adjacency(A)
    do i = 1,mol%nat
      l = 0
      do j = 1,mol%nat
        if (A(j,i) .eq. 0) cycle
        l = l+1
        if (l > self%maxnei) return
        if (self%neigh(l,i) .ne. j) return
      end do
      if (l < self%maxnei) then
        if (self%neigh(l+1,i) .ne. 0) return
      end if
    end do
    yesno = .true.
  end function same_graph

!========================================================================================!

  subroutine split_pair_ranks(self,mol)
!*********************************************************************
!* NMR magnetic-equivalency hack (only for explicit rank requests):  *
!* any canonical rank shared by EXACTLY two non-hydrogen atoms is    *
!* split into two distinct ranks. The second atom (in atom-index     *
!* order) receives a fresh rank above the current maximum, the first *
!* keeps the original rank.                                          *
!*                                                                   *
!* Must run AFTER the heavy-atom rank refinement but BEFORE the      *
!* hydrogen ranks are added, so the new distinction propagates to    *
!* attached H atoms through add_h_ranks. Reachable only via the      *
!* 'apsp+nmr' invtype, i.e. the get_canonical path; the iRMSD        *
!* matching routines never request it.                               *
!*                                                                   *
!* self : canonical_sorter (rank already heavy-atom refined)         *
!* mol  : the molecule (used only for atom types via hmap)           *
!*********************************************************************
    implicit none
    class(canonical_sorter),intent(inout) :: self
    type(coord),intent(in) :: mol
    integer :: i,r,maxrank,nmembers,second,newrank

    maxrank = maxval(self%rank(:),1)
    newrank = maxrank
    do r = 1,maxrank
! ── count non-H atoms carrying rank r, remembering the second one ───────────
      nmembers = 0
      second = 0
      do i = 1,self%hatms
        if (mol%at(self%hmap(i)) .eq. 1) cycle  !> never split H ranks
        if (self%rank(i) .ne. r) cycle
        nmembers = nmembers+1
        if (nmembers .eq. 2) second = i
      end do
! ── split only exact pairs into two unique ranks ────────────────────────────
      if (nmembers .eq. 2) then
        newrank = newrank+1
        self%rank(second) = newrank
      end if
    end do
  end subroutine split_pair_ranks

!========================================================================================!
!========================================================================================!

  subroutine debugprint(can,mol)
    implicit none
    type(canonical_sorter) :: can
    type(coord)  :: mol
    integer :: i,ii,ati
    write (stdout,'(a10,a5,a15,a10,a10)') 'heavy-atom','type','invariant','rank','prime'
    do i = 1,can%hatms
      ii = can%hmap(i)
      ati = mol%at(ii)
      write (stdout,'(i10,a5,i15,i10,i10)') i,i2e(ati,'nc'),can%invariants(i),can%rank(i),can%prime(i)
    end do
  end subroutine debugprint

  subroutine rankprint(can,mol)
    implicit none
    class(canonical_sorter) :: can
    type(coord)  :: mol
    integer :: i,ii,ati
    write (stdout,'(a10,a10,a12,a10,2x,a)') 'heavy-atom','type','invariant0','rank','neighbours'
    do i = 1,can%hatms
      ii = can%hmap(i)
      ati = mol%at(ii)
      write (stdout,'(i10,a10,i12,i10,2x,a)') i,i2e(ati,'nc'),can%invariants0(i),can%rank(i),print_neighbours(mol,can%neigh(:,ii))
    end do
  end subroutine rankprint

  function print_neighbours(mol,neigh) result(btmp)
    implicit none
    type(coord) :: mol
    integer,intent(in) :: neigh(:)
    character(len=:),allocatable :: btmp
    character(len=20) :: atmp
    integer :: i
    btmp = ''
    if (neigh(1) == 0) then
      btmp = ' ---'
      return
    end if
    do i = 1,size(neigh,1)
      if (neigh(i) .ne. 0) then
        write (atmp,'(a,a,a,i0,a)') ' ',trim(i2e(mol%at(neigh(i)))),'(',neigh(i),')'
        btmp = trim(btmp)//trim(atmp)
      end if
    end do
  end function print_neighbours

!========================================================================================!

  function unique_neighbours(n,arr) result(yesno)
    implicit none
    integer,intent(in) :: n
    integer,intent(in) :: arr(n)
    logical :: yesno
    integer :: i,j
    yesno = .true.
    do i = 1,n-1
      do j = i+1,n
        if (arr(i) == arr(j)) then
          yesno = .false.
          exit
        end if
      end do
      if (.not.yesno) exit
    end do
  end function unique_neighbours

!========================================================================================!

  function determineRS(coords) result(RS)
    implicit none
    integer :: RS !> 1 for "R", -1 for "S"
    real(wp),intent(inout) :: coords(3,4)
    real(wp) :: theta
    real(wp) :: vec(3),uec(3)
    integer :: k,l

    k = 4
    !> rotate the highest prio atom onto z axis (0,0,1)
    !> first into into xz plane, therefore get the -x unit vector first
    !> and the projection of the neighbour in the xy-plane
    vec = (/-1.0d0,0.0d0,0.0d0/)
    uec(1) = coords(1,1)
    uec(2) = coords(2,1)
    uec(3) = 0.0_wp
    theta = tangle(vec,uec)
    !> then rotate around z-axis
    vec = (/0.0d0,0.0d0,1.0d0/)
    if (uec(2) .lt. 0.0_wp) theta = -theta
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !> afterwards angle to the z-axis
    vec = (/0.0d0,0.0d0,1.0d0/)
    uec = coords(:,1)
    theta = tangle(vec,uec)
    !> then take this angle and rotate around y
    vec = (/0.0d0,1.0d0,0.0d0/)
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !> as a last step tak the lowest-prio neighbour and rotate it into the xy-plane
    vec = (/-1.0d0,0.0d0,0.0d0/)
    uec(1) = coords(1,k)
    uec(2) = coords(2,k)
    uec(3) = 0.0_wp
    theta = tangle(vec,uec)
    vec = (/0.0,0.0,1.0/)
    if (uec(2) .lt. 0.0_wp) theta = -theta
    do l = 1,k
      call rodrot(coords(:,l),vec,theta)
    end do

    !determine orientation
    if (coords(2,2) .gt. coords(2,3)) then
      RS = 1
    else
      RS = -1
    end if

  end function determineRS

!========================================================================================!
!========================================================================================!
end module canonical_mod
