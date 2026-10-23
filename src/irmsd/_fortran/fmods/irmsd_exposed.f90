module irmsd_exposed
  use,intrinsic :: iso_c_binding
  use,intrinsic :: iso_fortran_env,only:wp => real64
  use strucrd
  use crest_parameters
  use axis_module
  use irmsd_module
  use canonical_mod
  use irmsd_refine
  implicit none
contains

  subroutine get_irmsd_fortran(natoms1,types1_ptr,coords1_ptr,ranks1_ptr, &
                               natoms2,types2_ptr,coords2_ptr,ranks2_ptr, &
                               iinversion_c,rmsd_c,types_out1_ptr,coords_out1_ptr, &
                               types_out2_ptr,coords_out2_ptr) &
    bind(C,name="get_irmsd_fortran")
    !*********************************************************************
    !* Compute the iRMSD between two structures of equal atom count.     *
    !*                                                                   *
    !* If the caller supplies per-atom canonical ranks (ranks1/ranks2,   *
    !* all entries non-zero) and they pass checkranks(), those ranks are *
    !* used directly and the (expensive) canonical_sorter run is skipped *
    !* where possible. A zero in either array (the "not provided"        *
    !* sentinel) or a checkranks() rejection falls back to the regular   *
    !* canonical initialization.                                         *
    !*                                                                   *
    !* natoms1/types1/coords1/ranks1 : reference structure (ranks input) *
    !* natoms2/types2/coords2/ranks2 : mobile structure    (ranks input) *
    !* iinversion_c : 0=auto, 1=force inversion on, 2=force off          *
    !* rmsd_c       : (out) iRMSD value in Angstrom                      *
    !* *_out*       : (out) aligned types/coordinates of both structures *
    !*********************************************************************
    use,intrinsic :: iso_c_binding
    implicit none
    !> IN-/OUTPUTS
    integer(c_int),value :: natoms1,natoms2,iinversion_c
    type(c_ptr),value :: types1_ptr,coords1_ptr,ranks1_ptr
    type(c_ptr),value :: types2_ptr,coords2_ptr,ranks2_ptr
    type(c_ptr),value :: types_out1_ptr
    type(c_ptr),value :: coords_out1_ptr
    type(c_ptr),value :: types_out2_ptr
    type(c_ptr),value :: coords_out2_ptr
    real(c_double),intent(out) :: rmsd_c

    integer :: iinversion
    integer(c_int),pointer :: types_out1_c(:)
    real(c_double),pointer :: coords_out1_c(:)
    integer(c_int),pointer :: types_out2_c(:)
    real(c_double),pointer :: coords_out2_c(:)
    type(coord) :: mol,ref
    real(wp) :: rmsdval,tmpd(3),tmpdist
    integer :: i
    type(rmsd_cache) :: rcache
    logical :: mirror
    integer :: io
    real(wp) :: rmsdref

    logical,parameter :: debug = .false.

    !> Externally supplied canonical ranks ride along via the ranks pointers;
    !> C_to_mol stores them in ref%id / mol%id only if a complete (all
    !> non-zero) set is given, so setup_irmsd_ranks can decide per structure.
    call ref%C_to_mol(natoms1,types1_ptr,coords1_ptr,.true.,ranks1_ptr)
    call mol%C_to_mol(natoms2,types2_ptr,coords2_ptr,.true.,ranks2_ptr)

    if (natoms1 /= natoms2) then
      error stop 'both molecules need to have the same number of atoms'
    end if

    iinversion = iinversion_c

    !> move ref to CMA and align rotational axes
    call axis(ref%nat,ref%at,ref%xyz)

    !> allocate memory
    call rcache%allocate(ref%nat)

    !> determine the per-atom ranks (provided ids where available, otherwise
    !> recomputed) and the false-enantiomer flag, then apply inversion override
    call setup_irmsd_ranks(ref,mol,rcache%rank,rcache%stereocheck)
    select case (iinversion)
    case (0)  !> whatever the stereo check says
      mirror = .true.
    case (1)  !> force on
      mirror = .true.
      rcache%stereocheck = .true.
    case (2)  !> force off
      mirror = .false.
      rcache%stereocheck = .false.
    end select

    if (debug) write (stdout,*) 'allow inversion?:            ',mirror

    call min_rmsd(ref,mol,rcache=rcache,rmsdout=rmsdval,align=.true.,io=io)

    !> bond-preserving assignment: centre walk vs. min_rmsd, terminal atoms refined
    if (io == 0) then
      call refine_assignment(ref,mol,rcache%rank(:,1),rcache%stereocheck,rmsdref,io)
      if (io == 0) rmsdval = rmsdref
    end if

    if (debug) then
      do i = 1,mol%nat
        tmpd(:) = (mol%xyz(:,i)-ref%xyz(:,i))**2
        tmpdist = sqrt(sum(tmpd(:)))*autoaa
        if (tmpdist > 0.01_wp) then
          write (stdout,*) i,mol%at(i),tmpdist
        end if
      end do
    end if

    rmsdval = rmsdval*autoaa
    rmsd_c = real(rmsdval,c_double)

    call c_f_pointer(types_out1_ptr,types_out1_c, [natoms1])
    call c_f_pointer(coords_out1_ptr,coords_out1_c, [3*natoms1])

    call c_f_pointer(types_out2_ptr,types_out2_c, [natoms2])
    call c_f_pointer(coords_out2_ptr,coords_out2_c, [3*natoms2])

    call ref%mol_to_C(types_out1_c,coords_out1_c,.true.)
    call mol%mol_to_C(types_out2_c,coords_out2_c,.true.)

  end subroutine get_irmsd_fortran

!========================================================================================!

  subroutine bond_mismatch(ref,mol,nmis,atmis)
    !*******************************************************************
    !* Check whether the correspondence ref atom i <-> mol atom i
    !* preserves the bonding pattern. Bonds are detected separately for
    !* each structure from its own geometry, with the same covalent CN
    !* criterion the canonical ranks are built from.
    !*
    !* ref, mol : structures with equal atom count, in matched order
    !*            (e.g. the output of min_rmsd)
    !* nmis     : (out) number of atom pairs bonded in exactly one of
    !*            the two structures; 0 means the mapping is bond-preserving
    !* atmis    : (out) per-atom number of such pairs
    !*******************************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(out) :: nmis
    integer,intent(out) :: atmis(:)
    integer,allocatable :: aref(:,:),amol(:,:)
    integer :: i,j

    call ref%cn_to_adjacency(aref)
    call mol%cn_to_adjacency(amol)

    nmis = 0
    atmis(:) = 0
    do i = 1,ref%nat
      do j = 1,i-1
        if (aref(j,i) /= amol(j,i)) then
          nmis = nmis+1
          atmis(i) = atmis(i)+1
          atmis(j) = atmis(j)+1
        end if
      end do
    end do
  end subroutine bond_mismatch

  subroutine get_bond_mismatch_fortran(natoms,types1_ptr,coords1_ptr, &
                                       types2_ptr,coords2_ptr,nmis_c,atmis_ptr) &
    bind(C,name="get_bond_mismatch_fortran")
    !*******************************************************************
    !* C interface to bond_mismatch. Both structures must be given in
    !* matched atom order, e.g. the aligned output of get_irmsd_fortran.
    !*
    !* natoms        : number of atoms of both structures
    !* types*/coords*: atomic numbers and coordinates in Angstrom
    !* nmis_c        : (out) number of mismatched atom pairs
    !* atmis_ptr     : (out) per-atom number of mismatched pairs
    !*******************************************************************
    use,intrinsic :: iso_c_binding
    implicit none
    integer(c_int),value :: natoms
    type(c_ptr),value :: types1_ptr,coords1_ptr
    type(c_ptr),value :: types2_ptr,coords2_ptr
    integer(c_int),intent(out) :: nmis_c
    type(c_ptr),value :: atmis_ptr
    integer(c_int),pointer :: atmis_c(:)
    type(coord) :: ref,mol
    integer,allocatable :: atmis(:)
    integer :: nmis

    call ref%C_to_mol(natoms,types1_ptr,coords1_ptr,.true.)
    call mol%C_to_mol(natoms,types2_ptr,coords2_ptr,.true.)
    allocate (atmis(natoms))

    call bond_mismatch(ref,mol,nmis,atmis)

    nmis_c = int(nmis,c_int)
    call c_f_pointer(atmis_ptr,atmis_c, [natoms])
    atmis_c(:) = int(atmis(:),c_int)
  end subroutine get_bond_mismatch_fortran
end module irmsd_exposed
