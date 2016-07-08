module pdb_module
  implicit none
  integer :: iatom

  type pdb_class
    integer                       :: n_atoms, n_residues, n_chains
    integer         , allocatable :: AtomNum(:)
    character(len=4), allocatable :: AtomName(:)
    character(len=3), allocatable :: ResName(:)
    character(len=1), allocatable :: ChainId(:)
    integer         , allocatable :: ResNum(:)
    double precision, allocatable :: x(:),y(:),z(:)
    character(len=126) :: blnk=""
  contains
    procedure :: initialize
    procedure :: set_n_atoms
    procedure :: get_structure
    procedure :: output_ref
    procedure :: get_coords
    procedure, private :: count_atom_type
    procedure :: select_atoms
  end type
  
contains
  
  subroutine set_n_atoms(self, n_atoms) 
    class(pdb_class) :: self
    integer          :: n_atoms
    self%n_atoms = n_atoms
  end subroutine

  function get_n_atoms(self) result(n_atoms)
    class(pdb_class) :: self
    integer          :: n_atoms 

  end function

  subroutine initialize(self)
    class(pdb_class) :: self
    print*, "no. of atoms: ", self%n_atoms
    allocate( self%AtomNum(self%n_atoms), &
              self%AtomName(self%n_atoms),&
              self%ResName(self%n_atoms), &
              self%ChainId(self%n_atoms), &
              self%ResNum(self%n_atoms),  &
              self%x(self%n_atoms),self%y(self%n_atoms),self%z(self%n_atoms) )
    self%AtomNum(:) = 0
    self%ResNum(:)  = 0
    self%x(:) = 0.0;  self%y(:) = 0.0;  self%z(:) = 0.0
  end subroutine

  subroutine get_structure(self, filename)
    class(pdb_class) :: self
    integer          :: unit_no = 100
    character(len=*), intent(in) :: filename
    character(len=6) :: ATOM
    integer :: iatom

    write(*,*) "Unit No. is ", unit_no

    call self%initialize()

    open(unit_no, file = filename, status="old")
    do
      read(unit_no,*) ATOM

      select case(ATOM)
        case("ATOM")
          backspace(unit_no)
          do iatom = 1, self%n_atoms
            read(unit_no,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
                 ATOM       , &
                 self%AtomNum(iatom) , &
                 self%AtomName(iatom), &
                 self%ResName(iatom) , &
                 self%ChainId(iatom) , &
                 self%ResNum(iatom)  , &
                 self%x(iatom),        &
                 self%y(iatom),        &
                 self%z(iatom)
            !print*,self%AtomNum(iatom), self%AtomName(iatom), self%ResName(iatom), self%ChainId(iatom),&
            !       self%x(iatom), self%y(iatom), self%z(iatom)
          enddo
          exit
       case("REMARK") 
           print*, "THis is a REMARK line."
       case("HETATM")
           print*, "this is a HETATM line."
       case("TER")
           stop "A TER line is included, remove it."
       case default
           print*, "A weird line was etected:", ATOM
           stop "Remove it."
       end select
    enddo
    close(unit_no) 
  end subroutine

  subroutine output_ref(self)
      class(pdb_class) :: self
      integer :: unit = 111
      open(unit, file ="Reference.pdb", status="replace")

      do iatom = 1, size(self%AtomNum) 
        write(unit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
           "ATOM  ", &
           self%AtomNum(iatom), &
           self%AtomName(iatom),&
           self%ResName(iatom),&
           self%ChainId(iatom), &
           self%ResNum(iatom), &
           self%x(iatom),        &
           self%y(iatom),        &
           self%z(iatom),        & 
           self%blnk
      enddo
      write(*,"('>The reference PDB was outputed.')")
      close(unit)
  end subroutine

  function count_atom_type(self, atom_type, res_range) result(ncount)
    class(pdb_class)              :: self
    character(len=*), intent(in)  :: atom_type
    integer                       :: ncount, iatom_typ
    integer,intent(in)            :: res_range(2) 

    iatom_typ = 0
    do iatom = 1, size(self%x, 1)
      select case(atom_type)
      case("CA")
        if ( self%AtomName(iatom) == " CA "     .and. &
             self%ResNum(iatom) >= res_range(1) .and. &
             self%ResNum(iatom) <= res_range(2)       &
           ) iatom_typ = iatom_typ + 1
      case("BB")
        if ( (self%AtomName(iatom) == " CA " .or. &
              self%AtomName(iatom) == " N "  .or. &
              self%AtomName(iatom) == " C " )     &
              .and. & 
             (self%ResName(iatom)  /= "ACE" .and. &
              self%ResName(iatom)  /= "NME" .and. &
              self%ResName(iatom)  /= "NH2"  )    &
              .and. &
             (self%ResNum(iatom) >= res_range(1) .and. &
              self%ResNum(iatom) <= res_range(2))      &
              ) then 
          iatom_typ = iatom_typ + 1
        endif
      end select
    enddo
    ncount = iatom_typ
    !if (index(self%AtomName(iatom), atom_type) /= 0) then 
  end function


  !@@@@@make this!
  subroutine select_atoms(self, init_resi, last_resi)
    class(pdb_class)    :: self
    integer, intent(in) :: init_resi, last_resi  ! e.g.  1 15
    print*,"resi init-last: ", init_resi, last_resi
    if (init_resi >= last_resi) stop "resi(1) is more than resi(2)."
    do iatom = 1, size(self%ResNum,1)
      if (self%ResNum(iatom) >= init_resi .and. self%ResNum(iatom) <= last_resi ) then
        print*, "res no: ", self%AtomName(iatom), self%ResNum(iatom), self%ResName(iatom)
      endif
    enddo
  end subroutine

  function get_coords(self, atom_type, res_range) result(coords_atom_type)
    class(pdb_class)              :: self
    character(len=*) , intent(in) :: atom_type
    integer, optional, intent(in) :: res_range(2)
    real(8), allocatable          :: coords_atom_type(:,:)
    integer                       :: iatom_typ
    integer                       :: n_specified_atoms

    n_specified_atoms = self%count_atom_type(atom_type, res_range)
    allocate( coords_atom_type(3, n_specified_atoms))
    write(*,"('    SPECIFIED ATOM TYPE  : ', a5)") atom_type
    write(*,"('    NO OF COORDINATES    : ', i2)") size(coords_atom_type, 1)
    write(*,"('    NO OF SPECIFIED ATOMS: ', i7)") size(coords_atom_type, 2)
    coords_atom_type(:,:) = 0.0d0

    write(*,"('    # Fitting residue range', 2i7)") res_range
    write(*,"('    # Fitting residue list')")
    
    iatom_typ = 1
    do iatom = 1, size(self%x, 1)

      select case(atom_type)

      case("CA")

        if (self%AtomName(iatom) == " CA "     .and. & 
            self%ResNum(iatom) >= res_range(1) .and. &
            self%ResNum(iatom) <= res_range(2)       &
           ) then 
          coords_atom_type(1,iatom_typ) = self%x(iatom) 
          coords_atom_type(2,iatom_typ) = self%y(iatom)
          coords_atom_type(3,iatom_typ) = self%z(iatom)
          !write(*,"(7x,a3,1x,i7,1x,a3,1x,3f10.3)") self%AtomName(iatom),&
          !                                         self%ResNum(iatom)  ,&
          !                                         self%ResName(iatom) ,&
          !                                         coords_atom_type(:,iatom_typ)
          iatom_typ = iatom_typ + 1
        endif


      case("BB")

        if ( (self%AtomName(iatom) == " CA "     .or.  &
              self%AtomName(iatom) == " N "      .or.  &
              self%AtomName(iatom) == " C "  )   .and. &
             (self%ResName(iatom)  /= "ACE"      .and. &
              self%ResName(iatom)  /= "NME"      .and. &
              self%ResName(iatom)  /= "NH2"  )   .and. & 
             (self%ResNum(iatom) >= res_range(1) .and. &
              self%ResNum(iatom) <= res_range(2))      &
           ) then 
          !write(*,"('No', i10, 2a4)") iatom, self%ResName(iatom), self%AtomName(iatom)
          coords_atom_type(1,iatom_typ) = self%x(iatom) 
          coords_atom_type(2,iatom_typ) = self%y(iatom)
          coords_atom_type(3,iatom_typ) = self%z(iatom)
          iatom_typ = iatom_typ + 1
        endif

      case default
        stop 'atom_type is wrong.'

      end select
    enddo

    print*, "no of specified atoms: ", iatom_typ -1

  end function



end module
