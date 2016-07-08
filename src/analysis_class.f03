module analysis_module
  use variables
  use input_module
  use pdb_module
  use calculation
  use alignment_module
  implicit none
  character(len=6)             :: ATOM
  integer :: iii

  !???Inheritance is Ok????
  type, extends(pdb_class) :: analysis_class
  contains
    procedure :: analysing_in_loop
    procedure :: read_a_frame
  end type
contains
  
  subroutine read_a_frame(self, n_atoms, unit_no)
    class(analysis_class) :: self
    integer, intent(in)    :: n_atoms   
    integer :: iatom, unit_no
    do iatom = 1, n_atoms
      read(unit_no,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
                     ATOM                , &
                     self%AtomNum(iatom) , &
                     self%AtomName(iatom), &
                     self%ResName(iatom) , &
                     self%ChainId(iatom) , &
                     self%ResNum(iatom)  , &
                     self%x(iatom)       , &
                     self%y(iatom)       , &
                     self%z(iatom)
      !write(*,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") & 
      !     ATOM, self%AtomNum(iatom), self%AtomName(iatom), self%ResName(iatom), &
      !     self%ChainId(iatom), self%ResNum(iatom), &
      !     self%x(iatom), self%y(iatom), self%z(iatom)
    enddo
  end subroutine

  function convert_xyz_r(x,y,z) result(coords)
    integer i
    real(8) :: x(:),y(:),z(:)
    real(8), allocatable :: coords(:,:)
    allocate(coords(3,size(x,1)))
    do i = 1, size(x,1)
      coords(1,i) = x(i)
      coords(2,i) = y(i)
      coords(3,i) = z(i)
    enddo
   end function

  subroutine analysing_in_loop(self, filename, pdb, inputs)
    class(analysis_class)         :: self
    character(len=*) , intent(in) :: filename
    type(pdb_class)  , intent(in) :: pdb
    type(input_class), intent(in) :: inputs
    type(coordinates)             :: ref, specified, snapshot, specified_in_snapshot, rotated
    integer                       :: unit_no = 1122, iatom,imodel
    integer                       :: temp_var_location_eigen
    real(8)                       :: S(4,4), eigenval(4), eigenvec(4,4)
    real(8)                       :: Quartanion(0:3), Rot_mat(3,3)

    print*, "    #dbg no of atoms in analysing class: ", pdb%n_atoms

    write(*,*)
    write(*,'("INFORMATION> ANALYSIS")')
    
    write(*,'("    # READING REFERENCE")')
    !# Initialization
    call self%set_n_atoms(pdb%n_atoms)
    call self%initialize()
    S(:,:)     = 0.0d0; eigenval(:) = 0.0d0; eigenvec(:,:) = 0.0d0 
    Quartanion(:) = 0.0d0

    !# store coordinates of the reference structure
    specified%coords = pdb%get_coords(inputs%atom_type, inputs%res_fit_range)
    !call pdb%output_ref()

    !# translate coordinates of the reference to the origin.  
    write(*,"('      COM OF REFERENCE (BEFORE): ', 3f10.3)") & 
    COM(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))

    specified%coords  = Move_coords_to_origin(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))

    write(*,"('      COM Of REFERENCE (AFTER) :',3f10.3)") & 
    COM(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))
    !----------------------------------------------

    open(unit_no, file = filename, status="old")
    imodel = 0

    !***loop for models
    do 
      read(unit_no,*, end=1234) ATOM
      select case(ATOM)
      case("MODEL")
          imodel = imodel + 1
          write(*,*)
          write(*,"('MODEL NO>', i5 )") imodel

          call self%read_a_frame(pdb%n_atoms, unit_no)

          specified_in_snapshot%coords = self%get_coords(inputs%atom_type,inputs%res_fit_range)

          !***Analyze here***
          write(*,*) ""
          write(*,"('    INFORMATION> (0) RMSD CALCULATION')")
          call least_square_RMSD(specified%coords, specified_in_snapshot%coords)

          !PCA class will be added here.
          !Principal axis of moment of inertia will be added. 
          !******************
        case default
            print*, "not model line"
      end select
    enddo
    1234 print*, "Reading was finished." 
  end subroutine
end module
