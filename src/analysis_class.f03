module analysis_module
  use variables, only : coordinates
  use input_module
  use list_inp_module
  use pdb_module
  use calculation
  use alignment_module
  implicit none
  character(len=6)  :: ATOM
  integer :: iii

  !???Inheritance is Ok????
  type, extends(pdb_class) :: FortTrj!analysis_class
  contains
    procedure :: analysing_in_loop
    procedure :: read_a_frame
  end type
contains
  
  subroutine read_a_frame(self, n_atoms, unit_no)
    class(FortTrj) :: self
    !class(analysis_class) :: self
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


  subroutine analysing_in_loop(self, filename, pdb, inputs)
    class(FortTrj)                :: self
    character(len=*) , intent(in) :: filename
    character(len=120)            :: filename_in_list
    type(pdb_class)  , intent(in) :: pdb
    type(input_class), intent(in) :: inputs
    type(list_inp_class)          :: linp
    type(coordinates)             :: ref, specified, snapshot, specified_in_snapshot, rotated
    integer                       :: unit_no = 1122,unit_no_list = 1123, iatom,imodel
    integer                       :: temp_var_location_eigen
    real(8)                       :: S(4,4), eigenval(4), eigenvec(4,4)
    real(8)                       :: Quartanion(0:3), Rot_mat(3,3)
   
    integer :: n_lists, ifile

    write(*,*)
    write(*,'("INFORMATION> ANALYSIS")')
    
    !# Initialization
    call self%set_n_atoms(pdb%n_atoms)
    call self%initialize()
    S(:,:)     = 0.0d0; eigenval(:) = 0.0d0; eigenvec(:,:) = 0.0d0 
    Quartanion(:) = 0.0d0

    !if (RMSDcalc == .true.) then 
        ! store reference coordinates which are used in RMSD calc 
        specified%coords = pdb%get_coords(inputs%atom_type, inputs%res_fit_range)

        ! translate the coordinates to the origin.  
        write(*,"('      COM OF REFERENCE (BEFORE): ', 3f10.3)") & 
        COM(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))

        specified%coords  = Move_coords_to_origin(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))

        write(*,"('      COM Of REFERENCE (AFTER) :',3f10.3)") & 
        COM(specified%coords(1,:),specified%coords(2,:),specified%coords(3,:))
    !endif

    open(unit_no_list, file = filename, status="old")
    n_lists = linp%count(inputs%file_name_trj_list)
    print*, "AAAAA@@@",n_lists

    imodel = 0
    listloop: do ifile = 1, n_lists
        
        read(unit_no_list,'(a)') filename_in_list
        open(unit_no, file=filename_in_list)
        write(*,'("-------------------------------------------")')
        write(*,'("FILE NAME: ", a)') trim(filename_in_list)

        modelLoop: do 
            read(unit_no,*, end=1234) ATOM
            select case(ATOM)
            case("MODEL")
                imodel = imodel + 1
                write(*,*)
                write(*,"('MODEL NO>', i5 )") imodel

                call self%read_a_frame(pdb%n_atoms, unit_no)

                specified_in_snapshot%coords = self%get_coords(inputs%atom_type,inputs%res_fit_range)

                !--------Analyze here------
                write(*,*) ""
                write(*,"('    INFORMATION> (0) RMSD CALCULATION')")
                call least_square_RMSD(specified%coords, specified_in_snapshot%coords)

                !PCA class will be added here.
                !-------------------------
            case default
                !print*, "not model line", ATOM
            end select
        enddo modelLoop

        1234 print*, "Reading was finished." 
        close(unit_no)

    enddo listLoop
  end subroutine
end module
