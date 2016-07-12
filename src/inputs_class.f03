module input_module 
    implicit none

    type input_class
        logical             :: IsOutput
        integer             :: n_atoms, n_residues, n_chains
        real(8)             :: cellsize(1:3)
        character(len=80)   :: file_name_ref, file_name_trj, file_name_trj_list
        character(len=2)    :: atom_type
        integer             :: res_fit_range(2)

        character(len=80)   :: file_name_moment_ref

        logical            :: is_calc = .false.
        logical            :: is_trj  = .false.
        logical            :: is_trj_list  = .false.

    contains
        procedure :: set_inputs
    end type

contains

    !https://gist.github.com/ponderomotion/3833266
    SUBROUTINE split_string(instring, string1, string2, delim)
      CHARACTER(len=*) :: instring
      CHARACTER(len=*) :: delim
      CHARACTER(len=*),INTENT(OUT):: string1, string2
      INTEGER :: index

      instring = TRIM(instring)

      index = SCAN(instring,delim)
      string1 = adjustl(instring(1:index-1))
      string2 = adjustl(instring(index+1:))
    END SUBROUTINE split_string


    subroutine set_inputs(self)
      class(input_class) :: self
      character(len=100) :: line 
      character(len=10)  :: header
      character(len=10)  :: range(2)
      character(len=100) :: param, string  
      integer :: i, ios
      write(*,"('INFORMATION> INPUTS')")

      do 
        read(*,fmt='(a)', iostat = ios, end=123) line

        call split_string(line, header, param, ":")

!        print*, "Header : ", header
!        print*, "param  : ", param

        select case(header)

        case("natoms")
            read(param,*) self%n_atoms !translate char to int
            write(*, '("    No of atoms: ", i7)') self%n_atoms

        case("ref")
            self%file_name_ref = param
            write(*, '("    Reference : ", a)')  self%file_name_ref

        case("trj")
            if (self%is_trj_list .eqv. .true.) stop "trj"
            self%file_name_trj = param
            write(*, '("    Trajectry : ", a)')  self%file_name_trj 
            self%is_trj = .true.

        case("trj_list")
            if (self%is_trj .eqv. .true.) stop "trj"
            self%file_name_trj_list = param
            write(*, '("    Trajectry list : ", a)')  self%file_name_trj_list 
            self%is_trj_list = .true.

        case("rmsd")
            if (self%is_calc .eqv. .true.) stop "Multiple calculations are specified. Select only one calc."
            call split_string(param ,self%atom_type,string,",")
            call split_string(string,range(1),range(2),"-")
            read(range,*) self%res_fit_range
      
            write(*, '("    Atom type : ", a)')  self%atom_type
            write(*, '("    fit range : ", 2i7)')  self%res_fit_range
            self%is_calc = .true.

        case("moment_ref")
            if (self%is_calc .eqv. .true.) stop "Multiple calculations are specified. Select only one calc."
            self%file_name_moment_ref = param
            write(*, '("    Reference for moment of inertia: ", a)') self%file_name_moment_ref
            self%is_calc = .true.

        case default
            write(*, '("    An unrecognized line was detected &
                            which may not affect this analysis.")')
            cycle

        end select
      enddo
      123 continue
      !print*, "exited"
      !stop

      !write(*, fmt = '(a)', advance = 'no') "Make frames?(.true. / .false.): "
      !read(*, *)  self%IsOutput
      !write(*, *) self%IsOutput

      !write(*, fmt = '(a)', advance = 'no') "No. of residues               : "
      !read(*, *)  self%n_residues
      !write(*, *) self%n_residues

      !write(*, fmt = '(a)', advance = 'no') "No. of chains                 : "
      !read(*, *)  self%n_chains
      !write(*, *)  self%n_chains

      !write(*, fmt = '(a)', advance = 'no') "Cell size (xcell,ycell,zcell) : "
      !read(*, *)  self%cellsize(1:3)
      !write(*, *)  self%cellsize(1:3)


    end subroutine
end module
