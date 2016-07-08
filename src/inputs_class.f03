module input_module 
    implicit none

    type input_class
        logical           :: IsOutput
        integer           :: n_atoms, n_residues, n_chains
        real(8)           :: cellsize(1:3)
        character(len=80) :: file_name_ref, file_name_trj
        character(len=2)  :: atom_type
        integer           :: res_fit_range(2)
    contains
        procedure :: set_inputs
    end type

contains

    subroutine set_inputs(self)
      class(input_class) :: self
      character(len=30)   :: header  
      write(*,"('INFORMATION> INPUTS')")

      !do 
      !  read(*,*) header
      !  print*, header
      !  stop
      !  select case(header)
      !  case("--output-frame")
      !    write(*, fmt = '(a)', advance = 'no') "Make frames?(.true. / .false.): "
      !    read(*, *)  self%IsOutput
      !    write(*, *) self%IsOutput

      !  end select
      !enddo
      !write(*, fmt = '(a)', advance = 'no') "Make frames?(.true. / .false.): "
      !read(*, *)  self%IsOutput
      !write(*, *) self%IsOutput

      write(*, fmt = '(a)', advance = 'no') "No. of atoms                  : "
      read(*, *)  self%n_atoms
      write(*, *) self%n_atoms

      !write(*, fmt = '(a)', advance = 'no') "No. of residues               : "
      !read(*, *)  self%n_residues
      !write(*, *) self%n_residues

      !write(*, fmt = '(a)', advance = 'no') "No. of chains                 : "
      !read(*, *)  self%n_chains
      !write(*, *)  self%n_chains

      !write(*, fmt = '(a)', advance = 'no') "Cell size (xcell,ycell,zcell) : "
      !read(*, *)  self%cellsize(1:3)
      !write(*, *)  self%cellsize(1:3)

      write(*, fmt = '(a)', advance = 'no') "Reference PDB                 : "
      read(*, "(a80)")  self%file_name_ref
      write(*, *)  self%file_name_ref

      write(*, fmt = '(a)', advance = 'no') "Trajectory file               : "
      read(*, *)  self%file_name_trj
      write(*, *)  self%file_name_trj

      write(*, fmt = '(a)', advance = 'no') "Atom_type                     : "
      read(*, *)  self%atom_type
      write(*, *)  self%atom_type

      write(*, fmt = '(a)', advance = 'no') "Residue range for fitiing     : "
      read(*, *)  self%res_fit_range(1), self%res_fit_range(2)
      write(*, *)  self%res_fit_range
    end subroutine
end module
