program main
  use input_module
  use pdb_module
  use analysis_module
  type(input_class)    :: ic 
  type(pdb_class)      :: ref_pdb
  type(analysis_class) :: ana
 
  call ic%set_inputs()

  call ref_pdb%set_n_atoms(ic%n_atoms)
  call ref_pdb%get_structure(ic%file_name_ref)
  !call ref_pdb%output_ref()

  !***substitute instance var, pdb, into this.
  call ana%analysing_in_loop(ic%file_name_trj, ref_pdb, ic)

end program

