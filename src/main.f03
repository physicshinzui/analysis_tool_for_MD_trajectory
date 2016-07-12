program main
  use input_module
  use pdb_module
  use analysis_module
  type(input_class)    :: ic 
  type(pdb_class)      :: ref_pdb
  type(fortTrj)        :: ana 
  real(8)              :: t1, t2
 
  call cpu_time(t1)

  call ic%set_inputs()

  call ref_pdb%set_n_atoms(ic%n_atoms)
  call ref_pdb%get_structure(ic%file_name_ref)
  !call ref_pdb%output_ref()

!  write(*,*) "START> Moment of inertial section"
!  call calc_inertial_axis(ref_pdb, ic%res_fit_range) 
!  write(*,*) "END  > Moment of inertial section"
!  stop

  !***substituting ref_pdb, into analysing_in_loop function.
  call ana%analysing_in_loop(ic%file_name_trj_list, ref_pdb, ic)

  call cpu_time(t2)
  write(*,'("#######################################")')
  write(*,'("Compuational time: ", f10.3)') t2 -t1
  write(*,'("#######################################")')
end program

