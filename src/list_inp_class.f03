module list_inp_module
  implicit none

  type list_inp_class
  contains
      procedure :: count
  end type
contains

  function count(self, filename) result(n_list)
      class(list_inp_class)       :: self
      character(len=*),intent(in) :: filename
      integer                     :: n_list
      call system("wc -l " // filename //" > noOfLists")
      open(11, file="noOfLists", status="old")
      read(11, *) n_list
      call system("rm noOfLists")
  end function

  subroutine open(self,filename)
    class(list_inp_class)       :: self
    character(len=*),intent(in) :: filename



  end subroutine
end module
