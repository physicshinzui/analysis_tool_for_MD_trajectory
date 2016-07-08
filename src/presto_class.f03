module presto
  implicit none

    type var_snapshot
        real(4)                :: potential
        real(4), allocatable   :: x(:), y(:), z(:)
        integer                :: iconf
        real(8)                :: cellsize(3)
        integer(4)             :: istp,iyn15v,iyn15h
        real(4)                :: sitime, sec, totalene, kinetic, temperature, rmsf, rmsd
    end type

    type :: presto_class

    end type

    type, extends(presto_class) :: load_trj

    end type
contains

!  function set_snapshot() result(aaa) 
  
  
!  end function  
end module
