module calculation 
  implicit none
contains

  !#@@@@@move COM calc to the alignment module
  function COM(x,y,z) 
    integer :: NoOfData 
    real(8) :: COM(3)
    real(8),intent(in)  :: x(:),y(:),z(:)
    real(8) :: SumX,SumY,SumZ 
    NoOfData = size(x)
    SumX = sum(x) 
    SumY = sum(y) 
    SumZ = sum(z) 
    COM(1) = SumX/NoOfData
    COM(2) = SumY/NoOfData
    COM(3) = SumZ/NoOfData
  end function 

  function Move_coords_to_origin(x,y,z) result(r)
    double precision,intent(in) :: x(:), y(:), z(:)
    double precision            :: CenterOfMass(3)
    real(8), allocatable :: r(:,:) 
    allocate(r(3, size(x)) )
    CenterOfMass = COM(x,y,z)
    !print*, "#dbg            Total size is ", size(r(:,:))
    !print*, "#dbg    2nd component size is ", size(r(:,1))
    !print*, "#dbg    3rd component size is ", size(r(1,:))
    r(1,:) = x(:) - CenterOfMass(1)
    r(2,:) = y(:) - CenterOfMass(2)
    r(3,:) = z(:) - CenterOfMass(3)
  end function

  !------------------------------
  subroutine Calc_RadiusOfGyration(n,x,y,z,COM, Rg)
  !  1st(INPUT) : N of atoms 
  !  5th(INPUT) : center of mass
  !  6th(OUTPUT): Radius of gyration
  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n), z(n)
  double precision, intent(in) :: COM(3) 
  double precision, intent(out) :: Rg
  double precision :: xtmp(n), ytmp(n), ztmp(n)
  integer :: i, j
  double precision :: ti, tf

  call cpu_time(ti)
  !Initialize
  Rg      = 0.0d0
  xtmp(:) = 0.0d0
  ytmp(:) = 0.0d0
  ztmp(:) = 0.0d0

  !(1) COM is subtructed from each coordinate (r - <r>) 
  xtmp(:) = x(:) - COM(1)
  ytmp(:) = y(:) - COM(2)
  ztmp(:) = z(:) - COM(3)

  !(2) ∑ _i^Natom (ri - <r>)^2
  do i = 1, n  
    Rg = Rg + xtmp(i)**2 + ytmp(i)**2 + ztmp(i)**2
  enddo


  !(3) 1/N \sigma(ri - r_av)^2
  Rg = sqrt(Rg/n)

  call cpu_time(tf)

  end subroutine
!------------------------------

!------------------------------
  subroutine Calc_EndToEnd(NumberOfAtomPair,x,y,z,TotalAtoms,ETE)
    integer,          intent(in)  :: NumberOfAtomPair(2) !Specify pair of atom number 
    integer,          intent(in)  :: TotalAtoms          !total atoms you read 
    double precision, intent(in)  :: x(:),y(:),z(:)
    double precision, intent(out) :: ETE
    double precision :: dx, dy, dz
    

    !(1) Calculate difference bet. ri- rj
    dx = x(NumberOfAtomPair(1)) - x(NumberOfAtomPair(2)) 
    dy = y(NumberOfAtomPair(1)) - y(NumberOfAtomPair(2)) 
    dz = z(NumberOfAtomPair(1)) - z(NumberOfAtomPair(2)) 
    ETE = sqrt(dx**2 + dy**2 + dz**2) 
  end subroutine


end module
