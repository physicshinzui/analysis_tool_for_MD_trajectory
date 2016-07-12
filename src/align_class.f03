module alignment_module 
  use linear_algebra_module, only : jacobi
  use calculation          , only : COM, Move_coords_to_origin
  implicit none
contains
  
  function MakeSymMat(rA,rB) result(S)
    integer :: i, j, k
    real(8), intent(in) :: rA(:,:), rB(:,:)
    real(8)  :: S(4,4)          !Symmetric matrix
    real(8), allocatable :: a(:,:), b(:,:)
    integer ::  n  

    n = size(rA,2)
    print*,"# dbg  Array size (rA) [no of coords] :", n
    allocate(a(3,n) ,b(3,n))

    do i = 1, n
      a(1,i) = rA(1,i) + rB(1,i)
      a(2,i) = rA(2,i) + rB(2,i)
      a(3,i) = rA(3,i) + rB(3,i)

      b(1,i) = rB(1,i) - rA(1,i)
      b(2,i) = rB(2,i) - rA(2,i)
      b(3,i) = rB(3,i) - rA(3,i)
    enddo

    !***Construct symmetric matrix
    S(:,:) = 0.0d0
    do j = 1, n
      S(1,1) = S(1,1) + b(1,j)**2 + b(2,j)**2 + b(3,j)**2
      S(2,1) = S(2,1) + a(3,j)*b(2,j) - a(2,j)*b(3,j)
      S(3,1) = S(3,1) - a(3,j)*b(1,j) + a(1,j)*b(3,j)
      S(4,1) = S(4,1) + a(2,j)*b(1,j) - a(1,j)*b(2,j)

      S(1,2) = S(1,2) + a(3,j)*b(2,j) - a(2,j)*b(3,j)
      S(2,2) = S(2,2) + b(1,j)**2 + a(2,j)**2 + a(3,j)**2
      S(3,2) = S(3,2) + b(1,j)*b(2,j) - a(1,j)*a(2,j)
      S(4,2) = S(4,2) + b(1,j)*b(3,j) - a(1,j)*a(3,j)

      S(1,3) = S(1,3) - a(3,j)*b(1,j) + a(1,j)*b(3,j)
      S(2,3) = S(2,3) + b(1,j)*b(2,j) - a(1,j)*a(2,j)
      S(3,3) = S(3,3) + a(1,j)**2 + b(2,j)**2 + a(3,j)**2 
      S(4,3) = S(4,3) + b(2,j)*b(3,j) - a(2,j)*a(3,j)

      S(1,4) = S(1,4) + a(2,j)*b(1,j) - a(1,j)*b(2,j)
      S(2,4) = S(2,4) + b(1,j)*b(3,j) - a(1,j)*a(3,j)
      S(3,4) = S(3,4) + b(2,j)*b(3,j) - a(2,j)*a(3,j)
      S(4,4) = S(4,4) + a(1,j)**2 + a(2,j)**2 + b(3,j)**2
    enddo
    S(:,:) = S(:,:)/n

!*  **Judge weather symmetric matrix or not
    do i = 1, 4
      do j = 1, 4
        if(S(j,i) == S(i,j)) then
!          write(*,'("   Symmetry! ","S(j,i): ", 2i4, f10.3, ", S(i,j): ", 2i4, f10.3)') i, j, S(j,i), j, i, S(i,j) 
        else
          write(*,'("No Symmetry! ","S(j,i): ", 2i4, f10.3, ", S(i,j): ", 2i4, f10.3)') i, j, S(j,i), j, i, S(i,j) 
          stop
        endif
      enddo
    enddo

  end function

  subroutine MakeRotationMat(q, R)
    integer :: i
    real(8), intent(in)  :: q(0:3)     !quartanion
    real(8), intent(out) :: R(1:3,1:3) !Rotation matrix
    R(1,1) = 2*q(0)**2   + 2*q(1)**2 - 1
    R(2,1) = 2*q(1)*q(2) + 2*q(0)*q(3)
    R(3,1) = 2*q(1)*q(3) - 2*q(0)*q(2)
    R(1,2) = 2*q(1)*q(2) - 2*q(0)*q(3)
    R(2,2) = 2*q(0)**2   + 2*q(2)**2 - 1
    R(3,2) = 2*q(2)*q(3) + 2*q(0)*q(1)
    R(1,3) = 2*q(1)*q(3) + 2*q(0)*q(2)
    R(2,3) = 2*q(2)*q(3) - 2*q(0)*q(1)
    R(3,3) = 2*q(0)**2   + 2*q(3)**2 - 1 
  end subroutine

  function rotate(rot_mat, coords) result(rotated_coords)
    real(8) :: coords(:,:)
    real(8) :: rot_mat(:,:)
    real(8), allocatable :: rotated_coords(:,:)
    allocate(rotated_coords(size(coords,1),size(coords,2)))
    rotated_coords = matmul(rot_mat, coords)
  end function

  function CalcRMSD(rA, rB) result(rmsd)
    real(8), intent(in) :: rA(:,:), rB(:,:)  !Coordinates of conformation A(<=target) and B(<=trajectry)
    real(8)             :: rmsd       
    real(8) :: xdiff, ydiff, zdiff, rdiff
    integer :: i, j
    
    xdiff = 0.0d0
    ydiff = 0.0d0
    zdiff = 0.0d0

    do i = 1, size(rA, 2)
      !if ( abs(rB(1,i)-rA(1,i)) > 0.00001 ) then
      !    print*,"does not match" 
      !    stop
      !endif
      xdiff = xdiff + (rB(1,i) - rA(1,i))**2
      ydiff = ydiff + (rB(2,i) - rA(2,i))**2
      zdiff = zdiff + (rB(3,i) - rA(3,i))**2
    enddo
    rdiff = xdiff + ydiff + zdiff

    rmsd = sqrt( rdiff / dble(size(rA, 2)) )
  end function

  subroutine show_sym_mat(S)
    real(8) :: S(4,4)
    !print*, "    # dbg size of symmetry matrix: ", size(S, 1), size(S, 2)
    write(*,*) ""
    write(*,"('    INFORMATION> (1) SYMMETRIC MATRIX (S)')")
    write(*,"(4x, 4f10.3)") S(1,1),S(1,2),S(1,3),S(1,4)
    write(*,"(4x, 4f10.3)") S(2,1),S(2,2),S(2,3),S(2,4)
    write(*,"(4x, 4f10.3)") S(3,1),S(3,2),S(3,3),S(3,4)
    write(*,"(4x, 4f10.3)") S(4,1),S(4,2),S(4,3),S(4,4)
  end subroutine
  
  subroutine show_eigen(eigenval, eigenvec)
    integer :: i
    real(8) :: eigenval(4), eigenvec(4,4)
    write(*,*) ""
    write(*,"('    INFORMATION> (2) EIGENVALUES AND EIGENVECTORS')")
    write(*,"('          LOC MIN EIGENVALE : ', i3)") minloc(eigenval)
    write(*,"('          MIN EIGENVALE     : ', f10.3)") eigenval(minloc(eigenval))
    write(*,"('          THE EIGENVECTOR   : ', 4f10.3)") eigenvec(1:4, minloc(eigenval)) 
    write(*,"('          <LIST OF EIGENVALUES>')" )
    write(*,"('          ', i3, f10.3)") (i, eigenval(i), i = 1, size(eigenval,1))  
    write(*,*) ""
    write(*,"('          <LIST OF EIGENVECTORS>')" )
    write(*,"(5x, 4f10.3)") eigenvec(1,1), eigenvec(1,2), eigenvec(1,3),eigenvec(1,4)
    write(*,"(5x, 4f10.3)") eigenvec(2,1), eigenvec(2,2), eigenvec(2,3),eigenvec(2,4)
    write(*,"(5x, 4f10.3)") eigenvec(3,1), eigenvec(3,2), eigenvec(3,3),eigenvec(3,4)
    write(*,"(5x, 4f10.3)") eigenvec(4,1), eigenvec(4,2), eigenvec(4,3),eigenvec(4,4)
  end subroutine

  subroutine show_rot_mat(Rot_mat)
    real(8), intent(in) :: Rot_mat(:,:)
    write(*,*) ""
    write(*,"('          <ROTATION MATRIX>')") 
    write(*,"(5x, 3f10.3)") Rot_mat(1,1), Rot_mat(1,2), Rot_mat(1,3)
    write(*,"(5x, 3f10.3)") Rot_mat(2,1), Rot_mat(2,2), Rot_mat(2,3)
    write(*,"(5x, 3f10.3)") Rot_mat(3,1), Rot_mat(3,2), Rot_mat(3,3)
  end subroutine

  !@@@@@@@@@@@@make this for alignment
  subroutine align(ref_coords, tag_coords)
    real(8), intent(in)  :: ref_coords(:,:), tag_coords(:,:) 
    real(8), allocatable :: tag_coords_orig(:,:), rotated_coords(:,:)
    real(8)              :: S(4,4)
    real(8)              :: eigenval(4), eigenvec(4,4)
    real(8)              :: Quartanion(0:3), Rot_mat(3,3)

  end subroutine

  subroutine least_square_RMSD(ref_coords, tag_coords)
    !***
    !Inputs coordinates size must be the same size.
    !***
    real(8), intent(in)  :: ref_coords(:,:), tag_coords(:,:) 
    real(8), allocatable :: tag_coords_orig(:,:), rotated_coords(:,:)
    real(8)              :: S(4,4)
    real(8)              :: eigenval(4), eigenvec(4,4)
    real(8)              :: Quartanion(0:3), Rot_mat(3,3)
    integer              :: unit_no = 111 

    open(111,file = "rmsd.dat")
    !open(111,file = "rmsd.dat", status="replace")
    write(*,"('       SIZE OF INPUT TARGET COORDINATES:',2i7 )") size(tag_coords,1), size(tag_coords,2)
    allocate(tag_coords_orig(size(tag_coords,1), size(tag_coords,2))) 
    allocate(rotated_coords( size(tag_coords,1), size(tag_coords,2)))
    
    write(*,"('       COM OF SPECIFIED ATOMS  (BEFORE):', 3f10.3)") & 
    COM(tag_coords(1,:), tag_coords(2,:), tag_coords(3,:) )
    
    tag_coords_orig = Move_coords_to_origin(tag_coords(1,:),tag_coords(2,:),tag_coords(3,:))
    
    write(*,"('       COM OF SPECIFIED ATOMS  (AFTER) :',3f10.3)") &
    COM(tag_coords_orig(1,:), tag_coords_orig(2,:), tag_coords_orig(3,:) ) 

    S = MakeSymMat(tag_coords_orig,ref_coords)  ; call show_sym_mat(S)

    call Jacobi(S,eigenval,eigenvec,4,4)        ; call show_eigen(eigenval, eigenvec)

    Quartanion(0:3) = eigenvec(1:4, minloc(eigenval,1))

    call MakeRotationMat(Quartanion, Rot_mat)   ; call show_rot_mat(Rot_mat)
    rotated_coords = rotate(Rot_mat, tag_coords_orig)
    write(*,"('     *RMSD: ', f10.3)") CalcRMSD(ref_coords, rotated_coords)
    write(111,*) CalcRMSD(ref_coords, rotated_coords)

  end subroutine

end module

