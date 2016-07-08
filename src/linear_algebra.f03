module linear_algebra_module
   implicit none

contains

   subroutine Jacobi(a,e,v,n,nn)
     integer :: i, j, aaa
     integer,intent(in) :: n, nn              !dimension of array (a, e, and v)
     integer :: p, q, kaisuu                  !Subscripts of a 
     integer, parameter :: kmax = 30000       !Max val for calculation of Jacobi loop
     !real(8), intent(inout)  :: a(:,:)      !Symmetric matrix
     real(8), intent(inout)  :: a(nn,nn)      !Symmetric matrix
     real(8), intent(out) :: e(nn), v(n,n)    !eigen val and vector
     real(8) :: eps = 0.000001, bunbo,R,T,S,C
     real(8) :: apq, apqmax
     real(8) :: aip, aiq, apj, aqj, vip, viq
   
   !***Judge weather symmetric matrix or not
     do i = 1, n
       do j = 1, n
         if (a(j,i) /= a(i,j)) then 
          !print*,i,j, a(j,i), a(i,j)
           stop "???Non-symmetric matrix was inputed.???"
         endif
       enddo
     enddo
   
   !***Construct Unit matrix
     do i = 1, N
       do j = 1, N
         v(i,j) = 0.0d0
       enddo
       v(i,i) = 1.0d0
     enddo
   
   !***Iteration of Jacobi
     do kaisuu = 1, kmax
       do p = 1, N -1
         do q = p+1, N
   !***generate rotation angle
           bunbo = a(p,p) - a(q,q)
           if (bunbo /= 0.0d0) then
             R = 2.0d0 * a(p,q)/bunbo
             T = 0.5d0 * atan(R)
           else
             T = 0.78539818d0  !0.5 * (0.5*pi) <= lmit of arctan.   
           endif
           s = sin(T)
           c = cos(T)
   !***Multiplication from left
           do j = 1, n
             apj = a(p,j)
             aqj = a(q,j)
             a(p,j) = apj*c + aqj*s
             a(q,j) = -apj*s + aqj *c
           enddo
   !***Multiplication from right
           do i = 1, n
             aip = a(i,p)
             aiq = a(i,q)
             a(i,p) =  aip*c + aiq*s
             a(i,q) = -aip*s + aiq*c
             vip = v(i,p)
             viq = v(i,q)
             v(i,p) =  vip*c + viq*s
             v(i,q) = -vip*s + viq*c
           enddo
         enddo
       enddo
   !***Judgement of convergence
       apqmax = 0.0d0
       do p = 1, n-1
         do q = p+1, n
           apq = abs(a(p,q))
           if(apq > apqmax) apqmax = apq
         enddo
       enddo
       if(apqmax <= eps) goto 100 
     enddo
   100 continue
   !write(*,'("N of loops(Jacobi): ",i5)') kaisuu
   !***trace eigenvalue
     forall (i = 1:n) e(i) = a(i,i)

   end subroutine


end module
