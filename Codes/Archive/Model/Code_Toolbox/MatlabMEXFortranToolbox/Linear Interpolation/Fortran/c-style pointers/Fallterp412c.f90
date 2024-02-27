#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)    
    ! Declarations
	implicit none

    ! mexFunction argument
    integer*4, intent(in) :: nlhs, nrhs
    mwpointer, intent(in), dimension(*) :: prhs
    mwpointer, intent(out), dimension(*) :: plhs
    
    ! Function declarations
    integer, pointer :: mxGetPr
    mwsize :: mxCreateDoubleMatrix,mxGetM,mxGetN
    double precision :: mxGetScalar
    
    ! Inputs
	mwsize nx1,nx2,nx3,nx4,ne3,ne4
	integer, pointer :: x1,x2,x3,x4
	double precision x1i, x2i
    integer, pointer :: x3i,x4i
	integer, pointer :: r1
	
	! Outputs
	integer, pointer :: o1
           
    ! Load Inputs
    ! Number of Points
    nx1 = mxGetN(prhs(1))
    nx2 = mxGetN(prhs(2))
    nx3 = mxGetN(prhs(3))
    nx4 = mxGetN(prhs(4))
    
    ! Grid Values
    x1 => mxGetPr(prhs(1))
    x2 => mxGetPr(prhs(2))
    x3 => mxGetPr(prhs(3))
    x4 => mxGetPr(prhs(4))
    
    ! Point to evaluate
    x1i = mxGetScalar(prhs(5))
    x2i = mxGetScalar(prhs(6))
    ne3 = mxGetM(prhs(7))
    x3i => mxGetPr(prhs(7))
    ne4 = mxGetM(prhs(8))
    x4i => mxGetPr(prhs(8))
    
    ! Rules
    r1 => mxGetPr(prhs(9))
    
    !Create array for return argument
    plhs(1) = mxCreateDoubleMatrix(ne3,ne4,0)   
    o1 => mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allterp(o1,nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,r1,ne3,ne4)

end subroutine mexFunction

subroutine allterp(o1,nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,r1,ne3,ne4)

    implicit none
    mwsize nx1,nx2,nx3,nx4,ne3,ne4
    mwsize i3, i4
    double precision x1i,x2i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx1,nx2,nx3,nx4) :: r1
    double precision, dimension(ne3) :: x3i
    double precision, dimension(ne4) :: x4i
    double precision, dimension(ne3,ne4) :: o1

    ! Interpolate
    do i4 = 1,ne4
      do i3 = 1,ne3    
        o1(i3,i4) = fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i(i3),x4i(i4),r1)
      end do  
    end do
end subroutine allterp


function fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,z)

    implicit none
    mwsize nx1,nx2,nx3,nx4
    double precision :: x1i,x2i,x3i,x4i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx1,nx2,nx3,nx4) :: z
    
    double precision :: s1, s2, s3, s4
    double precision :: x1i_min, x2i_min, x3i_min, x4i_min
    mwsize loc1, loc2, loc3, loc4
    double precision, dimension(4) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3, w4
    mwsize m1, m2, m3, m4
    
    s1 = x1(2) - x1(1)
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    s2 = x2(2) - x2(1)
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    s3 = x3(2) - x3(1)
    x3i_min = x3i - x3(1)
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));
    
    s4 = x4(2) - x4(1)
    x4i_min = x4i - x4(1)
    loc4 = min(nx4-1,max(1,floor(x4i_min/s4) + 1));       
    
    xi = [x1i, x2i, x3i, x4i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3), x4(loc4)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1), x4(loc4+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    w4 = [w_1(4), w_2(4)]
    
    fastrap = 0
    
    do m4 = 0, 1
     do m3 = 0, 1
      do m2 = 0, 1
       do m1 = 0, 1
        fastrap = fastrap + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4)
       end do
      end do
     end do
    end do

end function fastrap