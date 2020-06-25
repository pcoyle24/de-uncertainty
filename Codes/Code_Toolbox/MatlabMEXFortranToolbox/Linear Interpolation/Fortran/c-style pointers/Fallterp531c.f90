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
    mwsize :: mxCreateDoubleMatrix, mxGetM, mxGetN
    double precision :: mxGetScalar
    
    ! Inputs
	mwsize nx1,nx2,nx3,nx4,nx5,ne5
	integer, pointer :: x1,x2,x3,x4,x5
	double precision x1i, x2i, x3i, x4i
    integer, pointer :: x5i
	integer, pointer :: r1,r2,r3
	
	! Outputs
	integer, pointer :: o1,o2,o3
           
    ! Load Inputs
    ! Number of Points
    nx1 = mxGetN(prhs(1))
    nx2 = mxGetN(prhs(2))
    nx3 = mxGetN(prhs(3))
    nx4 = mxGetN(prhs(4))               
    nx5 = mxGetN(prhs(5))    
    
    ! Grid Values
    x1 => mxGetPr(prhs(1))
    x2 => mxGetPr(prhs(2))
    x3 => mxGetPr(prhs(3))
    x4 => mxGetPr(prhs(4))    
    x5 => mxGetPr(prhs(5)) 
    
    ! Point to evaluate
    x1i = mxGetScalar(prhs(6))
    x2i = mxGetScalar(prhs(7))
    x3i = mxGetScalar(prhs(8))   
    x4i = mxGetScalar(prhs(9)) 
    ne5 = mxGetM(prhs(10))
    x5i => mxGetPr(prhs(10))
    
    ! Rules
    r1 => mxGetPr(prhs(11))
    r2 => mxGetPr(prhs(12))
    r3 => mxGetPr(prhs(13))      
    
    !Create array for return argument
    plhs(1) = mxCreateDoubleMatrix(ne5,1,0)
    plhs(2) = mxCreateDoubleMatrix(ne5,1,0)
    plhs(3) = mxCreateDoubleMatrix(ne5,1,0)    
    o1 => mxGetPr(plhs(1))
    o2 => mxGetPr(plhs(2))
    o3 => mxGetPr(plhs(3))      
    
    ! Call subroutine for assignment
    call allterp(o1,o2,o3,nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i,r1,r2,r3,ne5)

end subroutine mexFunction

subroutine allterp(o1,o2,o3,nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i,r1,r2,r3,ne5)

    implicit none
    mwsize nx1,nx2,nx3,nx4,nx5,ne5
    mwsize i5
    double precision x1i,x2i,x3i,x4i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4    
    double precision, dimension(nx5) :: x5
    double precision, dimension(nx1,nx2,nx3,nx4,nx5) :: r1,r2,r3
    double precision, dimension(ne5) :: x5i
    double precision, dimension(ne5) :: o1,o2,o3

    ! Interpolate
    do i5 = 1,ne5
        o1(i5) = fastrap(nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i(i5),r1)
        o2(i5) = fastrap(nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i(i5),r2)
        o3(i5) = fastrap(nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i(i5),r3)              
    end do
end subroutine allterp


function fastrap(nx1,nx2,nx3,nx4,nx5,x1,x2,x3,x4,x5,x1i,x2i,x3i,x4i,x5i,z)

    implicit none
    mwsize nx1,nx2,nx3,nx4,nx5
    double precision :: x1i,x2i,x3i,x4i,x5i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx5) :: x5
    double precision, dimension(nx1,nx2,nx3,nx4,nx5) :: z
    
    double precision :: s1, s2, s3, s4, s5
    double precision :: x1i_min, x2i_min, x3i_min, x4i_min, x5i_min
    mwsize loc1, loc2, loc3, loc4, loc5
    double precision, dimension(5) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3, w4, w5
    mwsize m1, m2, m3, m4, m5
    
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
    
    s5 = x5(2) - x5(1)
    x5i_min = x5i - x5(1)
    loc5 = min(nx5-1,max(1,floor(x5i_min/s5) + 1));        
    
    xi = [x1i, x2i, x3i, x4i, x5i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3), x4(loc4), x5(loc5)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1), x4(loc4+1), x5(loc5+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    w4 = [w_1(4), w_2(4)]
    w5 = [w_1(5), w_2(5)]
    
    fastrap = 0
    
    do m5 = 0, 1
     do m4 = 0, 1
      do m3 = 0, 1
       do m2 = 0, 1
        do m1 = 0, 1
         fastrap = fastrap + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5)
        end do
       end do
      end do
     end do
    end do

end function fastrap