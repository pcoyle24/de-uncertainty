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
	mwsize nx1,nx2,ne2
	integer, pointer :: x1,x2
	double precision x1i
    integer, pointer :: x2i
	integer, pointer :: pf1
	
	! Outputs
	integer, pointer :: o1
           
    ! Load Inputs
    ! Number of Points
    nx1 = mxGetN(prhs(1))
    nx2 = mxGetN(prhs(2))
    
    ! Grid Values
    x1 => mxGetPr(prhs(1))
    x2 => mxGetPr(prhs(2))    
    
    ! Point to evaluate
    x1i = mxGetScalar(prhs(3))    
    ne2 = mxGetM(prhs(4))
    x2i => mxGetPr(prhs(4))
    
    ! Rules
    pf1 => mxGetPr(prhs(5)) 
    
    !Create array for return argument
    plhs(1) = mxCreateDoubleMatrix(ne2,1,0) 
    o1 => mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allterp(o1,nx1,nx2,x1,x2,x1i,x2i,pf1,ne2)

end subroutine mexFunction

subroutine allterp(o1,nx1,nx2,x1,x2,x1i,x2i,pf1,ne2)

    implicit none
    mwsize nx1,nx2,ne2
    mwsize i2
    double precision x1i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2 
    double precision, dimension(nx1,nx2) :: pf1
    double precision, dimension(ne2) :: x2i
    double precision, dimension(ne2) :: o1,o2

    ! Interpolate
    do i2 = 1,ne2
        o1(i2) = fastrap(nx1,nx2,x1,x2,x1i,x2i(i2),pf1)
    end do
end subroutine allterp


function fastrap(nx1,nx2,x1,x2,x1i,x2i,z)

    implicit none
    mwsize nx1,nx2
    double precision :: x1i,x2i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx1,nx2) :: z
    
    double precision :: s1, s2
    double precision :: x1i_min, x2i_min
    mwsize loc1, loc2
    double precision, dimension(2) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2
    mwsize m1, m2
    
    s1 = x1(2) - x1(1)
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    s2 = x2(2) - x2(1)
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    xi = [x1i, x2i]
    xi_left = [x1(loc1), x2(loc2)]
    xi_right = [x1(loc1+1), x2(loc2+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    
    fastrap = 0
    
    do m2 = 0, 1
      do m1 = 0, 1
        fastrap = fastrap + w1(m1+1)*w2(m2+1)*z(loc1+m1,loc2+m2)
      end do
    end do

end function fastrap