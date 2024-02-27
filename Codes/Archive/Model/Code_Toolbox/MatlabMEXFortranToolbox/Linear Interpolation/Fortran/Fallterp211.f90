#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    ! Declarations
	implicit none

    ! mexFunction argument
    mwPointer plhs(*), prhs(*)    
    integer*4 nlhs, nrhs
    
    ! Function declarations
    mwSize mxGetN  
    mwpointer mxGetPr, mxCreateNumericArray, mxGetDimensions 
    double precision mxGetScalar  
    integer*4 mxClassIDFromClassName 
    
    ! Pointers to input/output mxArrays
    mwpointer e2_pr
    mwpointer x1_pr,x2_pr
    mwpointer x2i_pr
    mwpointer pf1_pr
    mwpointer o1_pr
    
    ! Array information
	mwSize nx1,nx2,nodes,e2
	integer*4 myclassid
	double precision, allocatable, dimension(:) :: x1,x2,x2i
	double precision x1i
	double precision, allocatable, dimension(:,:) :: pf1
	double precision, allocatable, dimension(:) :: o1
        
    ! Load Inputs
    ! Grids
    nx1 = mxGetN(prhs(1))  
    nx2 = mxGetN(prhs(2))  
    nodes = nx1*nx2
    allocate(x1(nx1))
    allocate(x2(nx2))
    x1_pr = mxGetPr(prhs(1))
    x2_pr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(x1_pr,x1,nx1)
    call mxCopyPtrToReal8(x2_pr,x2,nx2)
    
    ! Point to evaluate
    x1i = mxGetScalar(prhs(3))
    e2_pr = mxGetDimensions(prhs(4))  
    call mxCopyPtrToInteger4(e2_pr,e2,1)
    allocate(x2i(e2))
    x2i_pr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(x2i_pr,x2i,e2)
    
    ! Rules
    allocate(pf1(nx1,nx2))
    pf1_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(e2)) 
    plhs(1) = mxCreateNumericArray(1,e2,myclassid,0)
    o1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allterp(o1,nx1,nx2,x1,x2,x1i,x2i,pf1,e2)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,e2)
    
    ! Deallocate arrays
    deallocate(x1)
    deallocate(x2)  
    deallocate(x2i)     
    deallocate(pf1)   
    deallocate(o1)

end subroutine mexFunction

subroutine allterp(o1,nx1,nx2,x1,x2,x1i,x2i,pf1,e2)

    implicit none
    mwSize :: nx1,nx2,e2
    mwSize i2
    double precision x1i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx1,nx2) :: pf1
    double precision, dimension(e2) :: x2i
    double precision, dimension(e2) :: o1

    ! Interpolate
    do i2 = 1,e2
        o1(i2) = fastrap(nx1,nx2,x1,x2,x1i,x2i(i2),pf1)
    end do
end subroutine allterp


function fastrap(nx1,nx2,x1,x2,x1i,x2i,z)

    implicit none
    mwSize nx1,nx2
    double precision :: x1i,x2i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx1,nx2) :: z
    
    double precision :: s1, s2
    double precision :: x1i_min, x2i_min
    mwSize loc1, loc2
    double precision, dimension(2) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2
    mwSize m1, m2
    
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