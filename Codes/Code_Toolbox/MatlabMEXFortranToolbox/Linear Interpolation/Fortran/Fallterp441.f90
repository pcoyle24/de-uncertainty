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
    mwpointer e4_pr
    mwpointer x1_pr,x2_pr,x3_pr,x4_pr
    mwpointer x4i_pr
    mwpointer pf1_pr,pf2_pr,pf3_pr,pf4_pr
    mwpointer o1_pr,o2_pr,o3_pr,o4_pr
    
    ! Array information
	mwSize nx1,nx2,nx3,nx4,nodes,e4
	integer*4 myclassid
	double precision, allocatable, dimension(:) :: x1,x2,x3,x4,x4i
	double precision x1i,x2i,x3i
	double precision, allocatable, dimension(:,:,:,:) :: pf1,pf2,pf3,pf4
	double precision, allocatable, dimension(:) :: o1,o2,o3,o4       
       
    ! Load Inputs
    ! Grids
    nx1 = mxGetN(prhs(1))  
    nx2 = mxGetN(prhs(2))  
    nx3 = mxGetN(prhs(3))
    nx4 = mxGetN(prhs(4))
    nodes = nx1*nx2*nx3*nx4
    allocate(x1(nx1))
    allocate(x2(nx2))
    allocate(x3(nx3))
    allocate(x4(nx4))
    x1_pr = mxGetPr(prhs(1))
    x2_pr = mxGetPr(prhs(2))
    x3_pr = mxGetPr(prhs(3))
    x4_pr = mxGetPr(prhs(4))    
    call mxCopyPtrToReal8(x1_pr,x1,nx1)
    call mxCopyPtrToReal8(x2_pr,x2,nx2)
    call mxCopyPtrToReal8(x3_pr,x3,nx3)
    call mxCopyPtrToReal8(x4_pr,x4,nx4)  
    ! Point to evaluate
    x1i = mxGetScalar(prhs(5))
    x2i = mxGetScalar(prhs(6))
    x3i = mxGetScalar(prhs(7))  
    e4_pr = mxGetDimensions(prhs(8))  
    call mxCopyPtrToInteger4(e4_pr,e4,1)
    allocate(x4i(e4))
    x4i_pr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(x4i_pr,x4i,e4)
      
    ! Rules
    allocate(pf1(nx1,nx2,nx3,nx4))
    allocate(pf2(nx1,nx2,nx3,nx4))
    allocate(pf3(nx1,nx2,nx3,nx4))
    allocate(pf4(nx1,nx2,nx3,nx4))
    pf1_pr = mxGetPr(prhs(9))
    pf2_pr = mxGetPr(prhs(10))
    pf3_pr = mxGetPr(prhs(11))
    pf4_pr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    call mxCopyPtrToReal8(pf2_pr,pf2,nodes)
    call mxCopyPtrToReal8(pf3_pr,pf3,nodes)
    call mxCopyPtrToReal8(pf4_pr,pf4,nodes)
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(e4)) 
    allocate(o2(e4))
    allocate(o3(e4))
    allocate(o4(e4))
    plhs(1) = mxCreateNumericArray(1,e4,myclassid,0)
    plhs(2) = mxCreateNumericArray(1,e4,myclassid,0)
    plhs(3) = mxCreateNumericArray(1,e4,myclassid,0)
    plhs(4) = mxCreateNumericArray(1,e4,myclassid,0)    
    o1_pr = mxGetPr(plhs(1))
    o2_pr = mxGetPr(plhs(2))
    o3_pr = mxGetPr(plhs(3))
    o4_pr = mxGetPr(plhs(4))
    
    ! Call subroutine for assignment
    call allterp(o1,o2,o3,o4,nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,pf1,pf2,pf3,pf4,e4)
        
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,e4)
    call mxCopyReal8ToPtr(o2,o2_pr,e4)
    call mxCopyReal8ToPtr(o3,o3_pr,e4)
    call mxCopyReal8ToPtr(o4,o4_pr,e4)    
    
    ! Deallocate arrays
    deallocate(x1)
    deallocate(x2)  
    deallocate(x3) 
    deallocate(x4) 
    deallocate(x4i)     
    deallocate(pf1)  
    deallocate(pf2)  
    deallocate(pf3)
    deallocate(pf4)
    deallocate(o1)
    deallocate(o2)
    deallocate(o3)
    deallocate(o4)
    
end subroutine mexFunction

subroutine allterp(o1,o2,o3,o4,nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,pf1,pf2,pf3,pf4,e4)

    implicit none
    mwSize :: nx1,nx2,nx3,nx4,e4
    mwSize i4
    double precision x1i,x2i,x3i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx1,nx2,nx3,nx4) :: pf1,pf2,pf3,pf4
    double precision, dimension(e4) :: x4i
    double precision, dimension(e4) :: o1,o2,o3,o4

    ! Interpolate
    do i4 = 1,e4
        o1(i4) = fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i(i4),pf1)
        o2(i4) = fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i(i4),pf2)
        o3(i4) = fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i(i4),pf3)
        o4(i4) = fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i(i4),pf4)
    end do
end subroutine allterp

function fastrap(nx1,nx2,nx3,nx4,x1,x2,x3,x4,x1i,x2i,x3i,x4i,z)

    implicit none
    mwSize nx1,nx2,nx3,nx4
    double precision :: x1i,x2i,x3i,x4i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx4) :: x4
    double precision, dimension(nx1,nx2,nx3,nx4) :: z
    
    double precision :: s1, s2, s3, s4
    double precision :: x1i_min, x2i_min, x3i_min, x4i_min
    mwSize loc1, loc2, loc3, loc4
    double precision, dimension(4) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3, w4
    mwSize m1, m2, m3, m4
    
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