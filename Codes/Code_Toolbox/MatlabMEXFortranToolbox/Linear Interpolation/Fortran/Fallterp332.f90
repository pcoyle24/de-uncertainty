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
    mwpointer e2_pr,e3_pr
    mwpointer x1_pr,x2_pr,x3_pr
    mwpointer x2i_pr,x3i_pr
    mwpointer pf1_pr,pf2_pr,pf3_pr
    mwpointer o1_pr,o2_pr,o3_pr
    
    ! Array information
    mwSize nx1,nx2,nx3,nodes,e2,e3
    integer*4 myclassid
    double precision, allocatable, dimension(:) :: x1,x2,x3,x2i,x3i
    double precision x1i
    double precision, allocatable, dimension(:,:,:) :: pf1,pf2,pf3
    double precision, allocatable, dimension(:,:) :: o1,o2,o3   
    
    ! Load Inputs
    ! Grids
    nx1 = mxGetN(prhs(1))  
    nx2 = mxGetN(prhs(2))  
    nx3 = mxGetN(prhs(3))
    nodes = nx1*nx2*nx3
    allocate(x1(nx1))
    allocate(x2(nx2))
    allocate(x3(nx3))
    x1_pr = mxGetPr(prhs(1))
    x2_pr = mxGetPr(prhs(2))
    x3_pr = mxGetPr(prhs(3))    
    call mxCopyPtrToReal8(x1_pr,x1,nx1)
    call mxCopyPtrToReal8(x2_pr,x2,nx2)
    call mxCopyPtrToReal8(x3_pr,x3,nx3) 
    ! Point to evaluate
    x1i = mxGetScalar(prhs(4))     
    e2_pr = mxGetDimensions(prhs(5))   
    e3_pr = mxGetDimensions(prhs(6)) 
    call mxCopyPtrToInteger4(e2_pr,e2,1) 
    call mxCopyPtrToInteger4(e3_pr,e3,1)
    allocate(x2i(e2))
    allocate(x3i(e3))
    x2i_pr = mxGetPr(prhs(5))
    x3i_pr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(x2i_pr,x2i,e2)
    call mxCopyPtrToReal8(x3i_pr,x3i,e3)

    ! Rules
    allocate(pf1(nx1,nx2,nx3))
    allocate(pf2(nx1,nx2,nx3))
    allocate(pf3(nx1,nx2,nx3))
    pf1_pr = mxGetPr(prhs(7))
    pf2_pr = mxGetPr(prhs(8))
    pf3_pr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(pf1_pr,pf1,nodes)
    call mxCopyPtrToReal8(pf2_pr,pf2,nodes)
    call mxCopyPtrToReal8(pf3_pr,pf3,nodes)

    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(o1(e2,e3)) 
    allocate(o2(e2,e3)) 
    allocate(o3(e2,e3)) 
    plhs(1) = mxCreateNumericArray(2,[e2,e3],myclassid,0)
    plhs(2) = mxCreateNumericArray(2,[e2,e3],myclassid,0)
    plhs(3) = mxCreateNumericArray(2,[e2,e3],myclassid,0)
    o1_pr = mxGetPr(plhs(1))
    o2_pr = mxGetPr(plhs(2))
    o3_pr = mxGetPr(plhs(3))
                  
    ! Call subroutine for assignment
    call allterp(o1,o2,o3,nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,pf1,pf2,pf3,e2,e3)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(o1,o1_pr,e2*e3)
    call mxCopyReal8ToPtr(o2,o2_pr,e2*e3)
    call mxCopyReal8ToPtr(o3,o3_pr,e2*e3)
    
    ! Deallocate arrays
    deallocate(x1)
    deallocate(x2)  
    deallocate(x3) 
    deallocate(x2i)     
    deallocate(x3i)     
    deallocate(pf1)     
    deallocate(pf2)     
    deallocate(pf3)          
    deallocate(o1)           
    deallocate(o2)           
    deallocate(o3)   

end subroutine mexFunction

subroutine allterp(o1,o2,o3,nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,pf1,pf2,pf3,e2,e3)

    implicit none
    mwSize :: nx1,nx2,nx3,e2,e3
    mwSize i2,i3
    double precision x1i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx1,nx2,nx3) :: pf1,pf2,pf3
    double precision, dimension(e2) :: x2i
    double precision, dimension(e3) :: x3i
    double precision, dimension(e2,e3) :: o1,o2,o3

    ! Interpolate
    do i3 = 1,e3
      do i2 = 1,e2
        o1(i2,i3) = fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i(i2),x3i(i3),pf1)
        o2(i2,i3) = fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i(i2),x3i(i3),pf2)
        o3(i2,i3) = fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i(i2),x3i(i3),pf3)
      end do
    end do
end subroutine allterp


function fastrap(nx1,nx2,nx3,x1,x2,x3,x1i,x2i,x3i,z)

    implicit none
    mwSize nx1,nx2,nx3
    double precision :: x1i,x2i,x3i,x4i,fastrap
    double precision, dimension(nx1) :: x1
    double precision, dimension(nx2) :: x2
    double precision, dimension(nx3) :: x3
    double precision, dimension(nx1,nx2,nx3) :: z
    
    double precision :: s1, s2, s3
    double precision :: x1i_min, x2i_min, x3i_min
    mwSize loc1, loc2, loc3
    double precision, dimension(3) :: xi, xi_left, xi_right, w_2, w_1
    double precision, dimension(2) :: w1, w2, w3, w4
    mwSize m1, m2, m3
    
    s1 = x1(2) - x1(1)
    x1i_min = x1i - x1(1)
    loc1 = min(nx1-1,max(1,floor(x1i_min/s1) + 1));
    
    s2 = x2(2) - x2(1)
    x2i_min = x2i - x2(1)
    loc2 = min(nx2-1,max(1,floor(x2i_min/s2) + 1));
    
    s3 = x3(2) - x3(1)
    x3i_min = x3i - x3(1)
    loc3 = min(nx3-1,max(1,floor(x3i_min/s3) + 1));           

    xi = [x1i, x2i, x3i]
    xi_left = [x1(loc1), x2(loc2), x3(loc3)]
    xi_right = [x1(loc1+1), x2(loc2+1), x3(loc3+1)]

    w_2 = (xi - xi_left)/(xi_right - xi_left)
    w_1 = 1 - w_2
    w1 = [w_1(1), w_2(1)]
    w2 = [w_1(2), w_2(2)]
    w3 = [w_1(3), w_2(3)]
    
    fastrap = 0
   
    do m3 = 0, 1
     do m2 = 0, 1
      do m1 = 0, 1  
       fastrap = fastrap + w1(m1+1)*w2(m2+1)*w3(m3+1)*z(loc1+m1,loc2+m2,loc3+m3)
      end do  
     end do
    end do
    
end function fastrap
