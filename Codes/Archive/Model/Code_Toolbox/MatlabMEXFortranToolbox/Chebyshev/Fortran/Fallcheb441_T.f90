#include "fintrf.h"
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    ! Declarations
	implicit none

    ! mexFunction argument
    mwPointer plhs(*), prhs(*)    
    integer*4 nlhs, nrhs
    
    ! Function declarations
    mwSize mxGetNumberofDimensions  
    mwpointer mxGetPr, mxCreateNumericArray, mxGetDimensions 
    double precision mxGetScalar  
    integer*4 mxClassIDFromClassName     
    
    ! Pointers to input/output mxArrays
    mwpointer m_pr,n_pr,e4_pr
    mwpointer z1_bnd_pr,z2_bnd_pr,z3_bnd_pr,z4_bnd_pr
    mwpointer z4i_pr
    mwpointer A1_pr,A2_pr,A3_pr,A4_pr
    mwpointer T_pr,P_pr
    mwpointer f1_pr,f2_pr,f3_pr,f4_pr
    
    ! Array information
	mwSize e4,m(2),n(4)
	mwSize mdim,ndim,edim
    mwSize coefs
	integer*4 myclassid
	double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd,z4_bnd
	double precision :: z1i,z2i,z3i
	double precision, allocatable, dimension(:) :: z4i,f1,f2,f3,f4
	double precision, allocatable, dimension(:,:,:,:) :: A1,A2,A3,A4 
	double precision, allocatable, dimension(:,:) :: T,P
	       
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    z2_bnd_pr = mxGetPr(prhs(2))
    z3_bnd_pr = mxGetPr(prhs(3))
    z4_bnd_pr = mxGetPr(prhs(4))   
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    call mxCopyPtrToReal8(z2_bnd_pr,z2_bnd,2)
    call mxCopyPtrToReal8(z3_bnd_pr,z3_bnd,2)
    call mxCopyPtrToReal8(z4_bnd_pr,z4_bnd,2)   
    
    ! Points to evaluate
    z1i = mxGetScalar(prhs(5))
    z2i = mxGetScalar(prhs(6))
    z3i = mxGetScalar(prhs(7))   
    edim = 1
    e4_pr = mxGetDimensions(prhs(8))
    call mxCopyPtrToInteger4(e4_pr,e4,edim)
    allocate(z4i(e4))
    z4i_pr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(z4i_pr,z4i,e4)
  
    ! Least square weights
    ndim = 4
    n_pr = mxGetDimensions(prhs(9))
    call mxCopyPtrToInteger4(n_pr,n,ndim)
    coefs = n(1)*n(2)*n(3)*n(4)
    allocate(A1(n(1),n(2),n(3),n(4))) 
    allocate(A2(n(1),n(2),n(3),n(4))) 
    allocate(A3(n(1),n(2),n(3),n(4))) 
    allocate(A4(n(1),n(2),n(3),n(4))) 
    A1_pr = mxGetPr(prhs(9))
    A2_pr = mxGetPr(prhs(10))
    A3_pr = mxGetPr(prhs(11))
    A4_pr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(A1_pr,A1,coefs)
    call mxCopyPtrToReal8(A2_pr,A2,coefs)
    call mxCopyPtrToReal8(A3_pr,A3,coefs)
    call mxCopyPtrToReal8(A4_pr,A4,coefs)
    
    ! Chebyshev Polynomial Parameters
    mdim = 2
    m_pr = mxGetDimensions(prhs(13))  
    call mxCopyPtrToInteger4(m_pr,m,mdim)    
    allocate(T(m(1),m(2)))
    allocate(P(m(1),m(2)))
    T_pr = mxGetPr(prhs(13))
    P_pr = mxGetPr(prhs(14))
    call mxCopyPtrToReal8(T_pr,T,m(1)*m(2))
    call mxCopyPtrToReal8(P_pr,P,m(1)*m(2))
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(e4)) 
    allocate(f2(e4))
    allocate(f3(e4))
    allocate(f4(e4))
    plhs(1) = mxCreateNumericArray(edim,e4,myclassid,0)
    plhs(2) = mxCreateNumericArray(edim,e4,myclassid,0)
    plhs(3) = mxCreateNumericArray(edim,e4,myclassid,0)
    plhs(4) = mxCreateNumericArray(edim,e4,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    f2_pr = mxGetPr(plhs(2))
    f3_pr = mxGetPr(plhs(3))
    f4_pr = mxGetPr(plhs(4))
    
    ! Call subroutine for assignment
    call allcheb(f1,f2,f3,f4,e4,z1i,z2i,z3i,z4i,z1_bnd,z2_bnd,z3_bnd,z4_bnd,m(1),m(2),n(1),n(2),n(3),n(4),T,P,A1,A2,A3,A4)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,e4)
    call mxCopyReal8ToPtr(f2,f2_pr,e4)
    call mxCopyReal8ToPtr(f3,f3_pr,e4)
    call mxCopyReal8ToPtr(f4,f4_pr,e4)
    
    ! Deallocate arrays
    deallocate(z4i)
    deallocate(A1)  
    deallocate(A2)
    deallocate(A3)
    deallocate(A4)     
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)
    deallocate(f2)
    deallocate(f3)
    deallocate(f4)
    
end subroutine mexFunction

subroutine allcheb(f1,f2,f3,f4,e4,z1i,z2i,z3i,z4i,z1_bnd,z2_bnd,z3_bnd,z4_bnd,m1,m2,nn1,nn2,nn3,nn4,T,P,A1,A2,A3,A4)

    implicit none
    mwSize :: e4,m1,m2,nn1,nn2,nn3,nn4
    double precision :: z1i,z2i,z3i,x1i,x2i,x3i,x4i
    double precision, dimension(e4) :: z4i,f1,f2,f3,f4
    double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd,z4_bnd  
    double precision, dimension(m1,m2) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1
    double precision, dimension(nn2,nn2) :: T2,P2
    double precision, dimension(nn3,nn3) :: T3,P3
    double precision, dimension(nn4,nn4) :: T4,P4
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn2) :: vec2 
    double precision, dimension(nn3) :: vec3 
    double precision, dimension(nn4) :: vec4 
    double precision, dimension(nn1,nn2,nn3,nn4) :: A1,A2,A3,A4
    double precision :: temp1,temp2,temp3,temp4
    mwSize :: j,j1,j2,j3,j4
    
    T1 = T(1:nn1,1:nn1);
    T2 = T(1:nn2,1:nn2);
    T3 = T(1:nn3,1:nn3);
    T4 = T(1:nn4,1:nn4);
    
    P1 = P(1:nn1,1:nn1);
    P2 = P(1:nn2,1:nn2);
    P3 = P(1:nn3,1:nn3);
    P4 = P(1:nn4,1:nn4);
    
    ! Transform Z in [a,b] to X in [-1,1]
    x1i = 2*(z1i-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1;
    x2i = 2*(z2i-z2_bnd(1))/(z2_bnd(2) - z2_bnd(1)) - 1;
    x3i = 2*(z3i-z3_bnd(1))/(z3_bnd(2) - z3_bnd(1)) - 1;	
    do j = 1,e4
	    x4i = 2*(z4i(j)-z4_bnd(1))/(z4_bnd(2) - z4_bnd(1)) - 1;
	    
        ! Evaluate Chebyshev Polynomials
        vec1 = sum(T1*x1i**P1,2);
        vec2 = sum(T2*x2i**P2,2);
        vec3 = sum(T3*x3i**P3,2);
        vec4 = sum(T4*x4i**P4,2);
        
        ! Evaluate function
        temp1 = 0
        temp2 = 0
        temp3 = 0
        temp4 = 0
        do j4 = 1,nn4
          do j3 = 1,nn3
            do j2 = 1,nn2
              do j1 = 1,nn1                
                  temp1 = temp1 + A1(j1,j2,j3,j4)*vec1(j1)*vec2(j2)*vec3(j3)*vec4(j4)
                  temp2 = temp2 + A2(j1,j2,j3,j4)*vec1(j1)*vec2(j2)*vec3(j3)*vec4(j4)
                  temp3 = temp3 + A3(j1,j2,j3,j4)*vec1(j1)*vec2(j2)*vec3(j3)*vec4(j4)
                  temp4 = temp4 + A4(j1,j2,j3,j4)*vec1(j1)*vec2(j2)*vec3(j3)*vec4(j4)                
              end do
            end do    
          end do
        end do
        f1(j) = temp1
        f2(j) = temp2
        f3(j) = temp3
        f4(j) = temp4
    end do
    
end subroutine allcheb