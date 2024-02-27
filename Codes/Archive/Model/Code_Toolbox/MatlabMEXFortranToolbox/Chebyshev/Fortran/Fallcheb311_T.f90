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
    mwpointer m_pr,n_pr,e3_pr
    mwpointer z1_bnd_pr,z2_bnd_pr,z3_bnd_pr
    mwpointer z3i_pr
    mwpointer A1_pr
    mwpointer T_pr,P_pr
    mwpointer f1_pr
    
    ! Array information
	mwSize e3,m(2),n(3)
	mwSize mdim,ndim,edim
    mwSize coefs
	integer*4 myclassid
	double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd
	double precision :: z1i,z2i
	double precision, allocatable, dimension(:) :: z3i,f1
	double precision, allocatable, dimension(:,:,:) :: A1
	double precision, allocatable, dimension(:,:) :: T,P
	       
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    z2_bnd_pr = mxGetPr(prhs(2))
    z3_bnd_pr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    call mxCopyPtrToReal8(z2_bnd_pr,z2_bnd,2)
    call mxCopyPtrToReal8(z3_bnd_pr,z3_bnd,2)
    
    ! Points to evaluate
    z1i = mxGetScalar(prhs(4))
    z2i = mxGetScalar(prhs(5))
    edim = 1
    e3_pr = mxGetDimensions(prhs(6))
    call mxCopyPtrToInteger4(e3_pr,e3,edim)
    allocate(z3i(e3))
    z3i_pr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(z3i_pr,z3i,e3)
  
    ! Least square weights
    ndim = 3
    n_pr = mxGetDimensions(prhs(7))
    call mxCopyPtrToInteger4(n_pr,n,ndim)
    coefs = n(1)*n(2)*n(3)
    allocate(A1(n(1),n(2),n(3)))
    A1_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(A1_pr,A1,coefs)
    
    ! Chebyshev Polynomial Parameters
    mdim = 2
    m_pr = mxGetDimensions(prhs(8))  
    call mxCopyPtrToInteger4(m_pr,m,mdim)
    allocate(T(m(1),m(2)))
    allocate(P(m(1),m(2)))
    T_pr = mxGetPr(prhs(8))
    P_pr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(T_pr,T,m(1)*m(2))
    call mxCopyPtrToReal8(P_pr,P,m(1)*m(2))
    
    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(e3))
    plhs(1) = mxCreateNumericArray(edim,e3,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allcheb(f1,e3,z1i,z2i,z3i,z1_bnd,z2_bnd,z3_bnd,m(1),m(2),n(1),n(2),n(3),T,P,A1)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,e3)
    
    ! Deallocate arrays
    deallocate(z3i)
    deallocate(A1)  
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)

end subroutine mexFunction

subroutine allcheb(f1,e3,z1i,z2i,z3i,z1_bnd,z2_bnd,z3_bnd,m1,m2,nn1,nn2,nn3,T,P,A1)

    implicit none
    mwSize :: e3,m1,m2,nn1,nn2,nn3
    double precision :: z1i,z2i,x1i,x2i,x3i
    double precision, dimension(e3) :: z3i,f1
    double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd
    double precision, dimension(m1,m2) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1
    double precision, dimension(nn2,nn2) :: T2,P2
    double precision, dimension(nn3,nn3) :: T3,P3
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn2) :: vec2
    double precision, dimension(nn3) :: vec3
    double precision, dimension(nn1,nn2,nn3) :: A1
    double precision :: temp1
    mwSize :: j,j1,j2,j3
    
    T1 = T(1:nn1,1:nn1);
    T2 = T(1:nn2,1:nn2);
    T3 = T(1:nn3,1:nn3);
    
    P1 = P(1:nn1,1:nn1);
    P2 = P(1:nn2,1:nn2);
    P3 = P(1:nn3,1:nn3);
    
    ! Transform Z in [a,b] to X in [-1,1]
    x1i = 2*(z1i-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1;
    x2i = 2*(z2i-z2_bnd(1))/(z2_bnd(2) - z2_bnd(1)) - 1;
    do j = 1,e3
	    x3i = 2*(z3i(j)-z3_bnd(1))/(z3_bnd(2) - z3_bnd(1)) - 1;
	    
        ! Evaluate Chebyshev Polynomials
        vec1 = sum(T1*x1i**P1,2);
        vec2 = sum(T2*x2i**P2,2);
        vec3 = sum(T3*x3i**P3,2);
        
        ! Evaluate function
        temp1 = 0
        do j3 = 1,nn3
          do j2 = 1,nn2
            do j1 = 1,nn1
              temp1 = temp1 + A1(j1,j2,j3)*vec1(j1)*vec2(j2)*vec3(j3)
            end do
          end do
        end do
        f1(j) = temp1
    end do
    
end subroutine allcheb