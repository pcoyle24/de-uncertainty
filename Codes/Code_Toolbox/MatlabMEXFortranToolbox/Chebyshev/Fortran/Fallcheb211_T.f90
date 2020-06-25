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
    mwpointer m_pr,n_pr,e2_pr
    mwpointer z1_bnd_pr,z2_bnd_pr
    mwpointer z2i_pr
    mwpointer A1_pr,A2_pr
    mwpointer T_pr,P_pr
    mwpointer f1_pr,f2_pr
    
    ! Array information
	mwSize e2,m(2),n(2)
	mwSize mdim,ndim,edim
    mwSize coefs
	integer*4 myclassid
	double precision, dimension(2) :: z1_bnd,z2_bnd
	double precision :: z1i
	double precision, allocatable, dimension(:) :: z2i,f1
	double precision, allocatable, dimension(:,:) :: A1
	double precision, allocatable, dimension(:,:) :: T,P
	       
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    z2_bnd_pr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    call mxCopyPtrToReal8(z2_bnd_pr,z2_bnd,2)
    
    ! Points to evaluate
    z1i = mxGetScalar(prhs(3))
    edim = 1
    e2_pr = mxGetDimensions(prhs(4))
    call mxCopyPtrToInteger4(e2_pr,e2,edim)
    allocate(z2i(e2))
    z2i_pr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(z2i_pr,z2i,e2)
  
    ! Least square weights
    ndim = 2
    n_pr = mxGetDimensions(prhs(5))
    call mxCopyPtrToInteger4(n_pr,n,ndim)
    coefs = n(1)*n(2)
    allocate(A1(n(1),n(2)))
    A1_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(A1_pr,A1,coefs)
    
    ! Chebyshev Polynomial Parameters
    mdim = 2
    m_pr = mxGetDimensions(prhs(6))  
    call mxCopyPtrToInteger4(m_pr,m,mdim)
    allocate(T(m(1),m(2)))
    allocate(P(m(1),m(2)))
    T_pr = mxGetPr(prhs(6))
    P_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(T_pr,T,m(1)*m(2))
    call mxCopyPtrToReal8(P_pr,P,m(1)*m(2))
    
    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(e2))
    plhs(1) = mxCreateNumericArray(edim,e2,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allcheb(f1,e2,z1i,z2i,z1_bnd,z2_bnd,m(1),m(2),n(1),n(2),T,P,A1)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,e2)
    
    ! Deallocate arrays
    deallocate(z2i)
    deallocate(A1)  
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)

end subroutine mexFunction

subroutine allcheb(f1,e2,z1i,z2i,z1_bnd,z2_bnd,m1,m2,nn1,nn2,T,P,A1)

    implicit none
    mwSize :: e2,m1,m2,nn1,nn2
    double precision :: z1i,x1i,x2i
    double precision, dimension(e2) :: z2i,f1
    double precision, dimension(2) :: z1_bnd,z2_bnd
    double precision, dimension(m1,m2) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1
    double precision, dimension(nn2,nn2) :: T2,P2
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn2) :: vec2
    double precision, dimension(nn1,nn2) :: A1
    double precision :: temp1,temp2
    mwSize :: j,j1,j2
    
    T1 = T(1:nn1,1:nn1);
    T2 = T(1:nn2,1:nn2);
    
    P1 = P(1:nn1,1:nn1);
    P2 = P(1:nn2,1:nn2);
    
    ! Transform Z in [a,b] to X in [-1,1]
    x1i = 2*(z1i-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1;
    do j = 1,e2
	    x2i = 2*(z2i(j)-z2_bnd(1))/(z2_bnd(2) - z2_bnd(1)) - 1;
	    
        ! Evaluate Chebyshev Polynomials
        vec1 = sum(T1*x1i**P1,2);
        vec2 = sum(T2*x2i**P2,2);
        
        ! Evaluate function
        temp1 = 0
        do j2 = 1,nn2
          do j1 = 1,nn1
              temp1 = temp1 + A1(j1,j2)*vec1(j1)*vec2(j2) 
          end do
        end do
        f1(j) = temp1
    end do
    
end subroutine allcheb