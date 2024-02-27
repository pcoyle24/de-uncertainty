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
    mwpointer m_pr
    mwpointer f1_pr,f2_pr
    mwpointer T_pr,P_pr,X_pr
    mwpointer A1_pr,A2_pr
    
    ! Array information
	mwSize :: m(2),n1,m1,n2,m2,nn1,nn2,n(2)
	mwSize ndim,coefs,mdim,mmax,nodes
	integer*4 myclassid
    double precision, allocatable, dimension(:,:) :: f1
    double precision, allocatable, dimension(:,:) :: T,P,X
	double precision, allocatable, dimension(:,:) :: A1
     
	! Load Inputs
	! Polynomial orders
    ndim = 2
    n1 = mxGetScalar(prhs(1))
    n2 = mxGetScalar(prhs(2))
    nn1 = n1 + 1
    nn2 = n2 + 1
    n = [nn1, nn2]
    coefs = n(1)*n(2)
    
    ! Original Functions
    mdim = 2
    m1 = mxGetScalar(prhs(3))
    m2 = mxGetScalar(prhs(4))
    m = [m1, m2]
    mmax = max(m(1),m(2))+1
    nodes = m(1)*m(2)
    allocate(f1(m(1),m(2)))
    f1_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(f1_pr,f1,nodes) 
    
    ! Chebyshev Polynomial Parameters
    allocate(T(mmax,mmax))
    allocate(P(mmax,mmax))
    allocate(X(mmax,mmax))
    T_pr = mxGetPr(prhs(6))
    P_pr = mxGetPr(prhs(7))
    X_pr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(T_pr,T,mmax*mmax)
    call mxCopyPtrToReal8(P_pr,P,mmax*mmax)
    call mxCopyPtrToReal8(X_pr,X,mmax*mmax)
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(A1(n(1),n(2)))
    plhs(1) = mxCreateNumericArray(ndim,n,myclassid,0)
    A1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call chebweights(A1,n(1),n(2),m(1),m(2),mmax,f1,T,P,X)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(A1,A1_pr,coefs)
    
    ! Deallocate arrays
    deallocate(f1)
    deallocate(T)  
    deallocate(P)     
    deallocate(X)
    deallocate(A1)
    
end subroutine mexFunction

subroutine chebweights(A1,nn1,nn2,m1,m2,mmax,f1,T,P,X)

    implicit none
    mwSize :: nn1,nn2,m1,m2,mmax
    double precision, dimension(nn1,nn2) :: A1  
    double precision, dimension(m1,m2) :: f1
    double precision, dimension(mmax,mmax) :: T,P,X
    
    double precision, dimension(m1) :: x1_grid
    double precision, dimension(m2) :: x2_grid  
    double precision, dimension(nn1) :: w1
    double precision, dimension(nn2) :: w2
    double precision, dimension(mmax) :: T1,P1,X1
    double precision, dimension(mmax) :: T2,P2,X2
    double precision tt1,tt2,nestsum1,temp
    mwSize :: j1,j2,i1,i2
    
    ! Use zeros from desired order of approximation for grid
    x1_grid = X(m1,1:m1);
    x2_grid = X(m2,1:m2);
    
    ! Weights when calculating coefficients
    w1(1) = 1.0d0/m1;
    w1(2:nn1) = 2.0d0/m1;
    w2(1) = 1.0d0/m2;
    w2(2:nn2) = 2.0d0/m2;    
    
    ! Calculate coefficients of continuous least squares approx.
    do j2 = 1,nn2
       do j1 = 1,nn1
            nestsum1 = 0;
            T1 = T(j1,:);
            T2 = T(j2,:);
            P1 = P(j1,:);
            P2 = P(j2,:);
            do i2 = 1,m2
              do i1 = 1,m1
                ! Evaluate Chebyshev polynomials
                tt1 = sum(T1*x1_grid(i1)**P1);
                tt2 = sum(T2*x2_grid(i2)**P2);
                nestsum1 = nestsum1 + f1(i1,i2)*tt1*tt2;
              end do 
            end do
            temp = w1(j1)*w2(j2);
            A1(j1,j2) = temp*nestsum1;
        end do
    end do
    
end subroutine chebweights