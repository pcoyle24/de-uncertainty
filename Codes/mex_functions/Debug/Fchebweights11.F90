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
    mwpointer f1_pr
    mwpointer T_pr,P_pr,X_pr
    mwpointer A1_pr
    
    ! Array information
		mwSize :: m,n1,nn1,n
		mwSize ndim,coefs,mdim,mmax,nodes
		integer*4 myclassid
	    double precision, allocatable, dimension(:) :: f1
	    double precision, allocatable, dimension(:,:) :: T,P,X
		double precision, allocatable, dimension(:) :: A1
    
	! Load Inputs
	! Polynomial orders
    ndim = 1
    n1 = mxGetScalar(prhs(1))
    nn1 = n1 + 1
    n = nn1
    coefs = n
    
    ! Original Functions
    mdim = 1
    m = mxGetScalar(prhs(2))
    mmax = m + 1
    nodes = m
    !if (nodes == 11) then
    !	call mexErrMsgTxt ( 'nodes == 11' )
    !end if
    allocate(f1(m))
    f1_pr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(f1_pr,f1,nodes)
   
    
    
    ! Chebyshev Polynomial Parameters
    allocate(T(mmax,mmax))
    allocate(P(mmax,mmax))
    allocate(X(mmax,mmax))
    T_pr = mxGetPr(prhs(4))
    P_pr = mxGetPr(prhs(5))
    X_pr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(T_pr,T,mmax*mmax)
    call mxCopyPtrToReal8(P_pr,P,mmax*mmax)
    call mxCopyPtrToReal8(X_pr,X,mmax*mmax)
    
    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(A1(n))
    plhs(1) = mxCreateNumericArray(ndim,n,myclassid,0)
    A1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call chebweights(A1,n,m,mmax,f1,T,P,X)

    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(A1,A1_pr,coefs)
    
    ! Deallocate arrays
    deallocate(f1)
    deallocate(T)  
    deallocate(P)     
    deallocate(X)
    deallocate(A1)
    
end subroutine mexFunction

subroutine chebweights(A1,nn1,m1,mmax,f1,T,P,X)

    implicit none
    mwSize :: nn1,m1,mmax
    double precision, dimension(nn1) :: A1  
    double precision, dimension(m1) :: f1
    double precision, dimension(mmax,mmax) :: T,P,X
    
    double precision, dimension(m1) :: x1_grid
    double precision, dimension(nn1) :: w1
    double precision, dimension(mmax) :: T1,P1,X1
    double precision tt1,nestsum1,temp
    mwSize :: j1,i1
    
    ! Use zeros from desired order of approximation for grid
    x1_grid = X(m1,1:m1);
    
    ! Weights when calculating coefficients
    w1(1) = 1.0d0/m1;
    w1(2:nn1) = 2.0d0/m1; 
    
    ! Calculate coefficients of continuous least squares approx.
     do j1 = 1,nn1
          nestsum1 = 0;
          T1 = T(j1,:);
          P1 = P(j1,:);
          do i1 = 1,m1
            ! Evaluate Chebyshev polynomials
            tt1 = sum(T1*x1_grid(i1)**P1);
            nestsum1 = nestsum1 + f1(i1)*tt1;
          end do 
          temp = w1(j1);
          A1(j1) = temp*nestsum1;
      end do
    
end subroutine chebweights