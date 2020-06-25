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
    mwpointer f1_pr,f2_pr,f3_pr
    mwpointer T_pr,P_pr,X_pr
    mwpointer A1_pr,A2_pr,A3_pr
    
    ! Array information
	mwSize :: m(4),n1,n2,n3,n4,nn1,nn2,nn3,nn4,n(4)
	mwSize ndim,coefs,mdim,mmax,nodes
	integer*4 myclassid
    double precision, allocatable, dimension(:,:,:,:) :: f1,f2,f3
    double precision, allocatable, dimension(:,:) :: T,P,X
	double precision, allocatable, dimension(:,:,:,:) :: A1,A2,A3
      
	! Load Inputs
	! Polynomial orders
    ndim = 4
    n1 = mxGetScalar(prhs(1))
    n2 = mxGetScalar(prhs(2))
    n3 = mxGetScalar(prhs(3))
    n4 = mxGetScalar(prhs(4))
    nn1 = n1 + 1
    nn2 = n2 + 1
    nn3 = n3 + 1
    nn4 = n4 + 1
    n = [nn1,nn2,nn3,nn4]
    coefs = n(1)*n(2)*n(3)*n(4)
    
    ! Original Functions
    mdim = 4
    m_pr = mxGetDimensions(prhs(5))  
    call mxCopyPtrToInteger4(m_pr,m,mdim)
    mmax = max(m(1),m(2),m(3),m(4))+1
    nodes = m(1)*m(2)*m(3)*m(4)
    allocate(f1(m(1),m(2),m(3),m(4))) 
    allocate(f2(m(1),m(2),m(3),m(4))) 
    allocate(f3(m(1),m(2),m(3),m(4)))
    f1_pr = mxGetPr(prhs(5))
    f2_pr = mxGetPr(prhs(6))
    f3_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(f1_pr,f1,nodes)
    call mxCopyPtrToReal8(f2_pr,f2,nodes) 
    call mxCopyPtrToReal8(f3_pr,f3,nodes) 
    
    ! Chebyshev Polynomial Parameters
    allocate(T(mmax,mmax))
    allocate(P(mmax,mmax))
    allocate(X(mmax,mmax))
    T_pr = mxGetPr(prhs(8))
    P_pr = mxGetPr(prhs(9))
    X_pr = mxGetPr(prhs(10))
    call mxCopyPtrToReal8(T_pr,T,mmax*mmax)
    call mxCopyPtrToReal8(P_pr,P,mmax*mmax)
    call mxCopyPtrToReal8(X_pr,X,mmax*mmax)
    
    !Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(A1(n(1),n(2),n(3),n(4))) 
    allocate(A2(n(1),n(2),n(3),n(4))) 
    allocate(A3(n(1),n(2),n(3),n(4)))
    plhs(1) = mxCreateNumericArray(ndim,n,myclassid,0)
    plhs(2) = mxCreateNumericArray(ndim,n,myclassid,0)
    plhs(3) = mxCreateNumericArray(ndim,n,myclassid,0)
    A1_pr = mxGetPr(plhs(1))
    A2_pr = mxGetPr(plhs(2))
    A3_pr = mxGetPr(plhs(3))
    
    ! Call subroutine for assignment
    call chebweights(A1,A2,A3,n(1),n(2),n(3),n(4),m(1),m(2),m(3),m(4),mmax,f1,f2,f3,T,P,X)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(A1,A1_pr,coefs)
    call mxCopyReal8ToPtr(A2,A2_pr,coefs)
    call mxCopyReal8ToPtr(A3,A3_pr,coefs)
    
    ! Deallocate arrays
    deallocate(f1)    
    deallocate(f2)    
    deallocate(f3)
    deallocate(T)
    deallocate(P)
    deallocate(X)
    deallocate(A1)
    deallocate(A2)
    deallocate(A3)
    
end subroutine mexFunction

subroutine chebweights(A1,A2,A3,nn1,nn2,nn3,nn4,m1,m2,m3,m4,mmax,f1,f2,f3,T,P,X)

    implicit none
    mwSize :: nn1,nn2,nn3,nn4,mmax,m1,m2,m3,m4
    double precision, dimension(nn1,nn2,nn3,nn4) :: A1,A2,A3
    double precision, dimension(m1,m2,m3,m4) :: f1,f2,f3
    double precision, dimension(mmax,mmax) :: T,P,X
    
    double precision, dimension(m1) :: x1_grid
    double precision, dimension(m2) :: x2_grid 
    double precision, dimension(m3) :: x3_grid 
    double precision, dimension(m4) :: x4_grid 
    double precision, dimension(nn1) :: w1
    double precision, dimension(nn2) :: w2
    double precision, dimension(nn3) :: w3
    double precision, dimension(nn4) :: w4
    double precision, dimension(mmax) :: T1,P1,X1
    double precision, dimension(mmax) :: T2,P2,X2
    double precision, dimension(mmax) :: T3,P3,X3
    double precision, dimension(mmax) :: T4,P4,X4
    double precision tt1,tt2,tt3,tt4,nestsum1,nestsum2,nestsum3,temp
    mwSize :: j1,j2,j3,j4,i1,i2,i3,i4
    
    ! Use zeros from desired order of approximation for grid
    x1_grid = X(m1,1:m1);
    x2_grid = X(m2,1:m2);
    x3_grid = X(m3,1:m3);
    x4_grid = X(m4,1:m4);
    
    ! Weights when calculating coefficients
    w1(1) = 1.0d0/m1;
    w1(2:nn1) = 2.0d0/m1;
    w2(1) = 1.0d0/m2;
    w2(2:nn2) = 2.0d0/m2;
    w3(1) = 1.0d0/m3;
    w3(2:nn3) = 2.0d0/m3;    
    w4(1) = 1.0d0/m4;
    w4(2:nn4) = 2.0d0/m4;
    
    ! Calculate coefficients of continuous least squares approx.
    do j4 = 1,nn4
      do j3 = 1,nn3
        do j2 = 1,nn2
          do j1 = 1,nn1
            nestsum1 = 0;
            nestsum2 = 0;
            nestsum3 = 0;
            T1 = T(j1,:);
            T2 = T(j2,:);
            T3 = T(j3,:);
            T4 = T(j4,:);
            P1 = P(j1,:);
            P2 = P(j2,:);
            P3 = P(j3,:);
            P4 = P(j4,:);
            do i4 = 1,m4
             do i3 = 1,m3
              do i2 = 1,m2
               do i1 = 1,m1
                ! Evaluate Chebyshev polynomials
                tt1 = sum(T1*x1_grid(i1)**P1);
                tt2 = sum(T2*x2_grid(i2)**P2);
                tt3 = sum(T3*x3_grid(i3)**P3);
                tt4 = sum(T4*x4_grid(i4)**P4);
                nestsum1 = nestsum1 + f1(i1,i2,i3,i4)*tt1*tt2*tt3*tt4;
                nestsum2 = nestsum2 + f2(i1,i2,i3,i4)*tt1*tt2*tt3*tt4;
                nestsum3 = nestsum3 + f3(i1,i2,i3,i4)*tt1*tt2*tt3*tt4;
               end do
              end do
             end do 
            end do
            temp = w1(j1)*w2(j2)*w3(j3)*w4(j4);
            A1(j1,j2,j3,j4) = temp*nestsum1;
            A2(j1,j2,j3,j4) = temp*nestsum2;
            A3(j1,j2,j3,j4) = temp*nestsum3;
          end do
        end do
      end do
    end do
    
end subroutine chebweights