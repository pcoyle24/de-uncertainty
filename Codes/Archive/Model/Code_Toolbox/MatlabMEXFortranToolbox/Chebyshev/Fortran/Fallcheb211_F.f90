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
    mwpointer m_pr,n_pr,c_pr
    mwpointer z1_bnd_pr,z2_bnd_pr
    mwpointer z1i_pr,z2i_pr
    mwpointer A1_pr,A2_pr
    mwpointer T_pr,P_pr
    mwpointer f1_pr,f2_pr
    
    ! Array information
	mwSize :: m(3),n(2),c(2)
	mwSize mdim,ndim,mmax
    mwSize nodes,coefs,chebs
	integer*4 myclassid
	double precision, dimension(2) :: z1_bnd,z2_bnd
	double precision, allocatable, dimension(:,:,:) :: z1i,z2i	
	double precision, allocatable, dimension(:,:) :: A1 
	double precision, allocatable, dimension(:,:) :: T,P   
    double precision, allocatable, dimension(:,:,:) :: f1
    
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    z2_bnd_pr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    call mxCopyPtrToReal8(z2_bnd_pr,z2_bnd,2)
    
    ! Points to evaluate
    mdim = 3
    m_pr = mxGetDimensions(prhs(3))  
    call mxCopyPtrToInteger4(m_pr,m,mdim)
    nodes = m(1)*m(2)*m(3) 
    allocate(z1i(m(1),m(2),m(3))) 
    allocate(z2i(m(1),m(2),m(3)))        
    z1i_pr = mxGetPr(prhs(3))
    z2i_pr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(z1i_pr,z1i,nodes)
    call mxCopyPtrToReal8(z2i_pr,z2i,nodes) 

    ! Least square weights
    ndim = 2
    n_pr = mxGetDimensions(prhs(5))
    call mxCopyPtrToInteger4(n_pr,n,ndim)
    coefs = n(1)*n(2)
    allocate(A1(n(1),n(2)))
    A1_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(A1_pr,A1,coefs)    
    
    ! Chebyshev Polynomial Parameters
    c_pr = mxGetDimensions(prhs(6))  
    call mxCopyPtrToInteger4(c_pr,c,2)
    chebs = c(1)*c(2)
    allocate(T(c(1),c(2)))
    allocate(P(c(1),c(2)))
    T_pr = mxGetPr(prhs(6))
    P_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(T_pr,T,chebs)
    call mxCopyPtrToReal8(P_pr,P,chebs)

    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(m(1),m(2),m(3)))
    plhs(1) = mxCreateNumericArray(mdim,m,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allcheb(f1,nodes,z1i,z2i,z1_bnd,z2_bnd,c(1),c(2),n(1),n(2),m(1),m(2),m(3),T,P,A1)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,nodes)
    
    ! Deallocate arrays
    deallocate(z1i)    
    deallocate(z2i)
    deallocate(A1)     
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)
    
end subroutine mexFunction

subroutine allcheb(f1,nodes,z1i,z2i,z1_bnd,z2_bnd,c1,c2,nn1,nn2,m1,m2,m3,T,P,A1)

    implicit none
    mwSize :: nodes,c1,c2,nn1,nn2,m1,m2,m3
    double precision, dimension(m1,m2,m3) :: z1i,z2i,f1
    double precision, dimension(2) :: z1_bnd,z2_bnd
    double precision, dimension(c1,c2) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1
    double precision, dimension(nn2,nn2) :: T2,P2    
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn2) :: vec2
    double precision, dimension(nn1,nn2) :: A1
    double precision :: temp1,x1i,x2i
    mwSize :: i1,i2,i3,j1,j2
      
    T1 = T(1:nn1,1:nn1)
    T2 = T(1:nn2,1:nn2)   
    P1 = P(1:nn1,1:nn1)
    P2 = P(1:nn2,1:nn2)
    do i3 = 1,m3
     do i2 = 1,m2
      do i1 = 1,m1
        ! Transform Z in [a,b] to X in [-1,1]
        x1i = 2*(z1i(i1,i2,i3)-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1;
        x2i = 2*(z2i(i1,i2,i3)-z2_bnd(1))/(z2_bnd(2) - z2_bnd(1)) - 1;

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
        f1(i1,i2,i3) = temp1
      end do
     end do
    end do
    
end subroutine allcheb