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
    mwpointer z1_bnd_pr,z2_bnd_pr,z3_bnd_pr
    mwpointer z1i_pr,z2i_pr,z3i_pr
    mwpointer A1_pr
    mwpointer T_pr,P_pr
    mwpointer f1_pr
    
    ! Array information
	mwSize :: m1,m2,m3,m(3),n1,n2,n3,nn1,nn2,nn3,n(3)
	mwSize mdim,ndim,c
    mwSize nodes,coefs,chebs
	integer*4 myclassid
	double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd
	double precision, allocatable, dimension(:,:,:) :: z1i,z2i,z3i
	double precision, allocatable, dimension(:,:,:) :: A1 
	double precision, allocatable, dimension(:,:) :: T,P   
  double precision, allocatable, dimension(:,:,:) :: f1
    
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    z2_bnd_pr = mxGetPr(prhs(2))
    z3_bnd_pr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    call mxCopyPtrToReal8(z2_bnd_pr,z2_bnd,2)
    call mxCopyPtrToReal8(z3_bnd_pr,z3_bnd,2)
        
    ! Points to evaluate
    mdim = 3
    m1 = mxGetScalar(prhs(4)) 
    m2 = mxGetScalar(prhs(5)) 
    m3 = mxGetScalar(prhs(6)) 
    m = [m1, m2, m3]
    !mmax = max(m(1),m(2),m(3))+1
    nodes = m(1)*m(2)*m(3)
    allocate(z1i(m(1),m(2),m(3))) 
    allocate(z2i(m(1),m(2),m(3)))
    allocate(z3i(m(1),m(2),m(3)))       
    z1i_pr = mxGetPr(prhs(7))
    z2i_pr = mxGetPr(prhs(8))
    z3i_pr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(z1i_pr,z1i,nodes)
    call mxCopyPtrToReal8(z2i_pr,z2i,nodes)
    call mxCopyPtrToReal8(z3i_pr,z3i,nodes)

    ! Least square weights
    ndim = 3
    n1 = mxGetScalar(prhs(10)) 
    n2 = mxGetScalar(prhs(11)) 
    n3 = mxGetScalar(prhs(12)) 
    nn1 = n1 + 1
    nn2 = n2 + 1
    nn3 = n3 + 1
    n = [nn1, nn2, nn3]
    coefs = n(1)*n(2)*n(3)
    allocate(A1(n(1),n(2),n(3)))
    A1_pr = mxGetPr(prhs(13))
    call mxCopyPtrToReal8(A1_pr,A1,coefs)    
    
    ! Chebyshev Polynomial Parameters
    c = mxGetScalar(prhs(14))
    allocate(T(c,c))
    allocate(P(c,c))
    T_pr = mxGetPr(prhs(15))
    P_pr = mxGetPr(prhs(16))
    call mxCopyPtrToReal8(T_pr,T,c*c)
    call mxCopyPtrToReal8(P_pr,P,c*c)

    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(m(1),m(2),m(3)))
    plhs(1) = mxCreateNumericArray(mdim,m,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allcheb(f1,nodes,z1i,z2i,z3i,z1_bnd,z2_bnd,z3_bnd,n(1),n(2),n(3),m(1),m(2),m(3),c,T,P,A1)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,nodes)
    
    ! Deallocate arrays
    deallocate(z1i)    
    deallocate(z2i)    
    deallocate(z3i)
    deallocate(A1)     
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)
    
end subroutine mexFunction

subroutine allcheb(f1,nodes,z1i,z2i,z3i,z1_bnd,z2_bnd,z3_bnd,nn1,nn2,nn3,m1,m2,m3,mmax,T,P,A1)

    implicit none
    mwSize :: nodes,nn1,nn2,nn3,m1,m2,m3,mmax
    double precision, dimension(m1,m2,m3) :: z1i,z2i,z3i,f1
    double precision, dimension(2) :: z1_bnd,z2_bnd,z3_bnd
    double precision, dimension(mmax,mmax) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1
    double precision, dimension(nn2,nn2) :: T2,P2
    double precision, dimension(nn3,nn3) :: T3,P3  
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn2) :: vec2
    double precision, dimension(nn3) :: vec3
    double precision, dimension(nn1,nn2,nn3) :: A1
    double precision :: temp1,x1i,x2i,x3i
    mwSize :: i1,i2,i3,j1,j2,j3
      
    T1 = T(1:nn1,1:nn1)
    T2 = T(1:nn2,1:nn2)
    T3 = T(1:nn3,1:nn3)   
    P1 = P(1:nn1,1:nn1)
    P2 = P(1:nn2,1:nn2)
    P3 = P(1:nn3,1:nn3)
    do i3 = 1,m3
     do i2 = 1,m2
      do i1 = 1,m1
       ! Transform Z in [a,b] to X in [-1,1]
       x1i = 2*(z1i(i1,i2,i3)-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1;
       x2i = 2*(z2i(i1,i2,i3)-z2_bnd(1))/(z2_bnd(2) - z2_bnd(1)) - 1;
       x3i = 2*(z3i(i1,i2,i3)-z3_bnd(1))/(z3_bnd(2) - z3_bnd(1)) - 1;
    
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
       f1(i1,i2,i3) = temp1
      end do
     end do
    end do
    
end subroutine allcheb