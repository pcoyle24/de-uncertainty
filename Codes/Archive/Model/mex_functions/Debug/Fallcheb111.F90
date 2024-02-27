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
		mwSize :: m,n,c
		mwSize mdim,ndim,mmax,nn1
	  mwSize nodes,coefs,chebs
		integer*4 myclassid
		double precision, dimension(2) :: z1_bnd
		double precision, allocatable, dimension(:) :: z1i	
		double precision, allocatable, dimension(:) :: A1 
		double precision, allocatable, dimension(:,:) :: T,P   
	  double precision, allocatable, dimension(:) :: f1
	  
    
    ! Load Inputs
    ! Bounds
    z1_bnd_pr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(z1_bnd_pr,z1_bnd,2)
    
    ! Points to evaluate
    mdim = 1
    m = mxGetScalar(prhs(2)) 
    nodes = m
    mmax = m + 1
    allocate(z1i(m))  
    z1i_pr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(z1i_pr,z1i,nodes)

    ! Least square weights
    ndim = 1
    n = mxGetScalar(prhs(4))
    coefs = n
    nn1 = n + 1 
    !if (nn1 .eq. 5) then
		!	call mexErrMsgTxt ( 'nn1 is 5' )
    !end if
    !call mexErrMsgTxt ( 'nn1 is not 5' )
    allocate(A1(nn1))
    A1_pr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(A1_pr,A1,nn1)    
    
    ! Chebyshev Polynomial Parameters
    allocate(T(mmax,mmax))
    allocate(P(mmax,mmax))
    T_pr = mxGetPr(prhs(6))
    P_pr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(T_pr,T,mmax*mmax)
    call mxCopyPtrToReal8(P_pr,P,mmax*mmax)
    
    ! Create array for return argument
    myclassid = mxClassIDFromClassName('double')
    allocate(f1(m))
    plhs(1) = mxCreateNumericArray(mdim,m,myclassid,0)
    f1_pr = mxGetPr(plhs(1))
    
    ! Call subroutine for assignment
    call allcheb(f1,nodes,z1i,z1_bnd,nn1,m,mmax,T,P,A1)
    
    ! Load Fortran array to pointer (output to MATLAB)
    call mxCopyReal8ToPtr(f1,f1_pr,nodes)
    
    ! Deallocate arrays
    deallocate(z1i)    
    deallocate(A1)     
    deallocate(T)  
    deallocate(P)  
    deallocate(f1)
    
end subroutine mexFunction

subroutine allcheb(f1,nodes,z1i,z1_bnd,nn1,m1,mmax,T,P,A1)

    implicit none
    mwSize :: nodes,nn1,m1,mmax
    double precision, dimension(m1) :: z1i,f1
    double precision, dimension(2) :: z1_bnd
    double precision, dimension(mmax,mmax) :: T,P
    double precision, dimension(nn1,nn1) :: T1,P1 
    double precision, dimension(nn1) :: vec1
    double precision, dimension(nn1) :: A1
    double precision :: temp1,x1i
    integer i_stat
    mwSize :: i1,j1
    
      
    T1 = T(1:nn1,1:nn1)
    P1 = P(1:nn1,1:nn1)
    do i1 = 1,m1
      ! Transform Z in [a,b] to X in [-1,1]
      
      
      
      !if (z1_bnd(2) .eq. 1.025d0) then
			!	call mexErrMsgTxt ( 'they are equal' )
    	!end if
    	!call mexErrMsgTxt ( 'they are NOT equal' )
      
      
      !if (abs(z1i(i1) - 0.975254463952977d0) .le. 1e-15) then
			!	call mexErrMsgTxt ( 'they are equal!!' )
    	!end if
    	!call mexErrMsgTxt ( 'they are NOT equal' )
      
      
      x1i = 2d0*(z1i(i1)-z1_bnd(1))/(z1_bnd(2) - z1_bnd(1)) - 1d0;
			!if (abs(x1i - (-0.989821441880933d0)) .le. 1e-15) then
			!	call mexErrMsgTxt ( 'they are equal' )
    	!end if
    	!call mexErrMsgTxt ( 'they are NOT equal' )
      ! Evaluate Chebyshev Polynomials
      vec1 = sum(T1*x1i**P1,2);
      
      !open(unit = 3, file = 'temp.dat', status = 'replace', action = 'write', iostat = i_stat)
			!300 format(f25.15,2x,f25.15,2x,f25.15,2x)
			!do j1 = 1,nn1
      !  write(3,300)  dble(j1), vec1(j1)
      !end do
      !close(unit = 3)
       
      !if (abs(vec1(3) - (0.959492973614498d0)) .le. 1e-15) then
			!	call mexErrMsgTxt ( 'they are equal' )
    	!end if
    	!call mexErrMsgTxt ( 'they are NOT equal!!!' )
      ! Evaluate function
      temp1 = 0
      !open(unit = 3, file = 'temp.dat', status = 'replace', action = 'write', iostat = i_stat)
			!300 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)
      do j1 = 1,nn1
        temp1 = temp1 + A1(j1)*vec1(j1)
       ! write(3,300)  dble(j1), vec1(j1)*A1(j1), temp1
      end do
      !close(unit = 3)
      !call mexErrMsgTxt ( 'stop' )
      f1(i1) = temp1
      
    end do
    
end subroutine allcheb