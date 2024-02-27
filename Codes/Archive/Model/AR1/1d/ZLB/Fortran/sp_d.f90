! ************************************************************************
! Filename : sp_d.f90                                                    
!                                                                         
! Author : Philip Coyle                                                   
!                                                                         
! Date Created : August 17 2018                                      
!                                                                         
! Description : This program will use a method to minimize the sum 
! of squared residuals, following a Christiano-Fisher solition method
! approach to solve a simple New Keynesian economy with a standard Taylor Rule. 
!                                                                         
! Routine:                                                                
! cd /mq/philipprojects/RA_Work/Taisuke_Nakata/Zero_Lower_Bound/DeflationaryRegime/RSSInflation/Codes/AR1/1d/ZLB/Fortran
! sh                                                                      
! . /etc/profile.d/modules.sh
! module load comp/intel/14.0.5 lib/mkl/11.1.4 lib/imsl/7.1.0
! $FC $F90FLAGS -o sp_d sp_d.f90 $LINK_FNL -heap-arrays
! ./sp_d
! ************************************************************************

! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! module : global_params                                                         
!                                                                         
! Description : This module will form the foudation for our program. In it
! we will allocate space for all paramaters used in this program.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

module global_params

use Numerical_Libraries

implicit none

! -----------------------------------------------------------------------
! *******************DECLARATION OF PARAMETERS AND VARIABLES*************
! -----------------------------------------------------------------------
! Model Parameters
double precision, parameter :: 				cBET 							= 1d0/(1.0025d0) 
double precision, parameter :: 				cCHIc 						= 1d0 
double precision, parameter :: 				cCHIn 						= 1d0 
double precision, parameter :: 				cTHETA 						= 11d0 
double precision, parameter :: 				cTAU 							= 1d0/cTHETA 
double precision, parameter :: 				cIOTA 						= 1d0
double precision, parameter :: 				cVARPHI	  				= 500d0 
double precision, parameter :: 				cPHIpi 						= 16d0 
double precision, parameter :: 				cPHIy 						= 0d0 
double precision, parameter :: 				cRHO      				= 0.8d0 
double precision, parameter :: 				cSIGMAd	  				= (0.455d0/100d0) 

double precision, parameter :: 				cALPHA 						= 1d0
double precision, parameter :: 				cPItarg 					= (2d0/400d0) + 1d0



! Tolerance level for convergence and max iterations
double precision, parameter :: 				tol			 					= 1d-10
integer, parameter 					:: 				max_iter 					= 20000
integer, parameter 					:: 				max_func_eval 		= 100000
integer 										::				it		 		= 0 !iteration counter
integer 										:: 				converged		= 0

! the value for pi (for GH Integration)
double precision, parameter :: 				pi 								= acos(-1.0d0)

! -----------------------------------------------------------------------
! *************************STEADY STATE VALUES***************************
! -----------------------------------------------------------------------  
double precision            :: 				Cbar	
double precision            :: 				PIbar 
double precision            :: 				Nbar  			
double precision            :: 				Ybar  			
double precision            :: 				Wbar  			 			
double precision            :: 				Rbar  		  
double precision,parameter  :: 				Rbar_zlb					= 1d0 		
double precision            :: 				Vbar  			  		
double precision,parameter  :: 				DELbar  					= 1d0

! -----------------------------------------------------------------------
! ****************************GRID SET**********************************
! ----------------------------------------------------------------------- 
! Set up for the GH nodes and weights
integer, parameter 				  :: 				n_e 							= 10
double precision 						:: 				e_nodes(n_e), e_nodes_ed(n_e)
double precision 						:: 				e_weights(n_e), e_weights_ed(n_e)
integer						 				  :: 				i_e, it_e, it_ed


! Set up for discritizing the state space 
integer						 				  :: 				it_del
integer, parameter 				  :: 				n_del 				= 101
integer, parameter 				  :: 				griddim 					= n_del
double precision 						:: 				del_grid(n_del)
double precision, parameter			:: 				cSIGMAded = (cSIGMAd**2d0/(1d0-cRHO**2d0))**(0.5d0)
double precision, parameter :: 				min_del 					= DELbar - 4d0*cSIGMAded
double precision, parameter :: 				max_del 					= DELbar + 4d0*cSIGMAded
double precision, parameter,dimension(2,1) :: 				del_bound 				= (/min_del, max_del/)
double precision					  :: 				del_today
double precision					  :: 				del_today_zlb
double precision					  :: 				del_tomorrow(n_e)

integer						 				  :: 				it_del_eqm

integer 										::			  i_stat


! set up for sum of nonlinear solver
double precision, parameter				:: 				pi_yesterday = 1d0

double precision 						:: 				c_today
double precision 						:: 				pi_today
double precision 						:: 				pi_tilde_today
double precision 						:: 				n_today
double precision 						:: 				y_today
double precision 						:: 				w_today
double precision 						:: 				v_today
double precision 						:: 				r_today

double precision 						:: 				c_today_zlb
double precision 						:: 				pi_today_zlb
double precision 						:: 				pi_tilde_today_zlb
double precision 						:: 				n_today_zlb
double precision 						:: 				y_today_zlb
double precision 						:: 				w_today_zlb
double precision 						:: 				v_today_zlb
double precision						:: 				r_today_zlb 

double precision 						:: 				c_tomorrow
double precision 						:: 				pi_tomorrow
double precision 						:: 				pi_tilde_tomorrow
double precision 						:: 				n_tomorrow
double precision 						:: 				y_tomorrow
double precision 						:: 				w_tomorrow
double precision 						:: 				v_tomorrow
double precision 						:: 				r_tomorrow

double precision 						:: 				c_tomorrow_zlb
double precision 						:: 				pi_tomorrow_zlb
double precision 						:: 				pi_tilde_tomorrow_zlb
double precision 						:: 				n_tomorrow_zlb
double precision 						:: 				y_tomorrow_zlb
double precision 						:: 				v_tomorrow_zlb

double precision 						:: 				exp_ee
double precision 						:: 				exp_pc
double precision 						:: 				exp_v
double precision 						:: 				exp_ee_zlb
double precision 						:: 				exp_pc_zlb
double precision 						:: 				exp_v_zlb
double precision 						:: 				exp_ee_int
double precision 						:: 				exp_pc_int
double precision 						:: 				exp_v_int
double precision 						:: 				exp_ee_int_zlb
double precision 						:: 				exp_pc_int_zlb
double precision 						:: 				exp_v_int_zlb

! Allocating space for Policy Functions

double precision 						:: 				pf_c(n_del)
double precision 						:: 				pf_pi(n_del)
double precision 						:: 				pf_n(n_del)
double precision 						:: 				pf_y(n_del)
double precision 						:: 				pf_w(n_del)
double precision 						:: 				pf_v(n_del)
double precision 						:: 				pf_r(n_del)


double precision 						:: 				pf_c_zlb(n_del)
double precision 						:: 				pf_pi_zlb(n_del)
double precision 						:: 				pf_n_zlb(n_del)
double precision 						:: 				pf_y_zlb(n_del)
double precision 						:: 				pf_w_zlb(n_del)
double precision 						:: 				pf_v_zlb(n_del)
double precision 						:: 				pf_r_zlb(n_del)

!$OMP THREADPRIVATE(it_e, it_ed,it_del, del_today, del_tomorrow, c_today, c_today_zlb, c_tomorrow, c_tomorrow_zlb, pi_today, pi_today_zlb, pi_tomorrow, pi_tomorrow_zlb, pi_tilde_today, pi_tilde_today_zlb, pi_tilde_tomorrow, pi_tilde_tomorrow_zlb, n_today, n_today_zlb, n_tomorrow, n_tomorrow_zlb, y_today, y_today_zlb, y_tomorrow, y_tomorrow_zlb, w_today, w_today_zlb, r_today, r_tomorrow, v_today, v_today_zlb, v_tomorrow, v_tomorrow_zlb, exp_ee, exp_ee_zlb, exp_pc, exp_pc_zlb, exp_v, exp_v_zlb, exp_ee_int, exp_ee_int_zlb, exp_pc_int, exp_pc_int_zlb, exp_v_int, exp_v_int_zlb)
end module global_params

                                                            
! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : StandardTR                                                         
!                                                                         
! Description : This program will use time path iteration following a
! Christiano-Fisher solition method approach to solve a simple New
! Keynesian economy with a standard Taylor Rule.                              
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

program StandardTR

use Numerical_Libraries
use global_params
use omp_lib

implicit none

! allocating space for policy function updates

double precision 						:: 				pf_c_up(n_del)
double precision 						:: 				pf_pi_up(n_del)
double precision 						:: 				pf_n_up(n_del)
double precision 						:: 				pf_y_up(n_del)
double precision 						:: 				pf_w_up(n_del)
double precision 						:: 				pf_v_up(n_del)
double precision 						:: 				pf_r_up(n_del)

double precision 						:: 				pf_c_zlb_up(n_del)
double precision 						:: 				pf_pi_zlb_up(n_del)
double precision 						:: 				pf_n_zlb_up(n_del)
double precision 						:: 				pf_y_zlb_up(n_del)
double precision 						:: 				pf_w_zlb_up(n_del)
double precision 						:: 				pf_v_zlb_up(n_del)



integer, parameter          ::        n_eq_ss = 4
double precision, parameter	:: 				tol_ss = 1d-32
integer, parameter          ::        max_iter_ss = 10000
double precision 						:: 				x_guess_ss(n_eq_ss), x_out_ss(n_eq_ss), fnorm_ss
double precision 						:: 				residuals_ss(n_eq_ss)

integer, parameter 					:: 				n_eq 				= 2 ! number of equations for nonlinear solver
integer, parameter 					:: 				n_eq_zlb 		= 2 ! number of equations for nonlinear solver
double precision					  :: 				tol_eqm 		= 1d-32
integer, parameter 					:: 				max_iter_eqm= 10000
double precision 						:: 				x_guess(n_eq), x_out(n_eq), fnorm
double precision 						:: 				x_guess_zlb(n_eq_zlb), x_out_zlb(n_eq_zlb), fnorm_zlb
double precision 						:: 				residuals(n_eq),residuals_zlb(n_eq_zlb)
integer 										:: 				it_scale


double precision 						:: 				diff_c
double precision 						:: 				diff_pi
double precision 						:: 				diff_n
double precision 						:: 				diff_y
double precision 						:: 				diff_w
double precision 						:: 				diff_v 
double precision 						:: 				diff_r


double precision 						:: 				diff_c_zlb
double precision 						:: 				diff_pi_zlb
double precision 						:: 				diff_n_zlb
double precision 						:: 				diff_y_zlb
double precision 						:: 				diff_w_zlb
double precision 						:: 				diff_v_zlb

double precision 						:: 				max_diff

! Begin Computational Timer
integer 										::				beginning, end, rate

external Gauss_Hermite
external linspace
external lin_interp
external ss_solve
external eqm_nzlb
external eqm_zlb

call system_clock(beginning, rate)

! Discretizing the state space
write (*,*) "Discretizing the State Space"
write (*,*) ""
call linspace(min_del,max_del,n_del,del_grid)
!write(*,*) del_grid
!stop

! Determine Steady State Values
x_guess_ss = (/ 1d0, 1d0, 1d0, -200d0/)

write(*,*) "Computing the Steady State Values"
call erset (0,0,0)
call dneqnf(ss_solve,tol_ss,n_eq_ss,max_iter_ss,x_guess_ss,x_out_ss,fnorm_ss)

Cbar = x_out_ss(1)
PIbar = x_out_ss(2)
Nbar = x_out_ss(3)
Ybar = Nbar
Wbar = Cbar**cCHIc*Nbar**cCHIn
Rbar = cPItarg/cBET*((PIbar/cPItarg)**cPHIpi*(Ybar/Ybar)**(cPHIy))
Vbar = x_out_ss(4)




write (*,*) "**********************************************"
write (*,*) "STEADY STATE VALUES"
write (*,*) "Cbar = ", Cbar
write (*,*) "PIbar = ", 400d0*(PIbar-1d0)
write (*,*) "Nbar = ", Nbar
write (*,*) "Ybar = ", Ybar
write (*,*) "Wbar = ", Wbar
write (*,*) "Vbar = ", Vbar
write (*,*) "Rbar = ", 400d0*(Rbar-1d0)
write (*,*) "**********************************************"
write (*,*) ""


! Get GH nodes and weights for Quadrature
call Gauss_Hermite(n_e,e_nodes,e_weights,pi)

e_nodes_ed = 2d0**(0.5d0)*cSIGMAd*e_nodes

write(*,*) "*********Gaussian-Hermite Nodes*********"
write(*,*) e_nodes_ed
write(*,*) "**********Gaussian-Hermite Weights**********"
write(*,*) e_weights
write(*,*) " "

!write(*,*) "Setting up Policy Function guesses"

do it_del = 1,n_del
	pf_c(it_del) 			= Cbar
	pf_pi(it_del) 		= PIbar
	pf_n(it_del) 			= Nbar
	pf_y(it_del) 			= Ybar
	pf_w(it_del) 			= Wbar
	pf_v(it_del) 			= Vbar
	pf_r(it_del) 			= Rbar
	
	pf_c_zlb(it_del) 			= Cbar
	pf_pi_zlb(it_del) 		= PIbar
	pf_n_zlb(it_del) 			= Nbar
	pf_y_zlb(it_del) 			= Ybar
	pf_w_zlb(it_del) 			= Wbar
	pf_v_zlb(it_del) 			= Vbar 
end do


write(*,*)' -------------------------------------------------------'
write(*,*)' -------------------------------------------------------'
write(*,*)' -------------------------------------------------------'
write(*,*)' ------------ TIME ITERATION (LESS VALUE) --------------'
write(*,*)' -------------------------------------------------------'
write(*,*)' -------------------------------------------------------'
write(*,*)' -------------------------------------------------------'

converged = 0
it = 0

! Time Path Iteration
do while (converged == 0 .and. it < max_iter)
	call omp_set_num_threads(12)
	!$OMP PARALLEL DO PRIVATE(tol_eqm, x_guess, x_out, fnorm, x_guess_zlb, x_out_zlb)
		
	do it_del = 1,n_del

		! ******************************
		! Define today's state variables
		! ******************************
		del_today = del_grid(it_del)
		
		
		! Define tomorrow's state varibale
		do it_ed = 1,n_e
			del_tomorrow(it_ed) = cRHO*(del_today - 1d0) + 1d0 + e_nodes_ed(it_ed)
		end do 


		! ******************************
		! Solve for the non-binding case
		! ******************************
		
		x_guess(1) = pf_c(it_del)
		x_guess(2) = pf_pi(it_del)

		call erset (0,0,0)
		call dneqnf(eqm_nzlb,tol_eqm,n_eq,max_iter_eqm,x_guess,x_out,fnorm)

		c_today = x_out(1)
		pi_today = x_out(2)

		! Build out other PFs
		pi_tilde_today    = pi_today/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
		n_today     = (c_today/(1d0-(cVARPHI/2d0)*(pi_tilde_today-1d0)**2d0))
		y_today     = n_today
		w_today     = n_today**cCHIn*c_today**cCHIc
		r_today     = cPItarg/cBET*((pi_today/cPItarg)**(cPHIpi)*(y_today/Ybar)**(cPHIy))

		exp_ee_int = 0d0; 
    exp_pc_int = 0d0; 

    !Errors in EE and PC are now multiplied by R_t and \delta 
    do it_e = 1,n_e
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
        if (r_tomorrow >= 1d0) then
            ! Get Interpolated Values
            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c,c_tomorrow)
            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi,pi_tomorrow)

             ! Build out other PF Values
             pi_tilde_tomorrow  = pi_tomorrow/(cPItarg**cIOTA*pi_today**(1d0-cIOTA))**cALPHA;
             n_tomorrow         = (c_tomorrow/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow-1d0)**2d0));
             y_tomorrow         = n_tomorrow;

             exp_ee = c_tomorrow**(-cCHIc)*pi_tomorrow**(-1d0);
             exp_pc = (y_tomorrow/c_tomorrow**(cCHIc))*cVARPHI*(pi_tilde_tomorrow - 1d0)*pi_tilde_tomorrow;

            exp_ee_int = exp_ee_int + pi**(-0.5d0)*e_weights(it_e)*exp_ee; 
            exp_pc_int = exp_pc_int + pi**(-0.5d0)*e_weights(it_e)*exp_pc;
        else
            ! Get Interpolated Values
            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c_zlb,c_tomorrow_zlb)
            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi_zlb,pi_tomorrow_zlb)

             ! Build out other PF Values
             pi_tilde_tomorrow_zlb  = pi_tomorrow_zlb/(cPItarg**cIOTA*pi_today**(1d0-cIOTA))**cALPHA;
             n_tomorrow_zlb         = (c_tomorrow_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow_zlb-1d0)**2d0));
             y_tomorrow_zlb         = n_tomorrow_zlb;

             exp_ee_zlb= c_tomorrow_zlb**(-cCHIc)*pi_tomorrow_zlb**(-1d0);
             exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb**(cCHIc))*cVARPHI*(pi_tilde_tomorrow_zlb - 1d0)*pi_tilde_tomorrow_zlb;

            exp_ee_int = exp_ee_int + pi**(-0.5d0)*e_weights(it_e)*exp_ee_zlb; 
            exp_pc_int = exp_pc_int + pi**(-0.5d0)*e_weights(it_e)*exp_pc_zlb;
        end if
    end do

		! ******************************
		! Solve for the binding case
		! ******************************
		
			
		x_guess_zlb(1) = pf_c_zlb(it_del)
		x_guess_zlb(2) = pf_pi_zlb(it_del)
		
		call erset (0,0,0)
		call dneqnf(eqm_zlb,tol_eqm,n_eq_zlb,max_iter_eqm,x_guess_zlb,x_out_zlb,fnorm)

		c_today_zlb = x_out_zlb(1)
		pi_today_zlb = x_out_zlb(2)

		! Build out other PFs
        pi_tilde_today_zlb    = pi_today_zlb/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
        n_today_zlb     = (c_today_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_today_zlb-1d0)**2d0))
        y_today_zlb     = n_today_zlb
        w_today_zlb     = n_today_zlb**cCHIn*c_today_zlb**cCHIc
        r_today_zlb     = 1d0

        exp_ee_int_zlb = 0d0; 
        exp_pc_int_zlb = 0d0; 

        !Errors in EE and PC are now multiplied by R_t and \delta 
        do it_e = 1,n_e
            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
            if (r_tomorrow >= 1d0) then
                ! Get Interpolated Values
                call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c,c_tomorrow)
                call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi,pi_tomorrow)

                 ! Build out other PF Values
                 pi_tilde_tomorrow  = pi_tomorrow/(cPItarg**cIOTA*pi_today_zlb**(1d0-cIOTA))**cALPHA;
                 n_tomorrow         = (c_tomorrow/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow-1d0)**2d0));
                 y_tomorrow         = n_tomorrow;

                 exp_ee = c_tomorrow**(-cCHIc)*pi_tomorrow**(-1d0);
                 exp_pc = (y_tomorrow/c_tomorrow**(cCHIc))*cVARPHI*(pi_tilde_tomorrow - 1d0)*pi_tilde_tomorrow;

                exp_ee_int_zlb = exp_ee_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_ee; 
                exp_pc_int_zlb = exp_pc_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_pc;
            else
                ! Get Interpolated Values
                call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c_zlb,c_tomorrow_zlb)
                call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi_zlb,pi_tomorrow_zlb)

                 ! Build out other PF Values
                 pi_tilde_tomorrow_zlb  = pi_tomorrow_zlb/(cPItarg**cIOTA*pi_today_zlb**(1d0-cIOTA))**cALPHA;
                 n_tomorrow_zlb         = (c_tomorrow_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow_zlb-1d0)**2d0));
                 y_tomorrow_zlb         = n_tomorrow_zlb;

                 exp_ee_zlb= c_tomorrow_zlb**(-cCHIc)*pi_tomorrow_zlb**(-1d0);
                 exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb**(cCHIc))*cVARPHI*(pi_tilde_tomorrow_zlb - 1d0)*pi_tilde_tomorrow_zlb;
                 
                exp_ee_int_zlb = exp_ee_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_ee_zlb; 
                exp_pc_int_zlb = exp_pc_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_pc_zlb;
            end if
        end do

        ! ******************************
				! Populate Updates
				! ******************************

        pf_c_up(it_del) = c_today
        pf_pi_up(it_del) = pi_today
        pf_n_up(it_del) = n_today
        pf_y_up(it_del) = y_today
        pf_w_up(it_del) = w_today
        pf_r_up(it_del) = r_today

        pf_c_zlb_up(it_del) = c_today_zlb
        pf_pi_zlb_up(it_del) = pi_today_zlb
        pf_n_zlb_up(it_del) = n_today_zlb
        pf_y_zlb_up(it_del) = y_today_zlb
        pf_w_zlb_up(it_del) = w_today_zlb
    end do	
    !$OMP END PARALLEL DO

  ! ******************************
	! Policy Function Differences
	! ******************************

	diff_c = sum(abs(pf_c_up - pf_c))
	diff_pi = sum(abs(pf_pi_up - pf_pi))
	diff_n = sum(abs(pf_n_up - pf_n))
	diff_y = sum(abs(pf_y_up - pf_y))
	diff_w = sum(abs(pf_w_up - pf_w))
	diff_r = sum(abs(pf_r_up - pf_r))

	diff_c_zlb = sum(abs(pf_c_zlb_up - pf_c_zlb))
	diff_pi_zlb = sum(abs(pf_pi_zlb_up - pf_pi_zlb))
	diff_n_zlb = sum(abs(pf_n_zlb_up - pf_n_zlb))
	diff_y_zlb = sum(abs(pf_y_zlb_up - pf_y_zlb))
	diff_w_zlb = sum(abs(pf_w_zlb_up - pf_w_zlb))

	! Determine Max Difference
	max_diff = diff_c + diff_pi + diff_n + diff_y + diff_w + diff_r + diff_c_zlb + diff_pi_zlb + diff_n_zlb + diff_y_zlb + diff_w_zlb

	it = it + 1

	if (mod(it,500) == 0) then 
		write(*,*) ""	
		write(*,*) "********************************************"	
		write(*,*) "At Iteration = ", it
		write(*,*) "Max Difference = ", max_diff
		write(*,*) "********************************************"
		write(*,*) ""
	end if  

	if (max_diff < tol) then 
		converged = 1
		write(*,*) ""	
		write(*,*) "********************************************"	
		write(*,*) "At Iteration = ", it
		write(*,*) "Max Difference = ", max_diff
		write(*,*) "********************************************"
		write(*,*) ""
	end if

	pf_c = pf_c_up
  pf_pi = pf_pi_up
  pf_n = pf_n_up
  pf_y = pf_y_up
  pf_w = pf_w_up
  pf_r = pf_r_up

  pf_c_zlb = pf_c_zlb_up
  pf_pi_zlb = pf_pi_zlb_up
  pf_n_zlb = pf_n_zlb_up
  pf_y_zlb = pf_y_zlb_up
  pf_w_zlb = pf_w_zlb_up
end do


!write(*,*)' -------------------------------------------------------'
!write(*,*)' -------------------------------------------------------'
!write(*,*)' -------------------------------------------------------'
!write(*,*)' ------------- VALUE FUNCTION ITERATION ----------------'
!write(*,*)' -------------------------------------------------------'
!write(*,*)' -------------------------------------------------------'
!write(*,*)' -------------------------------------------------------'
!
!
!converged = 0
!it = 0
!
!! Time Path Iteration
!do while (converged == 0 .and. it < max_iter)
!	call omp_set_num_threads(12)
!	!$OMP PARALLEL DO PRIVATE(tol_eqm, x_guess, x_out, fnorm, x_guess_zlb, x_out_zlb)
!		
!	do it_del = 1,n_del
!
!		! ******************************
!		! Define today's state variables
!		! ******************************
!		del_today = del_grid(it_del)
!		
!		
!		! Define tomorrow's state varibale
!		do it_ed = 1,n_e
!			del_tomorrow(it_ed) = cRHO*(del_today - 1d0) + 1d0 + e_nodes_ed(it_ed)
!		end do 
!
!
!		! ******************************
!		! Solve for the non-binding case
!		! ******************************
!		
!		c_today = pf_c(it_del)
!		pi_today = pf_pi(it_del)
!
!		! Build out other PFs
!		pi_tilde_today    = pi_today/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
!		n_today     = (c_today/(1d0-(cVARPHI/2d0)*(pi_tilde_today-1d0)**2d0))
!		y_today     = n_today
!		w_today     = n_today**cCHIn*c_today**cCHIc
!		r_today     = cPItarg/cBET*((pi_today/cPItarg)**(cPHIpi)*(y_today/Ybar)**(cPHIy))
!
!    exp_v_int  = 0d0;
!
!    ! Errors in EE and PC are now multiplied by R_t and \delta 
!    do it_e = 1,n_e
!        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
!        if (r_tomorrow >= 1d0) then
!            ! Get Interpolated Values
!            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_v,v_tomorrow)
!            exp_v  = v_tomorrow
!            exp_v_int = exp_v_int + pi**(-0.5d0)*e_weights(it_e)*exp_v
!        else
!            ! Get Interpolated Values
!            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_v_zlb,v_tomorrow_zlb)
!						exp_v_zlb  = v_tomorrow_zlb;
!						exp_v_int = exp_v_int + pi**(-0.5d0)*e_weights(it_e)*exp_v_zlb;
!        end if
!    end do
!
!    if (cCHIc == 1d0) then
!        v_today = log(c_today) - n_today**(1d0+cCHIn)/(1d0+cCHIn) + cBET*del_today*exp_v_int;
!    else
!        v_today = c_today**(1d0-cCHIc)/(1d0-cCHIc) - n_today**(1d0+cCHIn)/(1d0+cCHIn) + cBET*del_today*exp_v_int;
!    end if
!
!
!		! ******************************
!		! Solve for the binding case
!		! ******************************
!		
!		c_today_zlb = pf_c_zlb(it_del)
!		pi_today_zlb = pf_pi_zlb(it_del)
!		
!		! Build out other PFs
!    pi_tilde_today_zlb    = pi_today_zlb/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
!    n_today_zlb     = (c_today_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_today_zlb-1d0)**2d0))
!    y_today_zlb     = n_today_zlb
!    w_today_zlb     = n_today_zlb**cCHIn*c_today_zlb**cCHIc
!    r_today_zlb     = 1d0
!    
!    exp_v_int_zlb  = 0d0;
!
!    ! Errors in EE and PC are now multiplied by R_t and \delta 
!    do it_e = 1,n_e
!        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
!        if (r_tomorrow >= 1d0) then
!            ! Get Interpolated Values
!            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_v,v_tomorrow)
!						exp_v  = v_tomorrow;
!						exp_v_int_zlb = exp_v_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_v;
!        else
!            ! Get Interpolated Values
!            call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_v_zlb,v_tomorrow_zlb)
!						exp_v_zlb  = v_tomorrow_zlb;
!						exp_v_int_zlb = exp_v_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_v_zlb;
!        end if
!    end do
!
!    if (cCHIc == 1d0) then
!        v_today_zlb = log(c_today_zlb) - n_today_zlb**(1d0+cCHIn)/(1d0+cCHIn) + cBET*del_today*exp_v_int_zlb;
!    else
!        v_today_zlb = c_today_zlb**(1d0-cCHIc)/(1d0-cCHIc) - n_today_zlb**(1d0+cCHIn)/(1d0+cCHIn) + cBET*del_today*exp_v_int_zlb;
!    end if
!
!    ! ******************************
!		! Populate Updates
!		! ******************************
!    pf_v_up(it_del) = v_today
!    pf_v_zlb_up(it_del) = v_today_zlb
!  end do	
!  !$OMP END PARALLEL DO
!
!  ! ******************************
!	! Policy Function Differences
!	! ******************************
!
!	diff_v = sum(abs(pf_v_up - pf_v))
!	diff_v_zlb = sum(abs(pf_v_zlb_up - pf_v_zlb))
!
!	! Determine Max Difference
!	max_diff = diff_v + diff_v_zlb
!
!	it = it + 1
!
!	if (mod(it,1000) == 0) then 
!		write(*,*) ""	
!		write(*,*) "********************************************"	
!		write(*,*) "At Iteration = ", it
!		write(*,*) "Max Difference = ", max_diff
!		write(*,*) "********************************************"
!		write(*,*) ""
!	end if  
!
!	if (max_diff < tol) then 
!		converged = 1
!		write(*,*) ""	
!		write(*,*) "********************************************"	
!		write(*,*) "At Iteration = ", it
!		write(*,*) "Max Difference = ", max_diff
!		write(*,*) "********************************************"
!		write(*,*) ""
!	end if
!	
!  pf_v = pf_v_up
!  pf_v_zlb = pf_v_zlb_up
!end do

!write(*,*) 'Simulating the model to find Expected Value'
!call simul()

call system_clock(end)
write(*,*) "******************************************************"
write(*,*) "Total elapsed time = ", real(end - beginning) / real(rate)," seconds"
write(*,*) "******************************************************"
write(*,*) ""

write(*,*) ""
write(*,*) "**************************************"
write(*,*) "Writing PFs in the following way:"
write(*,*) "DEL, C, PI, N, Y, W, V, R"
write(*,*) "**************************************"

open(unit = 2, file = 'PFs_sp_d_nzlb.dat', status = 'replace', action = 'write', iostat = i_stat)
open(unit = 3, file = 'PFs_sp_d_zlb.dat', status = 'replace', action = 'write', iostat = i_stat)
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15)
300 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15)

do it_del = 1,n_del
     write(2,200) del_grid(it_del), pf_c(it_del), 400*(pf_pi(it_del)-1d0), pf_n(it_del), pf_y(it_del), pf_w(it_del), pf_v(it_del), 400*(pf_r(it_del)-1d0)
end do

do it_del = 1,n_del
     write(3,300) del_grid(it_del), pf_c_zlb(it_del), 400*(pf_pi_zlb(it_del)-1d0), pf_n_zlb(it_del), pf_y_zlb(it_del), pf_w_zlb(it_del), pf_v_zlb(it_del), 0d0
end do

write(*,*) 'delta =  ', del_grid(51)
write(*,*) 'pf_r =  ', 400*(pf_r(51)-1d0)

write(*,*) ""
write(*,*) "**************************************"
write(*,*) "************END OF PROGRAM************"
write(*,*) "**************************************"
write(*,*) ""
end program StandardTR



! ************************************************************************
! ************************************************************************
! ************************************************************************
! ******************************SUBROUTINES*******************************
! ************************************************************************
! ************************************************************************
! ************************************************************************

! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : Gauss_Hermite
! 
! description : determines the nodes (x) and weights (w) for the Gauss-
! Hermite quadrature of order n>1, on the interval [-INF, +INF].
! This is adapted from a MATLAB code (c) Geert Van Damme, 2010.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine Gauss_Hermite(n_e,e_nodes,e_weights,pi)

use lin_eig_gen_int

implicit none

! Read in the inputs and prep for output
integer, intent(in)               :: n_e
double precision, intent(inout)   :: e_nodes(n_e)
double precision, intent(inout)   :: e_weights(n_e)
double precision, intent(in)      :: pi

! Prep the elements for the Companion Matrix (CM). CM is such that 
! det(xI - CM) = L_n(x), with L_n the Hermite polynomial under
! consideration. CM is also constructed in a way that it is symmetrical.
double precision    :: a_vec(n_e - 1)
integer             :: i_a
real(kind(1d0))      CM(n_e,n_e)
integer             :: i_CM,j_CM	  

! Prep the matrices for determining the abscissas and weights. The 
! abscissas are the roots of the characteristic polynomial (the
! eigenvalues of CM); the weights can then be derived from the
! corresponding eigenvectors
complex(kind(1d0))   V(n_e,n_e) ! To be populated with eigenvectors
complex(kind(1d0))   E(n_e) ! To be populated with the eigenvalues
integer :: i_e

! Used for sorting the eigenvalues in ascending order
!integer :: flag = 1
integer :: flag
complex(kind(1d0)) swap
complex(kind(1d0)) swap_vec(n_e)
integer :: i_swap

! Populate the a, and CM vectors/matrices accordingly
do i_a = 1,n_e - 1
    a_vec(i_a) = (dble(i_a)/2d0)**(0.5d0)
end do

do i_CM = 1,n_e
	do j_CM = 1,n_e
	   CM(i_CM,j_CM) = 0d0
	end do
end do

do i_CM = 1,n_e - 1
    CM(i_CM,i_CM+1) = a_vec(i_CM)
    CM(i_CM+1,i_CM) = a_vec(i_CM)
    !write(*,*) "CM(i_CM,i_CM+1)",CM(i_CM,i_CM+1)
end do

! Solve for the eigenvalues and eigenvectors
call lin_eig_gen(CM,E,v = V)

flag=1
! Sort the eigenvalues in ascending order; these are the nodes :)
do while (flag == 1)
    flag = 0
    do i_e = 1,n_e-1
        if (real(E(i_e))>real(E(i_e+1))) then 
            flag = 1

            swap     = E(i_e + 1)
            do i_swap = 1,n_e
                swap_vec(i_swap) = V(i_swap,i_e + 1)
            end do

            E(i_e + 1)  = E(i_e)
            E(i_e)      = swap

            V(:,i_e+1) = V(:,i_e)
            do i_swap = 1,n_e
                V(i_swap,i_e)   = swap_vec(i_swap)
            end do
        end if
    end do
end do

e_nodes = dble(E)

do i_e = 1,n_e
    e_weights(i_e) = (pi**(0.5d0) * dble(real(V(1,i_e)))**2d0)
end do

return

end subroutine Gauss_Hermite

! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : linspace
! 
! description : outputs grids with equally spaced points.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine linspace(a,b,n,y)

use global_params
use Numerical_Libraries

implicit none

double precision, intent(in)			:: a
double precision, intent(in)			:: b
integer, intent(in)								:: n
double precision, intent(out)			:: y(n)

double precision 									:: step
integer 													:: it_linspace
double precision 									:: x_grid(n)

step = (b - a)/(dble(griddim) - 1d0)

do it_linspace = 1,n
	y(it_linspace) = a + (dble(it_linspace) - 1d0)*step
end do

return

end subroutine linspace

! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : lin_interp
! 
! description : employs the methods of linear interpolation. Originally
! adapted from MATLAB code by (c) Richter, Throckmorton, and Walker, 2013
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine lin_interp(x1,x1i,nx1,pf,o1)

use global_params
use Numerical_Libraries

implicit none

integer, intent(in)             :: nx1
double precision, intent(in)    :: x1(nx1)
double precision, intent(in)    :: x1i
double precision, intent(in)    :: pf(nx1)
double precision, intent(out)   :: o1

double precision :: s1
double precision :: x1i_min
double precision :: loc1

double precision :: xi
double precision :: xi_left
double precision :: xi_right

double precision :: w_2
double precision :: w_1
double precision :: w1(2)

integer :: m1

o1 = 0

s1      = x1(2) - x1(1)
x1i_min = x1i - x1(1)
loc1    = min(nx1 - 1, max(1,floor(x1i_min/s1)+1))

xi = x1i
xi_left = x1(loc1)
xi_right = x1(loc1 + 1)

w_2 = (xi - xi_left)/(xi_right - xi_left)
w_1 = 1 - w_2

w1(1) = w_1
w1(2) = w_2

do m1 = 0,1
    o1 = o1 + w1(m1+1)*pf(loc1+m1)
end do

return

end subroutine lin_interp


! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : ss_solve
! 
! description : This subroutine solves the nonlinear system of equations
! at the steady state and returns the residuals of th nonlinear system.
! Output is the steady state values for various policy functions.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine ss_solve(x_guess_ss, residuals_ss, n_eq_ss)

use global_params
use Numerical_Libraries

implicit none

integer, intent(in) 					:: n_eq_ss
double precision, intent(in)  :: x_guess_ss(n_eq_ss)
double precision, intent(out) :: residuals_ss(n_eq_ss)

double precision 			  :: PItilde


Cbar = x_guess_ss(1)
PIbar = x_guess_ss(2)
Nbar = x_guess_ss(3)
Vbar = x_guess_ss(4)


PItilde = PIbar/(cPItarg)**cALPHA
Ybar = Nbar
Wbar = Cbar**cCHIc*Nbar**cCHIn
Rbar = cPItarg/cBET*((PIbar/cPItarg)**cPHIpi*(Ybar/Ybar)**(cPHIy))


residuals_ss(1) = Cbar**(-cCHIc) - cBET*DELbar*Rbar*Cbar**(-cCHIc)*PIbar**(-1d0)
residuals_ss(2) = cVARPHI*(PItilde - 1d0)*PItilde - (1d0 - cTHETA) - cTHETA*(1d0-cTAU)*Wbar - cBET*DELbar*cVARPHI*(PItilde - 1d0)*PItilde
residuals_ss(3) = Nbar - Cbar - cVARPHI/2d0*(PItilde - 1d0)**2d0*Nbar
if (cCHIc == 1d0) then
    residuals_ss(4) = Vbar - log(Cbar) + Nbar**(1d0+cCHIn)/(1d0+cCHIn) - cBET*DELbar*Vbar;
else
    residuals_ss(4) = Vbar - Cbar**(1d0-cCHIc)/(1d0-cCHIc) + Nbar**(1d0+cCHIn)/(1d0+cCHIn) - cBET*DELbar*Vbar;
end if

return

end subroutine ss_solve


! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : eqm_nzlb
! 
! description : This subroutine solves the nonlinear system of equations
! and returns the residuals of the expctational equations (usually the 
! consumption euler equation and the phillips curve). Unconstrained Case.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine eqm_nzlb(x_guess,residuals,n_eq)

use global_params
use Numerical_Libraries

implicit none

integer, intent(in) 					:: n_eq
double precision, intent(in)  :: x_guess(n_eq)
double precision, intent(out) :: residuals(n_eq)



external chebweights11
external allcheb111

! Solve for or extract variables from inputs
c_today = x_guess(1)
pi_today = x_guess(2)

! Build out other PFs
pi_tilde_today    = pi_today/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
n_today     = (c_today/(1d0-(cVARPHI/2d0)*(pi_tilde_today-1d0)**2d0))
y_today     = n_today
w_today     = n_today**cCHIn*c_today**cCHIc
r_today     = cPItarg/cBET*((pi_today/cPItarg)**(cPHIpi)*(y_today/Ybar)**(cPHIy))

exp_ee_int = 0d0; 
exp_pc_int = 0d0; 

do it_e = 1,n_e
    call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
    if (r_tomorrow >= 1d0) then
        ! Get Interpolated Values
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c,c_tomorrow)
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi,pi_tomorrow)

         ! Build out other PF Values
         pi_tilde_tomorrow  = pi_tomorrow/(cPItarg**cIOTA*pi_today**(1d0-cIOTA))**cALPHA;
         n_tomorrow         = (c_tomorrow/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow-1d0)**2d0));
         y_tomorrow         = n_tomorrow;

         exp_ee = c_tomorrow**(-cCHIc)*pi_tomorrow**(-1d0);
         exp_pc = (y_tomorrow/c_tomorrow**(cCHIc))*cVARPHI*(pi_tilde_tomorrow - 1d0)*pi_tilde_tomorrow;

        exp_ee_int = exp_ee_int + pi**(-0.5d0)*e_weights(it_e)*exp_ee; 
        exp_pc_int = exp_pc_int + pi**(-0.5d0)*e_weights(it_e)*exp_pc;
    else
        ! Get Interpolated Values
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c_zlb,c_tomorrow_zlb)
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi_zlb,pi_tomorrow_zlb)

         ! Build out other PF Values
         pi_tilde_tomorrow_zlb  = pi_tomorrow_zlb/(cPItarg**cIOTA*pi_today**(1d0-cIOTA))**cALPHA;
         n_tomorrow_zlb         = (c_tomorrow_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow_zlb-1d0)**2d0));
         y_tomorrow_zlb         = n_tomorrow_zlb;

         exp_ee_zlb= c_tomorrow_zlb**(-cCHIc)*pi_tomorrow_zlb**(-1d0);
         exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb**(cCHIc))*cVARPHI*(pi_tilde_tomorrow_zlb - 1d0)*pi_tilde_tomorrow_zlb;

        exp_ee_int = exp_ee_int + pi**(-0.5d0)*e_weights(it_e)*exp_ee_zlb; 
        exp_pc_int = exp_pc_int + pi**(-0.5d0)*e_weights(it_e)*exp_pc_zlb;
    end if
end do

residuals(1) = c_today**(-cCHIc) - cBET*del_today*r_today*exp_ee_int;
residuals(2) = (y_today/c_today**(cCHIc))*(cVARPHI*(pi_tilde_today-1)*pi_tilde_today-(1-cTHETA)-cTHETA*(1-cTAU)*w_today) - cBET*del_today*exp_pc_int;


return

end subroutine eqm_nzlb

! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : eqm_zlb
! 
! description : This subroutine solves the nonlinear system of equations
! and returns the residuals of the expctational equations (usually the 
! consumption euler equation and the phillips curve). Constrained Case.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine eqm_zlb(x_guess_zlb,residuals_zlb,n_eq_zlb)

use global_params
use Numerical_Libraries

implicit none

integer, intent(in) 					:: n_eq_zlb
double precision, intent(in)  :: x_guess_zlb(n_eq_zlb)
double precision, intent(out) :: residuals_zlb(n_eq_zlb)


external chebweights11
external allcheb111


! Solve for or extract variables from inputs
c_today_zlb = x_guess_zlb(1)
pi_today_zlb = x_guess_zlb(2)

! Build out other PFs
pi_tilde_today_zlb    = pi_today_zlb/(cPItarg**cIOTA*pi_yesterday**(1d0-cIOTA))**cALPHA
n_today_zlb     = (c_today_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_today_zlb-1d0)**2d0))
y_today_zlb     = n_today_zlb
w_today_zlb     = n_today_zlb**cCHIn*c_today_zlb**cCHIc
r_today_zlb     = 1d0

exp_ee_int_zlb = 0d0; 
exp_pc_int_zlb = 0d0; 

do it_e = 1,n_e
    call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_r,r_tomorrow)
    if (r_tomorrow >= 1d0) then
        ! Get Interpolated Values
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c,c_tomorrow)
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi,pi_tomorrow)

         ! Build out other PF Values
         pi_tilde_tomorrow  = pi_tomorrow/(cPItarg**cIOTA*pi_today_zlb**(1d0-cIOTA))**cALPHA;
         n_tomorrow         = (c_tomorrow/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow-1d0)**2d0));
         y_tomorrow         = n_tomorrow;

         exp_ee = c_tomorrow**(-cCHIc)*pi_tomorrow**(-1d0);
         exp_pc = (y_tomorrow/c_tomorrow**(cCHIc))*cVARPHI*(pi_tilde_tomorrow - 1d0)*pi_tilde_tomorrow;

        exp_ee_int_zlb = exp_ee_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_ee; 
        exp_pc_int_zlb = exp_pc_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_pc;
    else
        ! Get Interpolated Values
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_c_zlb,c_tomorrow_zlb)
        call lin_interp(del_grid,del_tomorrow(it_e),n_del, pf_pi_zlb,pi_tomorrow_zlb)

         ! Build out other PF Values
         pi_tilde_tomorrow_zlb  = pi_tomorrow_zlb/(cPItarg**cIOTA*pi_today_zlb**(1d0-cIOTA))**cALPHA;
         n_tomorrow_zlb         = (c_tomorrow_zlb/(1d0-(cVARPHI/2d0)*(pi_tilde_tomorrow_zlb-1d0)**2d0));
         y_tomorrow_zlb         = n_tomorrow_zlb;

         exp_ee_zlb= c_tomorrow_zlb**(-cCHIc)*pi_tomorrow_zlb**(-1d0);
         exp_pc_zlb = (y_tomorrow_zlb/c_tomorrow_zlb**(cCHIc))*cVARPHI*(pi_tilde_tomorrow_zlb - 1d0)*pi_tilde_tomorrow_zlb;

        exp_ee_int_zlb = exp_ee_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_ee_zlb; 
        exp_pc_int_zlb = exp_pc_int_zlb + pi**(-0.5d0)*e_weights(it_e)*exp_pc_zlb;
    end if
end do

residuals_zlb(1) = c_today_zlb**(-cCHIc) - cBET*del_today*r_today_zlb*exp_ee_int_zlb;
residuals_zlb(2) = (y_today_zlb/c_today_zlb**(cCHIc))*(cVARPHI*(pi_tilde_today_zlb-1)*pi_tilde_today_zlb-(1-cTHETA)-cTHETA*(1-cTAU)*w_today_zlb) - cBET*del_today*exp_pc_int_zlb;


return

end subroutine eqm_zlb


! ************************************************************************

! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : simul
! 
! description : This subroutine simulates the model to find the expected
! value. 
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine simul

use global_params
use Numerical_Libraries

implicit none

integer, parameter 		:: n_sim = 1001000
integer, parameter 		:: n_burn = 1000
integer					:: it_sim


double precision :: del_yesterday, ed_today
double precision :: ud1, ud2
double precision :: zlb_check
double precision :: sum_v 
double precision :: E_val



del_yesterday = 1d0
sum_v = 0d0

do it_sim = 1,n_sim
	call random_number(ud1)
	call random_number(ud2)

	ed_today = cSIGMAd*(-2d0*log(ud1))**(1d0/2d0)*cos(2d0*pi*ud2)
    del_today = cRHO*(del_yesterday - 1d0) + 1d0 + ed_today;

    if (it_sim .GE. n_burn) then
    	call lin_interp(del_grid,del_today,n_del, pf_r,zlb_check)
    	if (zlb_check .GE. 1d0) then
    		call lin_interp(del_grid,del_today,n_del, pf_v,v_today)
    		sum_v = sum_v + v_today
    	else
    		call lin_interp(del_grid,del_today,n_del, pf_v_zlb,v_today)
    		sum_v = sum_v + v_today
    	end if
    end if 

    del_yesterday = del_today
end do

E_val = sum_v / dble(n_sim - n_burn)

!open(unit = 4, file = 'exp_val.dat', status = 'replace', action = 'write', iostat = i_stat)
!400 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15)
!write(4,400) cALPHA, 400*(cPItarg-1d0), E_val, converged
write(*,*) cSIGMAd, cALPHA, 400*(cPItarg-1d0), E_val, converged

return

end subroutine

