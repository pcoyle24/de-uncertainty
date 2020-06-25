 ReadMe Documentation for "Deflationary Equilibrium Under Uncertainty" by Coyle, Nakata, and Schmidt. 

Codes to replecate the figures are found in \Draft\Figs. There are three subfolders:
1. ar1 - contains code files related to the model with an AR(1) demand shock specification
		- contains subfolder for codes that use an n-state Rouwenhorst discritization of demand shock.
2. iid - contains code files related to the model with an I.I.D. demand shock specification
	- contains subfolder for codes that use an n-state I.I.D. discritization of demand shock.
	- contains subfolder for codes that use a 3-state I.I.D. discritization of demand shock. 
3. Final - contains final figures that are published in the paper. All codes that produce figures are by default saved in this folder. 



******************************************
************* \iid\3state ****************
******************************************

********************
****** \RAFR *******
********************

All Matlab files in this folder are self contained. They all produce a different version of the Risk Adjusted Fisher Relation (RAFR) for what is found in the paper.
These codes solve for the RAFR analytically, as described in Appendix A. 

The codes are:  
--RAFR_cPHI_x.m and RAFR_cPHI_x_pitarg.m reproduce figrues 2-4 for x in {l (low),m (medium),h (high)}.
--RAFR_cPHI_lb.m and RAFR_cPHI_lb.m reproduce figure 11, which detail the infinitely many equilibria. 


***********************
****** \Moments *******
***********************

The main code files are: 
(1) max_shock.m (run this code first!)
--numerically finds the max level of uncertainty consistent with equilibrum existance for different (low, medium, high) values of \phi_{\pi} and outputs cSIGMAd_max.
(2) stylized_iid_moments_cPHIpi_l.m
--solves for various model moments along a grid of different levels of uncertainty for the low level of \phi_{\pi}. This code reproduces figure 9.
(3) stylized_iid_moments_cPHIpi_m.m
--solves for various model moments along a grid of different levels of uncertainty for the medium level of \phi_{\pi}. This code reproduces figure 8.
(4) stylized_iid_moments_cPHIpi_h.m
--solves for various model moments along a grid of different levels of uncertainty for the high level of \phi_{\pi}. This code reproduces figure 10.

The folder contains the following functions used in the above code files: 
-- iid.m: discritize state space based on IID shock.
-- eqmmat_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations assuming ZLB never or ZLB always binds respecively for TR and DR. 
-- eqmrefine_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations at each grid point if the ZLB does or does not bind.  
-- exp_val.m: Simulate expected value and ELB probablity of economy. 

All the main code files solve for the policy functions of the model. The setup is identical to the numerical solution method as decribed in Appendix B1, however, just looking at a 3state IID shock process. 


******************************************
************* \ar1\nstate ****************
******************************************

********************
****** \RAFR *******
********************

The main code files are: 
(1) max_shock.m (run this code first!)
--numerically finds the max level of uncertainty consistent with equilibrum existance for different (low, medium, high) values of \phi_{\pi} and outputs cSIGMAd_max.
(2) stylized_rouwenhorst_RAFR.m 
-- reproduces figure 1. 
(2) stylized_rouwenhorst_RAFR_cPHIpi_l.m
--solves for various model RAFR for 3 different levels of uncertainty for the low level of \phi_{\pi}. This code reproduces figure 6.
(3) stylized_rouwenhorst_RAFR_cPHIpi_m.m
--solves for various model RAFR for 3 different levels of uncertainty for the medium level of \phi_{\pi}. This code reproduces figure 5.
(4) stylized_rouwenhorst_RAFR_cPHIpi_h.m
--solves for various model RAFR for 3 different levels of uncertainty for the high level of \phi_{\pi}. This code reproduces figure 7.

The folder contains the following functions used in the above code files: 
-- rouwenhorst.m: discritize state space based on Rouwenhorst method.
-- eqmmat_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations assuming ZLB never or ZLB always binds respecively for TR and DR. 
-- eqmrefine_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations at each grid point if the ZLB does or does not bind.  

All the main code files solve for the policy functions of the model. The numerical solution method as decribed in Appendix B2.

*******************
****** \PFs *******
*******************

The main code files are: 
(1) max_shock.m (run this code first!)
--numerically finds the max level of uncertainty consistent with equilibrum existance for different (low, medium, high) values of \phi_{\pi} and outputs cSIGMAd_max.
(2) stylized_rouwenhorst_PFs_cPHIpi_l.m
--solves for various model PFs for 3 different levels of uncertainty for the low level of \phi_{\pi}. This code reproduces figure 13.
(3) stylized_rouwenhorst_PFs_cPHIpi_m.m
--solves for various model PFs for 3 different levels of uncertainty for the medium level of \phi_{\pi}. This code reproduces figure 12.
(4) stylized_rouwenhorst_PFs_cPHIpi_h.m
--solves for various model PFs for 3 different levels of uncertainty for the high level of \phi_{\pi}. This code reproduces figure 14.

The folder contains the following functions used in the above code files: 
-- rouwenhorst.m: discritize state space based on Rouwenhorst method.
-- eqmmat_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations assuming ZLB never or ZLB always binds respecively for TR and DR. 
-- eqmrefine_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations at each grid point if the ZLB does or does not bind.  

All the main code files solve for the policy functions of the model. The numerical solution method as decribed in Appendix B1.


***********************
****** \Moments *******
***********************

The main code files are: 
(1) max_shock.m (run this code first!)
--numerically finds the max level of uncertainty consistent with equilibrum existance for different (low, medium, high) values of \phi_{\pi} and outputs cSIGMAd_max.
(2) stylized_rouwenhorst_moments_cPHIpi_l.m
--solves for various model moments along a grid of different levels of uncertainty for the low level of \phi_{\pi}. This code reproduces figure 16.
(3) stylized_rouwenhorst_moments_cPHIpi_m.m
--solves for various model moments along a grid of different levels of uncertainty for the medium level of \phi_{\pi}. This code reproduces figure 15.
(4) stylized_rouwenhorst_moments_cPHIpi_h.m
--solves for various model moments along a grid of different levels of uncertainty for the high level of \phi_{\pi}. This code reproduces figure 17.

The folder contains the following functions used in the above code files: 
-- rouwenhorst.m: discritize state space based on Rouwenhorst method.
-- eqmmat_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations assuming ZLB never or ZLB always binds respecively for TR and DR. 
-- eqmrefine_x: solve model system of linear equations to obtain x in {TR (Target Regime), DR (Def. Regime)} allocations at each grid point if the ZLB does or does not bind.  
-- exp_val.m: Simulate expected value and ELB probablity of economy. 

All the main code files solve for the policy functions of the model. The numerical solution method as decribed in Appendix B1.
