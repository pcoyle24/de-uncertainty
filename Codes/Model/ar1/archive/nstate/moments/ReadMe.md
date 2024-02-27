## Codes to run
The main code files are:
1.  max_shock.m (run this code first!)
  - numerically finds the max level of uncertainty consistent with equilibrium existence for different low, medium, high.  values of \phi_{\pi} and outputs cSIGMAd_max.
2.  stylized_rouwenhorst_moments_cPHIpi_l.m
  - solves for various model moments along a grid of different levels of uncertainty for the low level of \phi_{\pi}. This code reproduces figure 16.
3.  stylized_rouwenhorst_moments_cPHIpi_m.m
  - solves for various model moments along a grid of different levels of uncertainty for the medium level of \phi_{\pi}. This code reproduces figure 15.
4.  stylized_rouwenhorst_moments_cPHIpi_h.m
  - solves for various model moments along a grid of different levels of uncertainty for the high level of \phi_{\pi}. This code reproduces figure 17.

The folder contains the following functions used in the above code files:
  - rouwenhorst.m: discritize state space based on Rouwenhorst method.
  - eqmmat_x: solve model system of linear equations to obtain x in {TR Target Regime. , DR Def. Regime. } allocations assuming ZLB never or ZLB always binds respecively for TR and DR.
  - eqmrefine_x: solve model system of linear equations to obtain x in {TR Target Regime. , DR Def. Regime. } allocations at each grid point if the ZLB does or respectively not bind.  
  - exp_val.m: Simulate expected value and ELB probability of economy.

All the main code files solve for the policy functions of the model. The numerical solution method as described in Appendix B1.
