## Codes to run
The main code files are:
1. max_shock.m (run this code first!)
  - numerically finds the max level of uncertainty consistent with equilibrium existence for different low, medium, high. values of \phi_{\pi} and outputs cSIGMAd_max.
2. RAFRFramework.m
  - reproduces figure 1.
2. RAFR_cPHIpi_l.m
  - solves for various model RAFR for 3 different levels of uncertainty for the low level of \phi_{\pi}.
3. RAFR_cPHIpi_l.m
  - solves for various model RAFR for 3 different levels of uncertainty for the medium level of \phi_{\pi}.
4. RAFR_cPHIpi_l.m
  - solves for various model RAFR for 3 different levels of uncertainty for the high level of \phi_{\pi}.
5. RAFR_cPHIpi_m_cPItarg.m
  - solves for various model RAFR for 3 different levels of uncertainty for the medium level of \phi_{\pi} with a positive inflation target. 

The folder contains the following functions used in the above code files:
  - rouwenhorst.m: discritize state space based on Rouwenhorst method.
  - eqmmat_x: solve model system of linear equations to obtain x in {TR Target Regime), DR Def. Regime)} allocations assuming ZLB never or ZLB always binds respectively for TR and DR.
  - eqmrefine_x: solve model system of linear equations to obtain x in {TR Target Regime), DR Def. Regime)} allocations at each grid point if the ZLB does or does not bind.  
