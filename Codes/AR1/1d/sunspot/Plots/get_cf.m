function [pf,zlb] = get_cf(P,O,S,G,C)

global pi_yesterday

O.del_fine = 1001;
G.delfine = chebspace(O.delbound(1),O.delbound(2),O.del_fine);
G.gridfine = size(ndgrid(G.delfine));

pf_c_s_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac_s,C.max,C.T,C.P);  
pf_inf_s_fine       = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf_s,C.max,C.T,C.P);  
pf_pitilde_s_fine   = pf_inf_s_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf_n_s_fine         = (pf_c_s_fine./(1-(P.varphi/2).*(pf_pitilde_s_fine-1).^2));
pf_y_s_fine         = pf_n_s_fine;
pf_w_s_fine         = pf_n_s_fine.^P.chin.*pf_c_s_fine.^P.chic;
pf_r_s_fine         = S.pi_targ/P.beta*((pf_inf_s_fine./S.pi_targ).^(P.phi_pi).*(pf_y_s_fine./S.y).^(P.phi_y));
pf_v_s_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av_s,C.max,C.T,C.P);  

pf_c_zlb_s_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac_zlb_s,C.max,C.T,C.P);  
pf_inf_zlb_s_fine       = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf_zlb_s,C.max,C.T,C.P);  
pf_pitilde_zlb_s_fine   = pf_inf_zlb_s_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf_n_zlb_s_fine         = (pf_c_zlb_s_fine./(1-(P.varphi/2).*(pf_pitilde_zlb_s_fine-1).^2));
pf_y_zlb_s_fine         = pf_n_zlb_s_fine;
pf_w_zlb_s_fine         = pf_n_zlb_s_fine.^P.chin.*pf_c_zlb_s_fine.^P.chic;
pf_r_zlb_s_fine         = ones(length(pf_r_s_fine),1);
pf_v_zlb_s_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av_zlb_s,C.max,C.T,C.P);  

zlb.s = find(pf_r_s_fine < 1,1); %cell where ZLB binds for increasing \delta levels

pf_c_d_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac_d,C.max,C.T,C.P);  
pf_inf_d_fine       = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf_d,C.max,C.T,C.P);  
pf_pitilde_d_fine   = pf_inf_d_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf_n_d_fine         = (pf_c_d_fine./(1-(P.varphi/2).*(pf_pitilde_d_fine-1).^2));
pf_y_d_fine         = pf_n_d_fine;
pf_w_d_fine         = pf_n_d_fine.^P.chin.*pf_c_d_fine.^P.chic;
pf_r_d_fine         = S.pi_targ/P.beta*((pf_inf_d_fine./S.pi_targ).^(P.phi_pi).*(pf_y_d_fine./S.y).^(P.phi_y)); 
pf_v_d_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av_s,C.max,C.T,C.P);  

zlb.d = find(pf_r_d_fine < 1,1); %cell where ZLB binds for increasing \delta levels

pf_c_zlb_d_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ac_zlb_d,C.max,C.T,C.P);  
pf_inf_zlb_d_fine       = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Ainf_zlb_d,C.max,C.T,C.P);  
pf_pitilde_zlb_d_fine   = pf_inf_zlb_d_fine/(S.pi_targ^P.iota*pi_yesterday^(1-P.iota))^P.alpha;
pf_n_zlb_d_fine         = (pf_c_zlb_d_fine./(1-(P.varphi/2).*(pf_pitilde_zlb_d_fine-1).^2));
pf_y_zlb_d_fine         = pf_n_zlb_d_fine;
pf_w_zlb_d_fine         = pf_n_zlb_d_fine.^P.chin.*pf_c_zlb_d_fine.^P.chic;
pf_r_zlb_d_fine         = ones(length(pf_r_d_fine),1);
pf_v_zlb_d_fine         = Fallcheb111(O.delbound,O.del_fine,G.delfine,O.n1,C.Av_zlb_d,C.max,C.T,C.P);  

pf.cf.c_s       = zeros(G.gridfine);
pf.cf.inf_s     = zeros(G.gridfine);
pf.cf.pitilde_s = zeros(G.gridfine);
pf.cf.n_s       = zeros(G.gridfine);
pf.cf.y_s       = zeros(G.gridfine);
pf.cf.w_s       = zeros(G.gridfine);
pf.cf.r_s       = zeros(G.gridfine);
pf.cf.v_s       = zeros(G.gridfine);

pf.cf.c_d       = zeros(G.gridfine);
pf.cf.inf_d     = zeros(G.gridfine);
pf.cf.pitilde_d = zeros(G.gridfine);
pf.cf.n_d       = zeros(G.gridfine);
pf.cf.y_d       = zeros(G.gridfine);
pf.cf.w_d       = zeros(G.gridfine);
pf.cf.r_d       = zeros(G.gridfine);
pf.cf.v_d       = zeros(G.gridfine);

for i = 1:O.del_fine
    if pf_r_s_fine(i) >= 1 
        pf.cf.c_s(i) = pf_c_s_fine(i);
        pf.cf.inf_s(i) = pf_inf_s_fine(i);
        pf.cf.pitilde_s(i) = pf_pitilde_s_fine(i);
        pf.cf.n_s(i) = pf_n_s_fine(i);
        pf.cf.y_s(i) = pf_y_s_fine(i);
        pf.cf.w_s(i) = pf_w_s_fine(i);
        pf.cf.r_s(i) = pf_r_s_fine(i);
        pf.cf.v_s(i) = pf_v_s_fine(i);
    else
        pf.cf.c_s(i) = pf_c_zlb_s_fine(i);
        pf.cf.inf_s(i) = pf_inf_zlb_s_fine(i);
        pf.cf.pitilde_s(i) = pf_pitilde_zlb_s_fine(i);
        pf.cf.n_s(i) = pf_n_zlb_s_fine(i);
        pf.cf.y_s(i) = pf_y_zlb_s_fine(i);
        pf.cf.w_s(i) = pf_w_zlb_s_fine(i);
        pf.cf.r_s(i) = pf_r_zlb_s_fine(i);
        pf.cf.v_s(i) = pf_v_zlb_s_fine(i);
    end
    if pf_r_d_fine(i) >= 1 
        pf.cf.c_d(i) = pf_c_d_fine(i);
        pf.cf.inf_d(i) = pf_inf_d_fine(i);
        pf.cf.pitilde_d(i) = pf_pitilde_d_fine(i);
        pf.cf.n_d(i) = pf_n_d_fine(i);
        pf.cf.y_d(i) = pf_y_d_fine(i);
        pf.cf.w_d(i) = pf_w_d_fine(i);
        pf.cf.r_d(i) = pf_r_d_fine(i);
        pf.cf.v_d(i) = pf_v_d_fine(i);
    else
        pf.cf.c_d(i) = pf_c_zlb_d_fine(i);
        pf.cf.inf_d(i) = pf_inf_zlb_d_fine(i);
        pf.cf.pitilde_d(i) = pf_pitilde_zlb_d_fine(i);
        pf.cf.n_d(i) = pf_n_zlb_d_fine(i);
        pf.cf.y_d(i) = pf_y_zlb_d_fine(i);
        pf.cf.w_d(i) = pf_w_zlb_d_fine(i);
        pf.cf.r_d(i) = pf_r_zlb_d_fine(i);
        pf.cf.v_d(i) = pf_v_zlb_d_fine(i);
    end
end 
