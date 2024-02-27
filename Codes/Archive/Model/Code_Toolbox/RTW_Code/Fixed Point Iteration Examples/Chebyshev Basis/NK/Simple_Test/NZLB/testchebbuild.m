function out = testchebbuild(delbound,Ac,T,P,e_pts,del_gr)

% del_today = del_gr(:,ones(1,1,e_pts));
del_today = del_gr(:);

[pf_c,basis] = allcheb111_F(delbound,del_today,Ac,T,P);

pf_c_alt = (Ac'*basis)';

pf_c - pf_c_alt

out = pf_c(:,1);
