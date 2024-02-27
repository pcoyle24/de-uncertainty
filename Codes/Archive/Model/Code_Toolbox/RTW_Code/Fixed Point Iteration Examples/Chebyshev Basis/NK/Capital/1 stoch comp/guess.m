function pf = guess(O,P,S,G)

% pf = guess(O,P,S,G) 
%   Sets the initial policy functions
% Inputs:
%   O : structure of options
%   P : structure of parameters
%   S : structure of steady state values 
%   G : structure of grids
% Outputs:
%  pf : structure of policy functions

if strcmp(O.loadpf,'guess')
    % Retrieve log-linear solution
    V = variables;
    [T,~,eu] = linmodel(P,S,V);
    disp(['eu:    ' mat2str(eu)])

    % Transform discretized state space to percent deviation from steady state
    b_gr_per = (G.l.b_gr - S.b) ./ S.b;
    r_gr_per = (G.l.r_gr - S.r) ./ S.r;
    m_gr_per = (G.l.m_gr - S.m) ./ S.m;    
    k_gr_per = (G.l.k_gr - S.k) ./ S.k;

    % Calculate linear policy functions on discretized state space
    linpf_n = zeros(G.l.griddim);
    linpf_pi = zeros(G.l.griddim);
    linpf_k = zeros(G.l.griddim);
    state = [b_gr_per(:) r_gr_per(:) m_gr_per(:) k_gr_per(:)]';
    linpf_n(:) = T(V.n,[V.b V.r V.m V.k])*state;
    linpf_pi(:) = T(V.pi,[V.b V.r V.m V.k])*state;
    linpf_k(:) = T(V.k,[V.b V.r V.m V.k])*state;     
    
    % Map policy functions to nonlinear state space
    linpf_n = squeeze(mean(linpf_n,2));
    linpf_pi = squeeze(mean(linpf_pi,2));
    linpf_k = squeeze(mean(linpf_k,2));

    % Convert back to levels    
    pf_n = P.n * (1 + linpf_n);
    pf_pi = P.pi * (1 + linpf_pi);
    pf_k = S.k * (1 + linpf_k);
    
    % Replicate over fiscal policy shock    
    pf.n = pf_n(:,:,:,ones(O.fp_pts,1));  
    pf.pi = pf_pi(:,:,:,ones(O.fp_pts,1));    
    pf.k = pf_k(:,:,:,ones(O.fp_pts,1)); 
elseif strcmp(O.loadpf,'current')    
    load('pf.mat')
end