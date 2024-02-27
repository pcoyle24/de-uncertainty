function diff = check_regime(params,numpts,y,pi,i,R,regime)

% diff = check_regime(params,regime)
% Construct Matrix and vector used to solve for solutions of model
%   Inputs:
%       params: a vector of parmaters
%       numpts: number of elements in discrictized shock vector
%       y: Vecor of Output Solution
%       pi: Vector of Inflation Solution
%       i: Vector of Interest Rate Solution
%       R: Structure related to Rowenhourst Approx Method
%       regime: string denoting which regime to check against (target or
%       deflationary)
%   Outputs:
%       diff: average difference between target and deflationary regime PFs

cPHIpi  = params(4);
cRstar  = params(5);

if strcmp(regime,'tr')  
    y_dr = y;
    pi_dr = pi;
    i_dr = i;
    
    % Get initial matrix and solution vector
    [A, b, out] = eqmmat_tr(params,R,numpts);

    y_tr = out(1:numpts);
    pi_tr = out(numpts+1:2*numpts);
    i_tr = out(2*numpts+1:3*numpts);


    % Refine Matrix to account for ZLB
    if sum(i_tr >= 0) == numpts
        converged = 1;
    else
        converged = 0;
    end

    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_tr(params,R,numpts,A,b);

        i_tr = out_up(2*numpts+1:3*numpts);

        if sum(i_tr >= 0) == numpts
            converged =1;

            y_tr = out_up(1:numpts);
            pi_tr = out_up(numpts+1:2*numpts);
            i_tr = out_up(2*numpts+1:3*numpts);
        end

        A = A_up;
        b = b_up;
    end
    
    diff_y = mean(abs(y_dr - y_tr));
    diff_pi = mean(abs(pi_dr - pi_tr));
    diff_i = mean(abs(i_dr - i_tr));
    
    diff = mean([diff_y,diff_pi,diff_i]);
    
    if sum(i_dr<0) > 0
        diff = 0;
    end
elseif strcmp(regime,'dr') 
    y_tr = y;
    pi_tr = pi;
    i_tr = i;
    
        
    % Get initial matrix and solution vector
    [A, b, out] = eqmmat_dr(params,R,numpts);

    y_dr = out(1:numpts);
    pi_dr = out(numpts+1:2*numpts);
    i_dr = out(2*numpts+1:3*numpts);

    i_real = cRstar + cPHIpi*pi_dr;

    % Refine Matrix to account for ZLB
    if sum(i_dr > 0) == sum(i_real > 0)
        converged = 1;
    else
        converged = 0;
    end
    
    while converged == 0
        [A_up, b_up, out_up] = eqmrefine_dr(params,R,numpts,A,b);

        pi_dr = out_up(numpts+1:2*numpts);
        i_dr = out_up(2*numpts+1:3*numpts);
        i_real = cRstar + cPHIpi*pi_dr;


        if sum(i_dr > 0) == sum(i_real > 0)
            converged =1;

            y_dr = out_up(1:numpts);
            pi_dr = out_up(numpts+1:2*numpts);
            i_dr = out_up(2*numpts+1:3*numpts);
        end

        A = A_up;
        b = b_up;
    end
    
    diff_y = mean(abs(y_dr - y_tr));
    diff_pi = mean(abs(pi_dr - pi_tr));
    diff_i = mean(abs(i_dr - i_tr));
    
    diff = mean([diff_y,diff_pi,diff_i]);
    
    if sum(i_dr<0) > 0
        diff = 0;
    end
end