function getConfIntervals(fname)
% function getConfIntervals(fname)
% 
% Gets confidence intervals associated with a dataset.
% 
% "fname" is the filename prefix for the dataset to analyze. 
% This function assumes you ALREADY ran the parameter optimization function
% runOctetFits_logp to obtain a best-fit parameter set (which was saved in
% the filename [fname '_logp.mat'].
% 
% This function saves its output in a file [fname '_conf.mat'], after
% systematically varying k_on, k_off and K_D and obtaining the
% sum-of-squared errors in each case.

%% load optimization results
opt_data = load(fname);

% fitting functions used to identify the optimal parameter values
F = opt_data.F;
Fplot = opt_data.Fplot;

% best-fit parameters from runOctetFits_logp.m
pfit = opt_data.pfit;

% set some optimization options, just as in runOctetFits_logp.m
options = optimoptions('fmincon');
options.Display = 'iter';
options.MaxFunctionEvaluations = 1e6;
options.OutputFcn = @(p,a,b) Fplot(p);

% optimization lower & upper bounds, just as in runOctetFits_logp.m
LB = pfit*0;
UB = pfit*1e3;
LB(1:3) = pfit(1:3)-3;
UB(1:3) = pfit(1:3)+3;

% What values of the parameter should we sweep through?
% param_ratio goes from 1/4x to 4x of the best-fit parameter value, and
% samples 21 values in this range:

% Nvals = 21;  % # of values to fix each parameter.
% foldchg = 4; % fold change in parameter values above & below
% param_ratio = logspace(log10(1/foldchg), log10(foldchg), Nvals);
scale = [1.05 1.1 1.2 1.4 1.8 2.2 2.6 3.0 3.5 4.0];
param_ratio  = [1./(scale(end:-1:1)) 1 scale]; 

%% STEP 1: Fix kon at various sub-optimal values and recompute SSE
kon_vec = pfit(1) + log10(param_ratio);

% now we are using EQUALITY constraints to fix kon at defined values
% The equality constraint is kon = kon_fixed, so it can be written in the
% form Aeq*p = beq, where Aeq is a vector with a 1 at position 1, and beq
% is the desired value of kon.
Aeq = zeros(1,length(pfit));
Aeq(1) = 1;

for i = 1:length(kon_vec)
    % For the equality constraint to fix k_on:
    beq = kon_vec(i);
    
    % Set initial parameter guess to the optimum ... 
    p0 = pfit;
    % ... but set initial guess at kon to satisfy the constraint:
    p0(1) = beq;
    
    % Perform the optimization!
    [P_kon(:,i) SSE_kon(i)] = fmincon(F, pfit, [], [], Aeq, beq, LB, UB, [], options);
end

%% Run fits with constrained koff to get confidence intervals
koff_vec = pfit(2) + log10(param_ratio);

% now we are using EQUALITY constraints to fix kon at defined values
% The equality constraint is kon = kon_fixed, so it can be written in the
% form Aeq*p = beq, where Aeq is a vector with a 1 at position 1, and beq
% is the desired value of kon.

Aeq = zeros(1,length(pfit));
Aeq(2) = 1;

for i = 1:length(koff_vec)
    beq = koff_vec(i);
    
    p0 = pfit;
    p0(2) = beq;
    
    [P_koff(:,i) SSE_koff(i)] = fmincon(F, pfit, [], [], Aeq, beq, LB, UB, [], options);
end

%% Run fits with constrained Kd to get confidence intervals
% What is Kd? Well, recall we are using LOG-TRANSFORMED koff and kon, so
% Kd = koff / kon -> log(Kd) = log(koff) - log(kon)
% That is how we get this equation for the log-transformed Kd:
Kd_vec = (pfit(2) - pfit(1)) + log10(param_ratio);

% now we are using EQUALITY constraints to fix Kd at defined values.
% The equality constraint is Kd = Kd_fixed. But because we are using LOG 
% parameters note:
% log(Kd) = log(koff) - log(kon)
% so log(Kd) = log(Kd_fixed) =>
% log(koff) - log(kon) = log(Kd_fixed)
% or [-1 1]*[log(kon); log(koff)] = log(Kd_fixed)
% so Aeq * p = beq for Aeq = [-1 1 0 0 0...] and beq = log(Kd_fixed)

Aeq = zeros(1,length(pfit));
Aeq([1 2]) = [-1 1];

for i = 1:length(Kd_vec)
    beq = Kd_vec(i);
    
    % make sure initial point is feasible
    p0 = pfit;
    p0(1) = p0(2) - beq; % set koff-kon = beq => kon = koff-beq
    
    [P_Kd(:,i) SSE_Kd(i)] = fmincon(F, pfit, [], [], Aeq, beq, LB, UB, [], options);
end

%% Save the output!
save([fname '_conf.mat'], 'pfit', 'F', 'Fplot', 'kon_vec', 'P_kon', 'SSE_kon', 'koff_vec', 'P_koff', 'SSE_koff', 'Kd_vec', 'P_Kd', 'SSE_Kd')
