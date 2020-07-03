function calculateConfidenceInterval(fname)
% Load data
conf_data = load([fname '_logp_conf.mat']);
opt_data = load([fname '_logp.mat']);

% Obtained from the Excel F.INV function using the approprate # of
% parameters and worst-case # of degrees of freedom (2800).
F_STAT = 1.43;

% Worst case number of data points in the OptoNB dataset.
n_data = 2800; 
n_params = length(opt_data.pfit); % # of parameters

DF = n_data - n_params; % # of degrees of freedom

F_thr = 1 + n_params / DF * F_STAT;

% Now interpolate and print out the confidence intervals for the 
% sum-of-squared error curves in kon, koff and Kd:
disp(fname)
print_CI(F_thr, conf_data.SSE_kon, conf_data.kon_vec, 'kon')
print_CI(F_thr, conf_data.SSE_koff, conf_data.koff_vec, 'koff')
print_CI(F_thr, conf_data.SSE_Kd, conf_data.Kd_vec, 'Kd')


function print_CI(F_thr, chi2, p_vec, name)

imin = find(chi2 == min(chi2)); % minimum chi2 is the best-fit point
chi2n = chi2 / min(chi2);       % change if you change sampling!

% All k_on, k_off params are log-transformed during optimization, so
% exponentiate to get the true value.
p_opt = 10.^p_vec(imin);

% Get the indep and dependent variable for finding parameter values where 
% the normalized chi2 crosses F_thr:
X = 10.^(p_vec); % independent variable: parameter value
Y = chi2n;       % dependent variable: normalized chi2

% Plot for fun if you want to see how SSE varies with parameter value:
plot(X, Y, '.-')
% drawnow
pause

% Now get the lower & upper bounds of the 95% confidence:
LB = interp1(Y(1:imin), X(1:imin), F_thr); % left-hand CI 
UB = interp1(Y(imin:end), X(imin:end), F_thr); % right-hand CI

% The 95% confidence interval for k_on:
CI = max(p_opt-LB, UB-p_opt);

fprintf('%s = %0.3g +/- %0.3g\n', name, p_opt, CI)