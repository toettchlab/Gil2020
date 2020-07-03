function pfit = runOctetFits(fpath, fname, p)
%% Load data
tmp = importdata([fpath '\' fname]);

%% Extract raw Octet traces from data file
t = tmp.data(:,1);
x = tmp.data(:,2:2:end);

x = x(:,any(x));

for i = 1:size(x,2)
    xh = imhmin(x, p.hmin);
end

%% Automatically identify switches between conditions 
% There is always a jump in BLI measurement value when the probe tips move
% between conditions, so we can identify the switches from the curve
% derivatives...

xm = max(abs(diff(xh)),[],2);
[v I] = sort(xm, 'descend');

for k = 1:20 % remove nearby ones
    for m = 1:50
        ii = find(I(k+1:end) == I(k)+m);
        I(k+ii) = nan;
    end
end
I = I(~isnan(I));

%% Plot traces to make sure condition switches are correct

i_s = sort(I(1:4))+1;
ts = t(i_s);

figure(1),clf
set(gcf, 'position', [0 420 662 190])
plot(t, x)
for i = 1:length(ts)
    line([1 1]*ts(i), [0 10])
end
title(fname(1:end-4), 'interpreter', 'none')
disp('Check that condition switches look good!')
pause(1)
% keyboard % <- just in case you need to manually adjust things...

%% Truncate full traces to just the binding and unbinding curves
% The 3rd & 4th (final) phases are association & dissociation, respectively
t1  = t(i_s(3):i_s(4)-1); 
x1  = x(i_s(3):i_s(4)-1,:); 

% Truncate the association phase to maxt1, which should be set to capture
% the full approach to steady state without being so long that the
% association phase is only captured by a small fraction of data points.
i1 = find(t1 <= t1(1) + p.maxt1);
t1 = t1(i1);
x1 = x1(i1,:);

t2 = t(i_s(4):end);
x2 = x(i_s(4):end,:);

% Truncate the dissociation phase to maxt2, which should be set to capture
% the full approach to steady state without being so long that the
% association phase is only captured by a small fraction of data points.
i2 = find(t2 <= t2(1) + p.maxt2);
t2 = t2(i2);
x2 = x2(i2,:);

% Combine to get the total assocation + dissociation trajectories if
% desired.
t12 = [t1; t2];
x12 = [x1; x2];

% Plot the results
figure(2),clf
set(gcf, 'position', [25 104 366 541])
subplot(3,1,1)
plot(t1, x1)
subplot(3,1,2)
plot(t2, x2)
subplot(3,1,3)
plot(t12, x12)
disp('This is the BLI dataset to analyze...')
pause(1)


%% Construct the objective function for fitting
% F can be used to evaluate the sum-of-squared error of a parameter set
% p.conc to the data in (t1,x1) and (t2,x2)
F     = construct_kinetics_fit_logp(p.conc, t1, x1, t2, x2);
% Fplot just plots the model curves over the data points for a particular 
% parameter set
Fplot = construct_kinetics_fit_logp(p.conc, t1, x1, t2, x2, 1);

%% Get all parameters
% p0 is the initial parameter guess

ka0   = 1e-3; % 1e-2 sometimes works better as initial guess
kd0   = 10^-0.5; 
k20   = 1e-3; % 1e-2 sometimes works better as initial guess

% Make a rough guess at initial offsets and amplitudes of BLI data
p0 = [0 0 0];
for i = 1:size(x1,2)
    bon0  = x1(1,i);
    aon0  = (x1(end,i)-bon0)*exp(k20*(t1(end)-t1(1)));
    boff0 = x2(end,i)*exp(k20*(t2(end)-t2(1)));
    aoff0 = x2(1,i)-boff0;
    p0 = [p0 aon0 bon0 aoff0 boff0];
end

% If any guesses are negative, set to 0.1
p0(p0 < 0) = 0.1;

% Now put initial guess of the log-kon, log-koff, and log-k2 in the p0
% vector
p0(1:3) = log10([ka0 kd0 k20]);

Fplot(p0)

%% Run the optimization
% Set some optimization options
options = optimoptions('fmincon');
options.Display = 'iter';
options.MaxFunctionEvaluations = 1e6;
options.OutputFcn = @(p,a,b) Fplot(p);

% Get lower and upper bounds for fit
% Assume most parameters are positive and are not more than 1,000-fold off
LB = p0*0;
UB = p0*1e3;

% Exception: log-kon, log-koff and log-k2 parameters CAN be negative but
% constrain to 3 orders of magnitude above & below
LB(1:3) = p0(1:3)-3;
UB(1:3) = p0(1:3)+3;

% Run the fit!
pfit = fmincon(F, p0, [], [], [], [], LB, UB, [], options);

% Save the output!
save([fname(1:end-4) '_logp.mat'], 'pfit', 't12', 'x12', 'F', 'Fplot');

