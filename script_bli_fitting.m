% **************************************
% NANOBODY BLI MEASUREMENTS
% **************************************

% filepath fp:
fp = '.\data_files';

%% WT LaM8
fn = 'LaM8_NB_Experiment_1.txt';
p.conc = [0.6 0.5 0.4 0.3 0.150 0.075 0.0375 0];
p.maxt1 = 80; 
p.maxt2 = 300;
p.hmin = 0.1; 
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% LaM8_GG15_wt_full measurement_in light
fn = 'LaM8_GG15_wt_full measurement_in light.txt';
p.conc = [10 6 3 1 0.6 0.4 0.2 0]; 
p.maxt1 = 20; 
p.maxt2 = 50;
p.hmin = 0.1; 
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% LaM8_GG15_C450V_Experiment_1_100ugml on sensor
fn = 'LaM8_GG15_C450V_Experiment_1_100ugml on sensor.txt';
p.conc = [3 1 0.6 0.4 0.2 0.1 0.05 0];
p.maxt1 = 50; 
p.maxt2 = 300;
p.hmin = 0.1; 
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% DONE: LaM8_GG15_I532E_A536E
fn = 'LaM8_GG15_I532E_A536E.txt';
p.conc = [10 6 3 1 0.6 0.4 0.2 0]; 
p.maxt1 = 20; 
p.maxt2 = 50;
p.hmin = 0.1; 
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% LaM8_AK74_C450V_Experiment_2_100ugml on sensor
% *************** THIS IS PROBLEMATIC. CANNOT DO IT HERE.
fn = 'LaM8_AK74_C450V_Experiment_2_100ugml on sensor.txt';
p.maxt1 = 20;
p.maxt2 = 50;
p.conc = [10 6 3 1 0.6 0.4 0.2 0];
p.hmin = 0.1;
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% LaM8_AK74_I532E_A536E
fn = 'LaM8_AK74_I532E_A536E.txt';
p.conc = [6 3 1 0.6 0.4 0.2 0.1 0];
p.maxt1 = 20;
p.maxt2 = 50;
p.hmin = 0.1;
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% LaM8_AK74_wt_full measurement_in light ANALYZED
fn = 'LaM8_AK74_wt_full measurement_in light ANALYZED.txt';
p.conc = [10 6 3 1 0.6 0.4 0.2 0];
p.maxt1 = 20;
p.maxt2 = 50;
p.hmin = 0.1;
pfit = runOctetFits_logp(fp, fn, p);
getConfIntervals([fn(1:end-4) '_logp'])

%% get all confidence intervals

calculateConfidenceInterval('LaM8_NB_Experiment_1')
calculateConfidenceInterval('LaM8_GG15_wt_full measurement_in light')
calculateConfidenceInterval('LaM8_GG15_C450V_Experiment_1_100ugml on sensor')
calculateConfidenceInterval('LaM8_GG15_I532E_A536E')
calculateConfidenceInterval('LaM8_AK74_C450V_Experiment_2_100ugml on sensor')
calculateConfidenceInterval('LaM8_AK74_I532E_A536E')
calculateConfidenceInterval('LaM8_AK74_wt_full measurement_in light ANALYZED')

