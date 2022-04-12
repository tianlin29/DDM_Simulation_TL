%% Simulation_switchcost
%
% This scirpt is based on 'Simulation1_fixed_and_RT_task'.
%
% If you only want to run 1 simualtion, then run '%% Run simulation of
% params 1'.
%
% If you want to run 2 simulations and compare the 2 results, please run '%%
% Run simulation of params 1' and '%% Run simulation of params 2 (need to
% run after simulation 1)' in sequence.
%
% Tianlin Luo, 2022/4/9

%% Simulation1_fixed_and_RT_task
%
% This script shows psychophysical kernels under fixed duration and reaction
% time task with and without non-decision time (corresponds to Fig. 3 in
% the paper).
%
% The original simulation is time consuming (runs 10^6 trials). As a
% workaround, this script runs a smaller number of trials (10^4) but applys
% a boxcar smoothing (50 ms) to reduce the noise in kernel. 
% No smoothing should be applied for the actual simulation.
%
% Copyright, Kiani Lab, NYU

clear;
%% Parameters
Niters = 3000; % number of simulated trials
smoothing_wid = 50; % boxcar smoothing (50 ms)

%% Parameters of switch cost task (RT task with non-decision time)
% set up the first set of params (non-switch trials)
clear p1;
p1.coh_list = [-0.96 -0.48 -0.24 -0.12 -0.06 0 0.06 0.12 0.24 0.48 0.96]; % stimulus strength
p1.iters = Niters;
p1.termination_rule = {'RT', NaN};
p1.t_max = 5000; % number of simulated time steps in each trial (ms)
p1.B = [-20 20]; % lower and upper bounds (1 x 2) or (t_max x 2)
p1.non_dec_time = 500; % average non-decision time
p1.non_dec_time_sd = 100; % SD of non-decision time
p1.cut_off_RT = 1000; % explicitly determine cut off RT (otherwise, median RT will be the cut off time)
p1.stim_noise = 0; % noise added to the stimulus fluctuation (sensory noise)
p1.dec_noise = 1; % noise added to decision variable (decision noise)
p1.k = 0.3;
% p1.k_max = 0.3; % k_max
% p1.t_k = 500; % time to reach k_max
% p1.k = [linspace(0,p1.k_max,p1.t_k), repmat(p1.k_max, [1,p1.t_max-p1.t_k])]; % k is a function of time

% set up the second set of params (switch trials)
clear p2;
p2.coh_list = [-0.96 -0.48 -0.24 -0.12 -0.06 0 0.06 0.12 0.24 0.48 0.96];
p2.iters = Niters;
p2.termination_rule = {'RT', NaN};
p2.t_max = 5000;
p2.B = [-20 20];
p2.non_dec_time = 500;
p2.non_dec_time_sd = 100;
p2.cut_off_RT = 1000;
p2.stim_noise = 1; % noise added to the stimulus fluctuation (sensory noise)
p2.dec_noise = 0;
p2.k = 0.3;
% p2.k_max = 0.3; % k_max
% p2.t_k = 500; % time to reach k_max
% p2.k = [linspace(0,p2.k_max,p2.t_k), repmat(p2.k_max, [1,p2.t_max-p2.t_k])]; % k is a function of time
% plot(p2.k)

%% Run simulation of params 1
resp = []; rt = []; coh = [];
for k = 1:length(p1.coh_list)
    fprintf('coh = %4.2f. The 1st simulation: running RT task with non-decision time...\n',p1.coh_list(k));
    p1.coh = p1.coh_list(k);
    RT_ndec_sim_1 = DDM_Kernel_Simulation(p1);
    
    resp = [resp; RT_ndec_sim_1.resp];
    rt = [rt; RT_ndec_sim_1.rt];
    coh = [coh; RT_ndec_sim_1.coh];
    
    if p1.coh_list(k)==0
        RT_ndec_sim_1_coh0 = RT_ndec_sim_1;        
    end
end

%% Show figure of simulation 1
% I = ~isnan(rt) & ~isnan(resp);
% fh1 = show_psych_fun(coh(I), resp(I));
% fh2 = show_chrono_fun(coh(I), rt(I));
% 
fh3 = figure;
xrange = [0 1200];
yrange = [-max(p1.k), 0, max(p1.k), 2*max(p1.k)]; % I change the yrange from [0 3] to [-k 2k]
subplot(1,2,1);
show_kernel(RT_ndec_sim_1_coh0, 'stim', smoothing_wid, xrange, yrange);
title('Simulation 1: RT task with non-decision time');
subplot(1,2,2);
show_kernel(RT_ndec_sim_1_coh0, 'resp', smoothing_wid, xrange, yrange);

%% Run simulation of params 2 (need to run after simulation 1)
for k = 1:length(p2.coh_list)
    fprintf('coh = %4.2f. The 2nd simulation: running RT task with non-decision time...\n',p1.coh_list(k));
    p2.coh = p2.coh_list(k);
    RT_ndec_sim_2 = DDM_Kernel_Simulation(p2);
    
    resp = [resp; RT_ndec_sim_2.resp];
    rt = [rt; RT_ndec_sim_2.rt];
    coh = [coh; RT_ndec_sim_2.coh];
    
    if p2.coh_list(k)==0
        RT_ndec_sim_2_coh0 = RT_ndec_sim_2;        
    end
end

%% Show figure of simulation 1&2
cond = [ones(length(coh)/2,1); 2*ones(length(coh)/2,1)];
opt.legend = {'non-switch', 'switch'};
I = ~isnan(rt) & ~isnan(resp);
fh4 = show_psych_fun_2cond(cond(I), coh(I), resp(I), opt);
fh5 = show_chrono_fun_2cond(cond(I), coh(I), rt(I), opt); format_panel()

fh6 = figure;
xrange = [0 1200];
yrange = [-max(p2.k), 0, max(p2.k), 2*max(p2.k)];
subplot(1,2,1);
show_kernel(RT_ndec_sim_2_coh0, 'stim', smoothing_wid, xrange, yrange);
title('Simulation 2: RT task with non-decision time');
subplot(1,2,2);
show_kernel(RT_ndec_sim_2_coh0, 'resp', smoothing_wid, xrange, yrange);

