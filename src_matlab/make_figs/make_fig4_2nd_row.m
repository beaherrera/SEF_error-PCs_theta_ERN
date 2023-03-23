%% create laminar power plots for sim lfps and cal stats

clear
clc

rng('default')

tspan = -500:1000;

%% electrodes

h = 150; % um, inter-electrode distance
Ne = 16; % number of electrodes
ze = 0:h:(Ne-1)*h;
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

%% path to data

path2data = fullfile("data/sim_data/processed_lfps");

path_simData_Go_L3 = ['D:\Theta_paper_sim\' ...
    'results_L3PCsPopMky\Go_trial_Aug30_mpi\' ...
    'neurons#1000_clustered_synp\StimDend#4_StimApic#4'];
path_simData_Go_L5 = ['D:\Theta_paper_sim\results_L5PCsPopMky\' ...
    'Go_trial_Oct10_bsinc4_dend_a0_5_spks2_2_dm70_120_sg140_250_apic_a2_spks1_d100_sg200_mpi\' ...
    'neurons#1000_clustered_synp\StimDend#4_StimOblq#0_StimApic#4'];

path_simData_NC_L3 = ['D:\Theta_paper_sim\results_L3PCsPopMky\' ...
    'NC_trial_Aug30_mpi\neurons#1000_clustered_synp\' ...
    'StimDend#4_StimApic#4'];
path_simData_NC_L5 = ['D:\Theta_paper_sim\results_L5PCsPopMky\' ...
    'NC_trial_Oct10_bsinc4_dend_a0_5_spks2_5_dm70_120_sg140_250_apic_a2_5_spks1_d100_280_sg200_mpi\' ...
    'neurons#1000_clustered_synp\StimDend#4_StimOblq#0_StimApic#4'];

file_name_L3 = 'Dend_r3.5_Apic_r2';
file_name_L5 = 'Dend_r2_Apic_r0.5';

%% load events time | ms
trials_number = 1:20;

load(fullfile(path_simData_Go_L3, ...
    ['events_timing_' file_name_L3 '.mat']))
saccade_times_Go = double(saccade_times(1,trials_number))';
target_times_Go = double(target_times(1,trials_number))';

load(fullfile(path_simData_NC_L3, ...
    ['events_timing_' file_name_L3 '.mat']))
saccade_times_NC = double(saccade_times(1,trials_number))';
target_times_NC = double(target_times(1,trials_number))';

% convert event times to sec
target_times_secs_Go = double(target_times_Go).*1e-3; % sec, target times
saccade_times_secs_Go = double(saccade_times_Go).*1e-3; % sec, saccade times
target_times_secs_NC = double(target_times_NC).*1e-3; % sec, target times
saccade_times_secs_NC = double(saccade_times_NC).*1e-3; % sec, saccade times

%% load data

load(fullfile(path2data, 'filtered_sim_lfps_Target.mat'), 'ft_lfpGo_theta_Target', ...
    'ft_lfpNC_theta_Target', 'ft_lfpGo_alpha_Target', ...
    'ft_lfpNC_alpha_Target', 'ft_lfpGo_beta_Target', ...
    'ft_lfpNC_beta_Target', 'ft_lfpGo_gamma_Target', ...
    'ft_lfpNC_gamma_Target')

load(fullfile(path2data, 'filtered_sim_lfps_Saccade.mat'), 'ft_lfpGo_theta_Saccade', ...
    'ft_lfpNC_theta_Saccade', 'ft_lfpGo_alpha_Saccade', ...
    'ft_lfpNC_alpha_Saccade', 'ft_lfpGo_beta_Saccade', ...
    'ft_lfpNC_beta_Saccade', 'ft_lfpGo_gamma_Saccade', ...
    'ft_lfpNC_gamma_Saccade')

%% calculate the Hilbert Transforms

fun_Hilbert = @(cell_array) cellfun(@(x) hilbert(x')', cell_array, ...
    'UniformOutput',false);

lfpGo_theta_Hilbert_Target = fun_Hilbert(ft_lfpGo_theta_Target.trial);
lfpGo_alpha_Hilbert_Target = fun_Hilbert(ft_lfpGo_alpha_Target.trial);
lfpGo_beta_Hilbert_Target = fun_Hilbert(ft_lfpGo_beta_Target.trial);
lfpGo_gamma_Hilbert_Target = fun_Hilbert(ft_lfpGo_gamma_Target.trial);

lfpGo_theta_Hilbert_Saccade = fun_Hilbert(ft_lfpGo_theta_Saccade.trial);
lfpGo_alpha_Hilbert_Saccade = fun_Hilbert(ft_lfpGo_alpha_Saccade.trial);
lfpGo_beta_Hilbert_Saccade = fun_Hilbert(ft_lfpGo_beta_Saccade.trial);
lfpGo_gamma_Hilbert_Saccade = fun_Hilbert(ft_lfpGo_gamma_Saccade.trial);

lfpNC_theta_Hilbert_Target = fun_Hilbert(ft_lfpNC_theta_Target.trial);
lfpNC_alpha_Hilbert_Target = fun_Hilbert(ft_lfpNC_alpha_Target.trial);
lfpNC_beta_Hilbert_Target = fun_Hilbert(ft_lfpNC_beta_Target.trial);
lfpNC_gamma_Hilbert_Target = fun_Hilbert(ft_lfpNC_gamma_Target.trial);

lfpNC_theta_Hilbert_Saccade = fun_Hilbert(ft_lfpNC_theta_Saccade.trial);
lfpNC_alpha_Hilbert_Saccade = fun_Hilbert(ft_lfpNC_alpha_Saccade.trial);
lfpNC_beta_Hilbert_Saccade = fun_Hilbert(ft_lfpNC_beta_Saccade.trial);
lfpNC_gamma_Hilbert_Saccade = fun_Hilbert(ft_lfpNC_gamma_Saccade.trial);

%% calculate amplitude and phase envelope

cal_mag = @(Hilbert_cell) cellfun(@abs, Hilbert_cell, ...
    'UniformOutput',false);

lfpGo_theta_amp_Target = cal_mag(lfpGo_theta_Hilbert_Target);
lfpNC_theta_amp_Target = cal_mag(lfpNC_theta_Hilbert_Target);
lfpGo_theta_amp_Saccade = cal_mag(lfpGo_theta_Hilbert_Saccade);
lfpNC_theta_amp_Saccade = cal_mag(lfpNC_theta_Hilbert_Saccade);

lfpGo_alpha_amp_Target = cal_mag(lfpGo_alpha_Hilbert_Target);
lfpNC_alpha_amp_Target = cal_mag(lfpNC_alpha_Hilbert_Target);
lfpGo_alpha_amp_Saccade = cal_mag(lfpGo_alpha_Hilbert_Saccade);
lfpNC_alpha_amp_Saccade = cal_mag(lfpNC_alpha_Hilbert_Saccade);

lfpGo_beta_amp_Target = cal_mag(lfpGo_beta_Hilbert_Target);
lfpNC_beta_amp_Target = cal_mag(lfpNC_beta_Hilbert_Target);
lfpGo_beta_amp_Saccade = cal_mag(lfpGo_beta_Hilbert_Saccade);
lfpNC_beta_amp_Saccade = cal_mag(lfpNC_beta_Hilbert_Saccade);

lfpGo_gamma_amp_Target = cal_mag(lfpGo_gamma_Hilbert_Target);
lfpNC_gamma_amp_Target = cal_mag(lfpNC_gamma_Hilbert_Target);
lfpGo_gamma_amp_Saccade = cal_mag(lfpGo_gamma_Hilbert_Saccade);
lfpNC_gamma_amp_Saccade = cal_mag(lfpNC_gamma_Hilbert_Saccade);

%% select window -.5 to 1 sec relative to saccade/target and store values

theta_trials_Saccade = {lock2event(lfpGo_theta_amp_Saccade, ft_lfpGo_theta_Saccade.time),...
    lock2event(lfpNC_theta_amp_Saccade, ft_lfpNC_theta_Saccade.time)};
alpha_trials_Saccade = {lock2event(lfpGo_alpha_amp_Saccade, ft_lfpGo_alpha_Saccade.time),...
    lock2event(lfpNC_alpha_amp_Saccade, ft_lfpNC_alpha_Saccade.time)};
beta_trials_Saccade = {lock2event(lfpGo_beta_amp_Saccade, ft_lfpGo_beta_Saccade.time),...
    lock2event(lfpNC_beta_amp_Saccade, ft_lfpNC_beta_Saccade.time)};
gamma_trials_Saccade = {lock2event(lfpGo_gamma_amp_Saccade, ft_lfpGo_gamma_Saccade.time),...
    lock2event(lfpNC_gamma_amp_Saccade, ft_lfpNC_gamma_Saccade.time)};

target_theta_trials = {lock2event(lfpGo_theta_amp_Target, ft_lfpGo_theta_Target.time),...
    lock2event(lfpNC_theta_amp_Target, ft_lfpNC_theta_Target.time)};
target_alpha_trials = {lock2event(lfpGo_alpha_amp_Target, ft_lfpGo_alpha_Target.time),...
    lock2event(lfpNC_alpha_amp_Target, ft_lfpNC_alpha_Target.time)};
target_beta_trials = {lock2event(lfpGo_beta_amp_Target, ft_lfpGo_beta_Target.time),...
    lock2event(lfpNC_beta_amp_Target, ft_lfpNC_beta_Target.time)};
target_gamma_trials = {lock2event(lfpGo_gamma_amp_Target, ft_lfpGo_gamma_Target.time),...
    lock2event(lfpNC_gamma_amp_Target, ft_lfpNC_gamma_Target.time)};

%% get the power and baseline correct it to 200ms pre target power
mod_data = @(x, x_bs) ((nan2zero(x).*1e6).^2 - ...
    mean((nan2zero(x_bs(:, tspan>=-200 & tspan<0, :)).*1e6).^2,2));

theta_Go_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    theta_trials_Saccade(:,1), target_theta_trials(:,1), 'UniformOutput',false);
theta_Go_cell = cat(3, theta_Go_cell{:});

theta_NC_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    theta_trials_Saccade(:,2), target_theta_trials(:,2), 'UniformOutput',false);
theta_NC_cell = cat(3, theta_NC_cell{:});

alpha_Go_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    alpha_trials_Saccade(:,1), target_alpha_trials(:,1), 'UniformOutput',false);
alpha_Go_cell = cat(3, alpha_Go_cell{:});

alpha_NC_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    alpha_trials_Saccade(:,2), target_alpha_trials(:,2), 'UniformOutput',false);
alpha_NC_cell = cat(3, alpha_NC_cell{:});

beta_Go_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    beta_trials_Saccade(:,1), target_beta_trials(:,1), 'UniformOutput',false);
beta_Go_cell = cat(3, beta_Go_cell{:});

beta_NC_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    beta_trials_Saccade(:,2), target_beta_trials(:,2), 'UniformOutput',false);
beta_NC_cell = cat(3, beta_NC_cell{:});

gamma_Go_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    gamma_trials_Saccade(:,1), target_gamma_trials(:,1), 'UniformOutput',false);
gamma_Go_cell = cat(3, gamma_Go_cell{:});

gamma_NC_cell = cellfun(@(x, x_bs) mod_data(x, x_bs), ...
    gamma_trials_Saccade(:,2), target_gamma_trials(:,2), 'UniformOutput',false);
gamma_NC_cell = cat(3, gamma_NC_cell{:});

%% prepare data for calculating tf stats using fieldtrip

% Go trials
TFR_Go_theta = [];
TFR_Go_theta.dimord = 'rpt_chan_freq_time';
TFR_Go_theta.time = tspan.*1e-3; % time span in seconds
TFR_Go_theta.freq = 1:16; % channels

TFR_Go_alpha = TFR_Go_theta;
TFR_Go_beta = TFR_Go_theta;
TFR_Go_gamma = TFR_Go_theta;

TFR_Go_theta.label = {'theta'}; 
TFR_Go_theta.powspctrm = zeros(size(theta_Go_cell, 3), length(TFR_Go_theta.label), ...
    size(theta_Go_cell, 1), size(theta_Go_cell, 2));
TFR_Go_theta.powspctrm(:, 1, :, :) = permute(theta_Go_cell, [3 1 2]);

TFR_Go_alpha.label = {'alpha'}; 
TFR_Go_alpha.powspctrm = zeros(size(alpha_Go_cell, 3), length(TFR_Go_alpha.label), ...
    size(alpha_Go_cell, 1), size(alpha_Go_cell, 2));
TFR_Go_alpha.powspctrm(:, 1, :, :) = permute(alpha_Go_cell, [3 1 2]);

TFR_Go_beta.label = {'beta'}; 
TFR_Go_beta.powspctrm = zeros(size(beta_Go_cell, 3), length(TFR_Go_beta.label), ...
    size(beta_Go_cell, 1), size(beta_Go_cell, 2));
TFR_Go_beta.powspctrm(:, 1, :, :) = permute(beta_Go_cell, [3 1 2]);

TFR_Go_gamma.label = {'gamma'}; 
TFR_Go_gamma.powspctrm = zeros(size(gamma_Go_cell, 3), length(TFR_Go_gamma.label), ...
    size(gamma_Go_cell, 1), size(gamma_Go_cell, 2));
TFR_Go_gamma.powspctrm(:, 1, :, :) = permute(gamma_Go_cell, [3 1 2]);
TFR_Go_gamma.cfg = [];

% NC trials
TFR_NC_theta = [];
TFR_NC_theta.dimord = 'rpt_chan_freq_time';
TFR_NC_theta.time = tspan.*1e-3; % time span in seconds
TFR_NC_theta.freq = 1:16; % channels

TFR_NC_alpha = TFR_NC_theta;
TFR_NC_beta = TFR_NC_theta;
TFR_NC_gamma = TFR_NC_theta;

TFR_NC_theta.label = {'theta'}; 
TFR_NC_theta.powspctrm = zeros(size(theta_NC_cell, 3), length(TFR_NC_theta.label), ...
    size(theta_NC_cell, 1), size(theta_NC_cell, 2));
TFR_NC_theta.powspctrm(:, 1, :, :) = permute(theta_NC_cell, [3 1 2]);

TFR_NC_alpha.label = {'alpha'}; 
TFR_NC_alpha.powspctrm = zeros(size(alpha_NC_cell, 3), length(TFR_NC_alpha.label), ...
    size(alpha_NC_cell, 1), size(alpha_NC_cell, 2));
TFR_NC_alpha.powspctrm(:, 1, :, :) = permute(alpha_NC_cell, [3 1 2]);

TFR_NC_beta.label = {'beta'}; 
TFR_NC_beta.powspctrm = zeros(size(beta_NC_cell, 3), length(TFR_NC_beta.label), ...
    size(beta_NC_cell, 1), size(beta_NC_cell, 2));
TFR_NC_beta.powspctrm(:, 1, :, :) = permute(beta_NC_cell, [3 1 2]);

TFR_NC_gamma.label = {'gamma'}; 
TFR_NC_gamma.powspctrm = zeros(size(gamma_NC_cell, 3), length(TFR_NC_gamma.label), ...
    size(gamma_NC_cell, 1), size(gamma_NC_cell, 2));
TFR_NC_gamma.powspctrm(:, 1, :, :) = permute(gamma_NC_cell, [3 1 2]);
TFR_NC_gamma.cfg = [];

%% calculate stats between conditions
p_thr = 0.01;
p_cluster_thr = 0.01;
cfg = [];
cfg.latency = 'all';
cfg.parameter = 'powspctrm';
cfg.method  = 'montecarlo';
cfg.correctm= 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.neighbours= [];
cfg.channel = 'all';
cfg.alpha = p_thr;
cfg.correcttail = 'alpha';
cfg.clusteralpha = p_cluster_thr; % culter-level thrshold (wether to include
% it in the cluster or not)
cfg.statistic = 'ft_statfun_indepsamplesT'; %'depsamplesT'; % indepsamplesT
cfg.ivar = 1; % number or list with indices indicating
% the independent variable(s)
cfg.numrandomization = 10000;
cfg.design = [ones(1, size(TFR_Go_theta.powspctrm,1)), ...
    2*ones(1,size(TFR_NC_theta.powspctrm,1))]; % design matrix

stats_NCvsGo_TFR_theta = ft_freqstatistics(cfg,TFR_NC_theta,TFR_Go_theta);
stats_NCvsGo_TFR_alpha = ft_freqstatistics(cfg,TFR_NC_alpha,TFR_Go_alpha);
stats_NCvsGo_TFR_beta = ft_freqstatistics(cfg,TFR_NC_beta,TFR_Go_beta);
stats_NCvsGo_TFR_gamma = ft_freqstatistics(cfg,TFR_NC_gamma,TFR_Go_gamma);

%% load monkey Eu's power

load(fullfile('data', 'monkeys_data', 'sess_avg_power.mat'), ...
    'power_Go_Eu_theta', 'power_NC_Eu_theta')

max_Eu_power = max(abs([power_Go_Eu_theta(:,tspan>=-500 & tspan < tspan_max) ...
    power_NC_theta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

%% create figure

power_Go_theta = 100.*squeeze(mean(squeeze(TFR_Go_theta.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;
power_NC_theta = 100.*squeeze(mean(squeeze(TFR_NC_theta.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;

power_Go_alpha = 100.*squeeze(mean(squeeze(TFR_Go_alpha.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;
power_NC_alpha = 100.*squeeze(mean(squeeze(TFR_NC_alpha.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;

power_Go_beta = 100.*squeeze(mean(squeeze(TFR_Go_beta.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;
power_NC_beta = 100.*squeeze(mean(squeeze(TFR_NC_beta.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;

power_Go_gamma = 100.*squeeze(mean(squeeze(TFR_Go_gamma.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;
power_NC_gamma = 100.*squeeze(mean(squeeze(TFR_NC_gamma.powspctrm), 1, ...
    "omitnan"))./max_Eu_power;

%% plot

p_thr_plot = 0.05;

c_scale = 200;

tspan_max = 500;
rows = 4;
cols = 4;

ERN = median([184 224]);
Pe = round(median([302 327]));

figure('Units', 'inches','Position',[0 0 6 7]);
tiledlayout(4, 3,'TileSpacing','Compact','Padding','Compact');

% -- theta --
ii = 1;
nexttile(ii);
caxis_lim = [-1 1].*max(abs([power_Go_theta(:,tspan>=-500 & tspan < tspan_max) ...
    power_NC_theta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

TFR_plot(tspan, ze, power_Go_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
xticklabels({' ', ' ', ' '})
hold on;
% text(-450, 200, '\theta', 'FontSize', 16, 'FontWeight','bold')
% ylabel('Cortical Depth (\mum)');
% title('Correct')
% set(axi, 'Position')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticks([-500 0 500 1000])
% % xticklabels({'-.5', '0', '.5', '1'})
% xticklabels({})
hold on;
% text(-450, 200, '\theta', 'FontSize', 16, 'FontWeight','bold')
% title('Error')
% c = colorbar;
% c.Label.String = {'\muV^2'};

ii = ii+1;
ax1 = nexttile(ii);
caxis_lim_diff = [-1 1].*max(abs(power_NC_theta(:,tspan>=-500 & tspan < tspan_max) ...
    - power_Go_theta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');

TFR_stats_plot(tspan, ze, power_NC_theta - power_Go_theta, ...
    squeeze(stats_NCvsGo_TFR_theta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
% text(-450, 200, '\theta', 'FontSize', 16, 'FontWeight','bold')
% title('Error - Correct')
% c = colorbar;
% c.Label.String = {'\muV^2'};

% -- alpha --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_alpha(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

TFR_plot(tspan, ze, power_Go_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
% text(-450, 200, '\alpha', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
% text(-450, 200, '\alpha', 'FontSize', 16, 'FontWeight','bold')
% c = colorbar;
% c.Label.String = {'\muV^2'};

ii = ii+1;
ax2 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_alpha(:,tspan>=-500 & tspan < tspan_max)), [], 'all');

TFR_stats_plot(tspan, ze, power_NC_alpha - power_Go_alpha, ...
    squeeze(stats_NCvsGo_TFR_alpha.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- beta --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_beta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

TFR_plot(tspan, ze, power_Go_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
% text(-450, 200, '\beta', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

ii = ii+1;
ax3 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_beta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');

TFR_stats_plot(tspan, ze, power_NC_beta - power_Go_beta, ...
    squeeze(stats_NCvsGo_TFR_beta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- gamma --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_gamma(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
% 
TFR_plot(tspan, ze, power_Go_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
xticklabels({' ', ' ', ' '})
hold on;

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticks([-500 0 500 1000])

ii = ii+1;
ax4 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_gamma(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
% 
TFR_stats_plot(tspan, ze, power_NC_gamma - power_Go_gamma, ...
    squeeze(stats_NCvsGo_TFR_gamma.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
% set colormap for diff plots
[map,~,~,~] = brewermap(256,'PiYG');
colormap(ax4, (map));
colormap(ax1, (map));
colormap(ax2, (map));
colormap(ax3, (map));

%% subfunctions

function TFR_plot(tspan, ze, power, caxis_lim, z_depth, tspan_max)

q = length(tspan); % number of interpolating points
[zeq, tspanq] = meshgrid(linspace(ze(1),ze(end), q), tspan);
power_interp = interp2(tspan, ze, power, tspanq, zeq, 'spline');

imagesc(tspan, zeq(1,:)', power_interp')
xlim([-500 tspan_max])
caxis(caxis_lim)
colormap(jet);
ax = gca; % current axes
ax.FontSize = 9;
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
set(gca,'linewidth',1.5,'fontsize',10,'fontweight','bold','TickDir','out')

end

function TFR_stats_plot(tspan, ze, power, stats_matrix, caxis_lim, z_depth, tspan_max)

q = length(tspan); % number of interpolating points eq(1,:)'
[zeq, tspanq] = meshgrid(linspace(ze(1),ze(end), q), tspan);

stats_interp = interp2(tspan, ze, stats_matrix, tspanq, zeq, 'nearest');
power_interp = interp2(tspan, ze, power, tspanq, zeq, 'spline');

ft_plot_matrix(tspan, zeq(1,:)', power_interp', 'highlight', stats_interp', ...
    'highlightstyle', 'outline')

xlim([-500 tspan_max])
caxis(caxis_lim)
colormap(jet);
ax = gca; % current axes
ax.FontSize = 9;
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
set(gca,'linewidth',1.5,'fontsize',10,'fontweight','bold','TickDir','out')

end
