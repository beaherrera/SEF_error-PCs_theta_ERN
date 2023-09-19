%% calculate stats for laminar power maps and create Fig 4-top and its Suppl Fig (Fig S2)

clear
clc

rng('default')

tspan = -500:1000;

%% path 2 data

path2data = fullfile('data', 'monkeys_data'); % IMPORTANT: you
% have to download and save the laminar data in this folder
path2tfdata_target = fullfile(path2data, 'laminar_power_target'); % path 2
% laminar power relative to the target
path2tfdata_saccade = fullfile(path2data, 'laminar_power_saccade');

%% perpendicular sessions

% sessions number
Sess_EuP1 = 14:19; % Eu, site P1
Sess_XP2P3 = [20:25 26:29]; % X, 20-25 site P2, 26-29 site P3
SessNumb = [Sess_EuP1, Sess_XP2P3];

%% electrodes

h = 150; % um, inter-electrode distance
Ne = 16; % number of electrodes
ze = 0:h:(Ne-1)*h;
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

%% frequencies of interest

thetaBand = [5 8]; % Hz
alphaBand = [9 14]; % Hz
betaBand = [15 29]; % Hz
gammaBand = [30 80]; % Hz

%% load amp envelope for each frequency band

% target
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_theta.mat'), 'target_theta_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_alpha.mat'), 'target_alpha_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_beta.mat'), 'target_beta_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_gamma.mat'), 'target_gamma_trials')

% saccade
load(fullfile(path2tfdata_saccade, ...
    'sessions_amp_phase_theta.mat'), 'theta_trials')
load(fullfile(path2tfdata_saccade, ...
    'sessions_amp_phase_alpha.mat'), 'alpha_trials')
load(fullfile(path2tfdata_saccade, ...
    'sessions_amp_phase_beta.mat'), 'beta_trials')
load(fullfile(path2tfdata_saccade, ...
    'sessions_amp_phase_gamma.mat'), 'gamma_trials')

%% get the power and baseline correct it to 200ms pre target power

mod_data = @(x, x_bs) ((nan2zero(x).*1e6).^2 - ...
    mean((nan2zero(x_bs(:, tspan>=-200 & tspan<0, :)).*1e6).^2,2));

theta_Go_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    theta_trials(:,1), target_theta_trials(:,1), 'UniformOutput',false);
theta_Go_cell = cat(3, theta_Go_cell{:});

theta_NC_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    theta_trials(:,2), target_theta_trials(:,2), 'UniformOutput',false);
theta_NC_cell = cat(3, theta_NC_cell{:});

alpha_Go_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    alpha_trials(:,1), target_alpha_trials(:,1), 'UniformOutput',false);
alpha_Go_cell = cat(3, alpha_Go_cell{:});

alpha_NC_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    alpha_trials(:,2), target_alpha_trials(:,2), 'UniformOutput',false);
alpha_NC_cell = cat(3, alpha_NC_cell{:});

beta_Go_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    beta_trials(:,1), target_beta_trials(:,1), 'UniformOutput',false);
beta_Go_cell = cat(3, beta_Go_cell{:});

beta_NC_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    beta_trials(:,2), target_beta_trials(:,2), 'UniformOutput',false);
beta_NC_cell = cat(3, beta_NC_cell{:});

gamma_Go_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    gamma_trials(:,1), target_gamma_trials(:,1), 'UniformOutput',false);
gamma_Go_cell = cat(3, gamma_Go_cell{:});

gamma_NC_cell = cellfun(@(x, x_bs) mean(mod_data(x, x_bs), 3), ...
    gamma_trials(:,2), target_gamma_trials(:,2), 'UniformOutput',false);
gamma_NC_cell = cat(3, gamma_NC_cell{:});

%% prepare data for calculating tf stats using fieldtrip
% instead of having the
% channels in this dimension, I'll have the frequency bands of interest.
% I'll have the channels in the frequency dimension because I want to
% calculate the stats across channels and time instead of frequencies.

%% Eu

% Go trials
TFR_Go_theta_Eu = [];
TFR_Go_theta_Eu.powspctrmdimord = 'rpt_chan_freq_time';
TFR_Go_theta_Eu.time = tspan.*1e-3; % time span in seconds
TFR_Go_theta_Eu.freq = 1:16; % channels

TFR_Go_alpha_Eu = TFR_Go_theta_Eu;
TFR_Go_beta_Eu = TFR_Go_theta_Eu;
TFR_Go_gamma_Eu = TFR_Go_theta_Eu;

TFR_Go_theta_Eu.label = {'theta'};
TFR_Go_theta_Eu.powspctrm = zeros(size(theta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_Go_theta_Eu.label), ...
    size(theta_Go_cell, 1), size(theta_Go_cell, 2));
TFR_Go_theta_Eu.powspctrm(:, 1, :, :) = permute(theta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_Go_alpha_Eu.label = {'alpha'};
TFR_Go_alpha_Eu.powspctrm = zeros(size(alpha_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_Go_alpha_Eu.label), ...
    size(alpha_Go_cell, 1), size(alpha_Go_cell, 2));
TFR_Go_alpha_Eu.powspctrm(:, 1, :, :) = permute(alpha_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_Go_beta_Eu.label = {'beta'};
TFR_Go_beta_Eu.powspctrm = zeros(size(beta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_Go_beta_Eu.label), ...
    size(beta_Go_cell, 1), size(beta_Go_cell, 2));
TFR_Go_beta_Eu.powspctrm(:, 1, :, :) = permute(beta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_Go_gamma_Eu.label = {'gamma'};
TFR_Go_gamma_Eu.powspctrm = zeros(size(gamma_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_Go_gamma_Eu.label), ...
    size(gamma_Go_cell, 1), size(gamma_Go_cell, 2));
TFR_Go_gamma_Eu.powspctrm(:, 1, :, :) = permute(gamma_Go_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);
TFR_Go_gamma_Eu.cfg = [];

% NC trials
TFR_NC_theta_Eu = [];
TFR_NC_theta_Eu.powspctrmdimord = 'rpt_chan_freq_time';
TFR_NC_theta_Eu.time = tspan.*1e-3; % time span in seconds
TFR_NC_theta_Eu.freq = 1:16; % channels

TFR_NC_alpha_Eu = TFR_NC_theta_Eu;
TFR_NC_beta_Eu = TFR_NC_theta_Eu;
TFR_NC_gamma_Eu = TFR_NC_theta_Eu;

TFR_NC_theta_Eu.label = {'theta'};
TFR_NC_theta_Eu.powspctrm = zeros(size(theta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_NC_theta_Eu.label), ...
    size(theta_NC_cell, 1), size(theta_NC_cell, 2));
TFR_NC_theta_Eu.powspctrm(:, 1, :, :) = permute(theta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_NC_alpha_Eu.label = {'alpha'};
TFR_NC_alpha_Eu.powspctrm = zeros(size(alpha_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_NC_alpha_Eu.label), ...
    size(alpha_NC_cell, 1), size(alpha_NC_cell, 2));
TFR_NC_alpha_Eu.powspctrm(:, 1, :, :) = permute(alpha_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_NC_beta_Eu.label = {'beta'};
TFR_NC_beta_Eu.powspctrm = zeros(size(beta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_NC_beta_Eu.label), ...
    size(beta_NC_cell, 1), size(beta_NC_cell, 2));
TFR_NC_beta_Eu.powspctrm(:, 1, :, :) = permute(beta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);

TFR_NC_gamma_Eu.label = {'gamma'};
TFR_NC_gamma_Eu.powspctrm = zeros(size(gamma_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), 3), length(TFR_NC_gamma_Eu.label), ...
    size(gamma_NC_cell, 1), size(gamma_NC_cell, 2));
TFR_NC_gamma_Eu.powspctrm(:, 1, :, :) = permute(gamma_NC_cell(:,:, ...
    ismember(SessNumb, Sess_EuP1)), [3 1 2]);
TFR_NC_gamma_Eu.cfg = [];

%% X Sess_XP2P3

% Go trials
TFR_Go_theta_X = [];
TFR_Go_theta_X.powspctrmdimord = 'rpt_chan_freq_time';
TFR_Go_theta_X.time = tspan.*1e-3; % time span in seconds
TFR_Go_theta_X.freq = 1:16; % channels

TFR_Go_alpha_X = TFR_Go_theta_X;
TFR_Go_beta_X = TFR_Go_theta_X;
TFR_Go_gamma_X = TFR_Go_theta_X;

TFR_Go_theta_X.label = {'theta'};
TFR_Go_theta_X.powspctrm = zeros(size(theta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_Go_theta_X.label), ...
    size(theta_Go_cell, 1), size(theta_Go_cell, 2));
TFR_Go_theta_X.powspctrm(:, 1, :, :) = permute(theta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_Go_alpha_X.label = {'alpha'};
TFR_Go_alpha_X.powspctrm = zeros(size(alpha_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_Go_alpha_X.label), ...
    size(alpha_Go_cell, 1), size(alpha_Go_cell, 2));
TFR_Go_alpha_X.powspctrm(:, 1, :, :) = permute(alpha_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_Go_beta_X.label = {'beta'};
TFR_Go_beta_X.powspctrm = zeros(size(beta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_Go_beta_X.label), ...
    size(beta_Go_cell, 1), size(beta_Go_cell, 2));
TFR_Go_beta_X.powspctrm(:, 1, :, :) = permute(beta_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_Go_gamma_X.label = {'gamma'};
TFR_Go_gamma_X.powspctrm = zeros(size(gamma_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_Go_gamma_X.label), ...
    size(gamma_Go_cell, 1), size(gamma_Go_cell, 2));
TFR_Go_gamma_X.powspctrm(:, 1, :, :) = permute(gamma_Go_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);
TFR_Go_gamma_X.cfg = [];

% NC trials
TFR_NC_theta_X = [];
TFR_NC_theta_X.powspctrmdimord = 'rpt_chan_freq_time';
TFR_NC_theta_X.time = tspan.*1e-3; % time span in seconds
TFR_NC_theta_X.freq = 1:16; % channels

TFR_NC_alpha_X = TFR_NC_theta_X;
TFR_NC_beta_X = TFR_NC_theta_X;
TFR_NC_gamma_X = TFR_NC_theta_X;

TFR_NC_theta_X.label = {'theta'};
TFR_NC_theta_X.powspctrm = zeros(size(theta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_NC_theta_X.label), ...
    size(theta_NC_cell, 1), size(theta_NC_cell, 2));
TFR_NC_theta_X.powspctrm(:, 1, :, :) = permute(theta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_NC_alpha_X.label = {'alpha'};
TFR_NC_alpha_X.powspctrm = zeros(size(alpha_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_NC_alpha_X.label), ...
    size(alpha_NC_cell, 1), size(alpha_NC_cell, 2));
TFR_NC_alpha_X.powspctrm(:, 1, :, :) = permute(alpha_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_NC_beta_X.label = {'beta'};
TFR_NC_beta_X.powspctrm = zeros(size(beta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_NC_beta_X.label), ...
    size(beta_NC_cell, 1), size(beta_NC_cell, 2));
TFR_NC_beta_X.powspctrm(:, 1, :, :) = permute(beta_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);

TFR_NC_gamma_X.label = {'gamma'};
TFR_NC_gamma_X.powspctrm = zeros(size(gamma_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), 3), length(TFR_NC_gamma_X.label), ...
    size(gamma_NC_cell, 1), size(gamma_NC_cell, 2));
TFR_NC_gamma_X.powspctrm(:, 1, :, :) = permute(gamma_NC_cell(:,:, ...
    ismember(SessNumb, Sess_XP2P3)), [3 1 2]);
TFR_NC_gamma_X.cfg = [];

%% calculate stats
p_thr = 0.01; % 0.05 or 0.01
p_cluster_thr = 0.05;
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

%% Eu
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = [];
Nsubj = length(Sess_EuP1);
design = zeros(2, Nsubj*2);
design(2,:) = [1:Nsubj 1:Nsubj];
design(1,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
cfg.uvar   = 2;
cfg.ivar   = 1; % number or list with indices indicating
% the independent variable(s)
cfg.numrandomization='all';

stats_NCvsGo_Eu_TFR_theta = ft_freqstatistics(cfg,TFR_NC_theta_Eu,TFR_Go_theta_Eu);
stats_NCvsGo_Eu_TFR_beta = ft_freqstatistics(cfg,TFR_NC_beta_Eu,TFR_Go_beta_Eu);
stats_NCvsGo_Eu_TFR_alpha = ft_freqstatistics(cfg,TFR_NC_alpha_Eu,TFR_Go_alpha_Eu);
stats_NCvsGo_Eu_TFR_gamma = ft_freqstatistics(cfg,TFR_NC_gamma_Eu,TFR_Go_gamma_Eu);

%% X

cfg.statistic = 'ft_statfun_depsamplesT';
cfg.design = [];
Nsubj = length(Sess_XP2P3);
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2; % number or list with indices indicating
% the independent variable(s)
cfg.numrandomization='all';

stats_NCvsGo_X_TFR_theta = ft_freqstatistics(cfg,TFR_NC_theta_X,TFR_Go_theta_X);
stats_NCvsGo_X_TFR_alpha = ft_freqstatistics(cfg,TFR_NC_alpha_X,TFR_Go_alpha_X);
stats_NCvsGo_X_TFR_beta = ft_freqstatistics(cfg,TFR_NC_beta_X,TFR_Go_beta_X);
stats_NCvsGo_X_TFR_gamma = ft_freqstatistics(cfg,TFR_NC_gamma_X,TFR_Go_gamma_X);

%% create figure

mean_fun = @(TFR) squeeze(mean(squeeze(TFR.powspctrm), 1, "omitnan"));

power_Go_Eu_theta = mean_fun(TFR_Go_theta_Eu.powspctrm);
power_Go_X_theta = mean_fun(TFR_Go_theta_X.powspctrm);

power_NC_Eu_theta = mean_fun(TFR_NC_theta_Eu.powspctrm);
power_NC_X_theta = mean_fun(TFR_NC_theta_X.powspctrm);

power_Go_Eu_alpha = mean_fun(TFR_Go_alpha_Eu.powspctrm);
power_Go_X_alpha = mean_fun(TFR_Go_alpha_X.powspctrm);

power_NC_Eu_alpha = mean_fun(TFR_NC_alpha_Eu.powspctrm);
power_NC_X_alpha = mean_fun(TFR_NC_alpha_X.powspctrm);

power_Go_Eu_beta = mean_fun(TFR_Go_beta_Eu.powspctrm);
power_Go_X_beta = mean_fun(TFR_Go_beta_X.powspctrm);

power_NC_Eu_beta = mean_fun(TFR_NC_beta_Eu.powspctrm);
power_NC_X_beta = mean_fun(TFR_NC_beta_X.powspctrm);

power_Go_Eu_gamma = mean_fun(TFR_Go_gamma_Eu.powspctrm);
power_Go_X_gamma = mean_fun(TFR_Go_gamma_X.powspctrm);

power_NC_Eu_gamma = mean_fun(TFR_NC_gamma_Eu.powspctrm);
power_NC_X_gamma = mean_fun(TFR_NC_gamma_X.powspctrm);

%% save

save(fullfile(path2data, 'sess_avg_power.mat'), 'power_NC_*_*', ...
    'power_Go_*_*')

%% plot

ERN = median([184 224]); % ERN peak time
Pe = round(median([302 327])); % Pe peak time

%% ==================== Eu ====================

p_thr_plot = 0.05;

c_scale = 200;

tspan_max = 500;
rows = 4;
cols = 4;

figure('Units', 'inches','Position',[0 0 6 7]);
tiledlayout(4, 3,'TileSpacing','Compact','Padding','Compact');

% -- theta --
ii = 1;
nexttile(ii);
caxis_lim = [-1 1].*max(abs([power_Go_Eu_theta(:,tspan>=-500 & tspan < tspan_max) ...
    power_NC_theta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

TFR_plot(tspan, ze, power_Go_Eu_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
xticklabels({' ', ' ', ' '})
hold on;
text(-450, 200, '\theta', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_Eu_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticks([-500 0 500 1000])
xticklabels({})
hold on;

ii = ii+1;
ax1 = nexttile(ii);
caxis_lim_diff = [-1 1].*max(abs(power_NC_Eu_theta(:,tspan>=-500 & tspan < tspan_max) ...
    - power_Go_Eu_theta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
TFR_stats_plot(tspan, ze, power_NC_Eu_theta - power_Go_Eu_theta, ...
    squeeze(stats_NCvsGo_Eu_TFR_theta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- alpha --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_Eu_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_Eu_alpha(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
TFR_plot(tspan, ze, power_Go_Eu_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
text(-450, 200, '\alpha', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_Eu_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

ii = ii+1;
ax2 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_Eu_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_Eu_alpha(:,tspan>=-500 & tspan < tspan_max)), [], 'all');

TFR_stats_plot(tspan, ze, power_NC_Eu_alpha - power_Go_Eu_alpha, ...
    squeeze(stats_NCvsGo_Eu_TFR_alpha.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- beta --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_Eu_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_Eu_beta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');

TFR_plot(tspan, ze, power_Go_Eu_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;
text(-450, 200, '\beta', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_Eu_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)

yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

ii = ii+1;
ax3 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_Eu_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_Eu_beta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');

TFR_stats_plot(tspan, ze, power_NC_Eu_beta - power_Go_Eu_beta, ...
    squeeze(stats_NCvsGo_Eu_TFR_beta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)

yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- gamma --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_Eu_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_Eu_gamma(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
%
TFR_plot(tspan, ze, power_Go_Eu_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
xticklabels({' ', ' ', ' '})
hold on;
text(-450, 200, '\gamma', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_Eu_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)

yticklabels(ytick_labels)
% hold on;
xticks([-500 0 500 1000])
%
% text(-450, 200, '\gamma', 'FontSize', 16, 'FontWeight','bold')
% xlabel('Time relative to saccade (ms)')

% c.Label.String = {' \muV^2'};
% colorbar_pos = c.Position;
% c.Position = [0.4543 colorbar_pos(2) colorbar_pos(3) colorbar_pos(4)];
% set(axi, 'Position', [axi.Position(1) axi.Position(2) 0.1079 axi.Position(4)])

ii = ii+1;
ax4 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_Eu_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_Eu_gamma(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
%
TFR_stats_plot(tspan, ze, power_NC_Eu_gamma - power_Go_Eu_gamma, ...
    squeeze(stats_NCvsGo_Eu_TFR_gamma.mask), caxis_lim_diff, z_depth, tspan_max)
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

%% ==================== X ====================

figure('Units', 'inches','Position',[0 0 6 7]);
tiledlayout(4, 3,'TileSpacing','Compact','Padding','Compact');

% -- theta --
ii = 1;
nexttile(ii);
caxis_lim = [-1 1].*max(abs([power_Go_X_theta(:,tspan>=-500 & tspan < tspan_max) ...
    power_NC_theta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
TFR_plot(tspan, ze, power_Go_X_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
text(-450, 200, '\theta', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_X_theta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})

ii = ii+1;
ax1 = nexttile(ii);
caxis_lim_diff = [-1 1].*max(abs(power_NC_X_theta(:,tspan>=-500 & tspan < tspan_max) ...
    - power_Go_X_theta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
TFR_stats_plot(tspan, ze, power_NC_X_theta - power_Go_X_theta, ...
    squeeze(stats_NCvsGo_X_TFR_theta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})

% -- alpha --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_X_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_X_alpha(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
TFR_plot(tspan, ze, power_Go_X_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
text(-450, 200, '\alpha', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_X_alpha, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})

ii = ii+1;
ax2 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_X_alpha(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_X_alpha(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
TFR_stats_plot(tspan, ze, power_NC_X_alpha - power_Go_X_alpha, ...
    squeeze(stats_NCvsGo_X_TFR_alpha.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})

% -- beta --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_X_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_X_beta(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
TFR_plot(tspan, ze, power_Go_X_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
text(-450, 200, '\beta', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
nexttile(ii);
TFR_plot(tspan, ze, power_NC_X_beta, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})

ii = ii+1;
ax3 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_X_beta(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_X_beta(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
TFR_stats_plot(tspan, ze, power_NC_X_beta - power_Go_X_beta, ...
    squeeze(stats_NCvsGo_X_TFR_beta.mask), caxis_lim_diff, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticklabels({' ', ' ', ' '})
hold on;

% -- gamma --
ii = ii+1;
nexttile(ii);
% caxis_lim = [-1 1].*max(abs([power_Go_X_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     power_NC_X_gamma(:,tspan>=-500 & tspan < tspan_max)]), [], 'all');
TFR_plot(tspan, ze, power_Go_X_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
xticklabels({' ', ' ', ' '})
text(-450, 200, '\gamma', 'FontSize', 16, 'FontWeight','bold')

ii = ii+1;
axi = nexttile(ii);
TFR_plot(tspan, ze, power_NC_X_gamma, caxis_lim, z_depth, tspan_max)
hold on;
xline(ERN, '--k', 'LineWidth',2)
xline(Pe, '-.b', 'LineWidth',2)
yticklabels(ytick_labels)
xticks([-500 0 500 1000])

ii = ii+1;
ax4 = nexttile(ii);
% caxis_lim = [-1 1].*max(abs(power_NC_X_gamma(:,tspan>=-500 & tspan < tspan_max) ...
%     - power_Go_X_gamma(:,tspan>=-500 & tspan < tspan_max)), [], 'all');
TFR_stats_plot(tspan, ze, power_NC_X_gamma - power_Go_X_gamma, ...
    squeeze(stats_NCvsGo_X_TFR_gamma.mask), caxis_lim_diff, z_depth, tspan_max)
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

% % [map_stats,~,~,~] = brewermap(256,'Accent');
% map_stats = zeros(size(bone)); %bone; % hot, zeros(size(bone))
% loc = 156;
% map_stats(loc:end,:) = ones(length(loc:size(map_stats,1)), 3);
% colormap(ax_stats1, (map_stats));
% colormap(ax_stats2, (map_stats));
% colormap(ax_stats3, (map_stats));
% colormap(ax_stats4, (map_stats));

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
