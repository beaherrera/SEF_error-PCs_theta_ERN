%% calculate the laminar power relative to saccade

clear
clc

noShowFigures = 1; % if =1; figures don't pop out
saveFigures = 1; % if =1; save figures (jpeg)

%% path 2 filtered lfps

path2data = fullfile('data', 'monkeys_data','laminar_data'); % IMPORTANT: you
% have to download and save the laminar data in this folder 
path2tfdata_target = fullfile(path2data, 'laminar_power_target'); % path 2
% laminar power relative to the target

%% perpendicular sessions

% sessions number
Sess_EuP1 = 14:19; % Eu, site P1
Sess_XP2P3 = [20:25 26:29]; % X, 20-25 site P2, 26-29 site P3
SessNumb = [Sess_EuP1, Sess_XP2P3];

% load electrodes' co-registration across sessions and monkeys
load(fullfile('data', 'monkeys_data','eleAlignment.mat'))

h = 150; % um, inter-electrode distance
Ne = 16; % number of electrodes
ze = 0:h:(Ne-1)*h;
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

%% frequencies of interest

thetaBand = [5 8]; % Hz [5 8]
alphaBand = [9 14]; % Hz [9 14]
betaBand = [15 29]; % Hz [15 29]
gammaBand = [30 80]; % Hz []

%% create saving directory

savePath = fullfile(path2data, 'laminar_power_saccade');
if ~exist(savePath, 'dir') % checks if the folder already exists
    mkdir(savePath);  % creates a folder named 'file'
end

savePathFigs = fullfile(savePath, 'Figs');
if ~exist(savePathFigs, 'dir') % checks if the folder already exists
    mkdir(savePathFigs);  % creates a folder named 'file'
end

%% load target laminar power data

load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_theta.mat'), 'target_theta_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_alpha.mat'), 'target_alpha_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_beta.mat'), 'target_beta_trials')
load(fullfile(path2tfdata_target, ...
    'sessions_amp_phase_gamma.mat'), 'target_gamma_trials')

%% pre-allocating memory

tspan = -500:1000;

theta_Go = cell(length(SessNumb), 2);
theta_NC = cell(length(SessNumb), 2);
alpha_Go = cell(length(SessNumb), 2);
alpha_NC = cell(length(SessNumb), 2);
beta_Go = cell(length(SessNumb), 2);
beta_NC = cell(length(SessNumb), 2);
gamma_Go = cell(length(SessNumb), 2);
gamma_NC = cell(length(SessNumb), 2);

theta_trials = cell(length(SessNumb), 2);
alpha_trials = cell(length(SessNumb), 2);
beta_trials = cell(length(SessNumb), 2);
gamma_trials = cell(length(SessNumb), 2);

%% loop for the sessions

for n = 1:length(SessNumb)

    sprintf('session %d', SessNumb(n))

    %% load filtered data
    load(fullfile(path2data,['Session_' num2str(SessNumb(n)) '_Saccade_filt.mat']),...
        'ft_lfpGo_theta_Saccade', 'ft_lfpNC_theta_Saccade', ...
        "ft_lfpNC_alpha_Saccade", "ft_lfpNC_beta_Saccade", ...
        "ft_lfpNC_gamma_Saccade", "ft_lfpGo_gamma_Saccade", ...
        "ft_lfpGo_beta_Saccade", "ft_lfpGo_alpha_Saccade")

    %% calculate the Hilbert Transforms

    lfpGo_theta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_theta_Saccade.trial, ...
        'UniformOutput',false);
    lfpGo_alpha_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_alpha_Saccade.trial, ...
        'UniformOutput',false);
    lfpGo_beta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_beta_Saccade.trial, ...
        'UniformOutput',false);
    lfpGo_gamma_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_gamma_Saccade.trial, ...
        'UniformOutput',false);

    lfpNC_theta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_theta_Saccade.trial, ...
        'UniformOutput',false);
    lfpNC_alpha_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_alpha_Saccade.trial, ...
        'UniformOutput',false);
    lfpNC_beta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_beta_Saccade.trial, ...
        'UniformOutput',false);
    lfpNC_gamma_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_gamma_Saccade.trial, ...
        'UniformOutput',false);

    %% calculate amplitude and phase envelope

    lfpGo_theta_amp = cellfun(@abs, lfpGo_theta_Hilbert, ...
        'UniformOutput',false);
    lfpGo_theta_phase = cellfun(@angle, lfpGo_theta_Hilbert, ...
        'UniformOutput',false);
    lfpNC_theta_amp = cellfun(@abs, lfpNC_theta_Hilbert, ...
        'UniformOutput',false);
    lfpNC_theta_phase = cellfun(@angle, lfpNC_theta_Hilbert, ...
        'UniformOutput',false);

    lfpGo_alpha_amp = cellfun(@abs, lfpGo_alpha_Hilbert, ...
        'UniformOutput',false);
    lfpGo_alpha_phase = cellfun(@angle, lfpGo_alpha_Hilbert, ...
        'UniformOutput',false);
    lfpNC_alpha_amp = cellfun(@abs, lfpNC_alpha_Hilbert, ...
        'UniformOutput',false);
    lfpNC_alpha_phase = cellfun(@angle, lfpNC_alpha_Hilbert, ...
        'UniformOutput',false);

    lfpGo_beta_amp = cellfun(@abs, lfpGo_beta_Hilbert, ...
        'UniformOutput',false);
    lfpGo_beta_phase = cellfun(@angle, lfpGo_beta_Hilbert, ...
        'UniformOutput',false);
    lfpNC_beta_amp = cellfun(@abs, lfpNC_beta_Hilbert, ...
        'UniformOutput',false);
    lfpNC_beta_phase = cellfun(@angle, lfpNC_beta_Hilbert, ...
        'UniformOutput',false);

    lfpGo_gamma_amp = cellfun(@abs, lfpGo_gamma_Hilbert, ...
        'UniformOutput',false);
    lfpGo_gamma_phase = cellfun(@angle, lfpGo_gamma_Hilbert, ...
        'UniformOutput',false);
    lfpNC_gamma_amp = cellfun(@abs, lfpNC_gamma_Hilbert, ...
        'UniformOutput',false);
    lfpNC_gamma_phase = cellfun(@angle, lfpNC_gamma_Hilbert, ...
        'UniformOutput',false);

    %% store session values

    theta_trials(n, :) = {lock2event(lfpGo_theta_amp, ft_lfpGo_theta_Saccade.time),...
        lock2event(lfpNC_theta_amp, ft_lfpNC_theta_Saccade.time)};
    alpha_trials(n, :) = {lock2event(lfpGo_alpha_amp, ft_lfpGo_alpha_Saccade.time),...
        lock2event(lfpNC_alpha_amp, ft_lfpNC_alpha_Saccade.time)};
    beta_trials(n, :) = {lock2event(lfpGo_beta_amp, ft_lfpGo_beta_Saccade.time),...
        lock2event(lfpNC_beta_amp, ft_lfpNC_beta_Saccade.time)};
    gamma_trials(n, :) = {lock2event(lfpGo_gamma_amp, ft_lfpGo_gamma_Saccade.time),...
        lock2event(lfpNC_gamma_amp, ft_lfpNC_gamma_Saccade.time)};

    theta_Go(n, :) = {mean(theta_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_theta_phase, ft_lfpGo_theta_Saccade.time), [], 3)};
    theta_NC(n, :) = {mean(theta_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_theta_phase, ft_lfpNC_theta_Saccade.time), [], 3)};

    alpha_Go(n, :) = {mean(alpha_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_alpha_phase, ft_lfpGo_alpha_Saccade.time), [], 3)};
    alpha_NC(n, :) = {mean(alpha_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_alpha_phase, ft_lfpNC_alpha_Saccade.time), [], 3)};

    beta_Go(n, :) = {mean(beta_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_beta_phase, ft_lfpGo_beta_Saccade.time), [], 3)};
    beta_NC(n, :) = {mean(beta_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_beta_phase, ft_lfpNC_beta_Saccade.time), [], 3)};

    gamma_Go(n, :) = {mean(gamma_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_gamma_phase, ft_lfpGo_gamma_Saccade.time), [], 3)};
    gamma_NC(n, :) = {mean(gamma_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_gamma_phase, ft_lfpNC_gamma_Saccade.time), [], 3)};

    %% plot and save figure

    min_amp_theta = min(abs([theta_Go{n, 1} theta_NC{n, 1}]).*1e6, [], 'all');
    min_amp_alpha = min(abs([alpha_Go{n, 1} alpha_NC{n, 1}]).*1e6, [], 'all');
    min_amp_beta = min(abs([beta_Go{n, 1} beta_NC{n, 1}]).*1e6, [], 'all');
    min_amp_gamma = min(abs([gamma_Go{n, 1} gamma_NC{n, 1}]).*1e6, [], 'all');

    max_amp_theta = max(abs([theta_Go{n, 1} theta_NC{n, 1}]).*1e6, [], 'all');
    max_amp_alpha = max(abs([alpha_Go{n, 1} alpha_NC{n, 1}]).*1e6, [], 'all');
    max_amp_beta = max(abs([beta_Go{n, 1} beta_NC{n, 1}]).*1e6, [], 'all');
    max_amp_gamma = max(abs([gamma_Go{n, 1} gamma_NC{n, 1}]).*1e6, [], 'all');

    %% Go
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
    end

    sgtitle('Go Trials')

    subplot(2,2,1)
    plot_amp_phase_depth(theta_Go{n, 1}.*1e6, theta_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')

    subplot(2,2,2)
    plot_amp_phase_depth(alpha_Go{n, 1}.*1e6, alpha_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')

    subplot(2,2,3)
    plot_amp_phase_depth(beta_Go{n, 1}.*1e6, beta_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')

    subplot(2,2,4)
    plot_amp_phase_depth(gamma_Go{n, 1}.*1e6, gamma_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_gamma max_amp_gamma])
    title('\gamma')

    if noShowFigures == 1
        % Set CreateFcn callback
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end

    % Save Fig file
    if saveFigures==1
        saveas(hFig, fullfile(savePathFigs, sprintf('Sess%d_Go.jpg', SessNumb(n))));
    end

    %% NC
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
    end

    sgtitle('NC Trials')

    tspan = -500:1000;

    subplot(2,2,1)
    plot_amp_phase_depth(theta_NC{n, 1}.*1e6, theta_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')

    subplot(2,2,2)
    plot_amp_phase_depth(alpha_NC{n, 1}.*1e6, alpha_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')

    subplot(2,2,3)
    plot_amp_phase_depth(beta_NC{n, 1}.*1e6, beta_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')

    subplot(2,2,4)
    plot_amp_phase_depth(gamma_NC{n, 1}.*1e6, gamma_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_gamma max_amp_gamma])
    title('\gamma')

    if noShowFigures == 1
        % Set CreateFcn callback
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end

    % Save Fig file
    if saveFigures==1
        saveas(hFig, fullfile(savePathFigs, sprintf('Sess%d_NC.jpg', SessNumb(n))));
    end

    close all

end

%% mod data
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

%% mean across sessions

mean_theta_amp_Go = mean(theta_Go_cell, 3, "omitnan");
mean_theta_phi_Go = circ_mean(cat(3, theta_Go{:,2}), [], 3);

mean_alpha_amp_Go = mean(alpha_Go_cell, 3, "omitnan");
mean_alpha_phi_Go = circ_mean(cat(3, alpha_Go{:,2}), [], 3);

mean_beta_amp_Go = mean(beta_Go_cell, 3, "omitnan");
mean_beta_phi_Go = circ_mean(cat(3, beta_Go{:,2}), [], 3);

mean_gamma_amp_Go = mean(gamma_Go_cell, 3, "omitnan");
mean_gamma_phi_Go = circ_mean(cat(3, gamma_Go{:,2}), [], 3);

mean_theta_amp_NC = mean(theta_NC_cell, 3, "omitnan");
mean_theta_phi_NC = circ_mean(cat(3, theta_NC{:,2}), [], 3);

mean_alpha_amp_NC = mean(alpha_NC_cell, 3, "omitnan");
mean_alpha_phi_NC = circ_mean(cat(3, alpha_NC{:,2}), [], 3);

mean_beta_amp_NC = mean(beta_NC_cell, 3, "omitnan");
mean_beta_phi_NC = circ_mean(cat(3, beta_NC{:,2}), [], 3);

mean_gamma_amp_NC = mean(gamma_NC_cell, 3, "omitnan");
mean_gamma_phi_NC = circ_mean(cat(3, gamma_NC{:,2}), [], 3);

%% stats 
max_amp_theta = max(abs([mean_theta_amp_Go mean_theta_amp_NC]), [], 'all');
max_amp_alpha = max(abs([mean_alpha_amp_Go mean_alpha_amp_NC]), [], 'all');
max_amp_beta = max(abs([mean_beta_amp_Go mean_beta_amp_NC]), [], 'all');
max_amp_gamma = max(abs([mean_gamma_amp_Go mean_gamma_amp_NC]), [], 'all');

min_amp_theta = -max_amp_theta; 
min_amp_alpha = -max_amp_alpha;
min_amp_beta = -max_amp_beta;
min_amp_gamma = -max_amp_gamma;

stats_cal_plots

%% store data

save(fullfile(savePath, 'sessions_amp_phase_theta.mat'), 'theta_NC', ...
    'theta_Go', 'mean_theta_amp_Go', 'mean_theta_amp_NC', ...
    'mean_theta_phi_Go', 'mean_theta_phi_NC', 'theta_trials')
save(fullfile(savePath, 'sessions_amp_phase_alpha.mat'), 'alpha_NC', ...
    'alpha_Go', 'mean_alpha_amp_Go', 'mean_alpha_amp_NC', ...
    'mean_alpha_phi_Go', 'mean_alpha_phi_NC', 'alpha_trials')
save(fullfile(savePath, 'sessions_amp_phase_beta.mat'), 'beta_NC', ...
    'beta_Go', 'mean_beta_amp_Go', 'mean_beta_amp_NC', ...
    'mean_beta_phi_Go', 'mean_beta_phi_NC', 'beta_trials')
save(fullfile(savePath, 'sessions_amp_phase_gamma.mat'), 'gamma_NC', ...
    'gamma_Go', 'mean_gamma_amp_Go', 'mean_gamma_amp_NC', ...
    'mean_gamma_phi_Go', 'mean_gamma_phi_NC', 'gamma_trials')

%% plot and save figure
%% Go
if noShowFigures == 1
    hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
else
    hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
end

sgtitle('Go Trials')

tspan = -500:1000;

subplot(2,2,1)
plot_amp_phase_depth(mean_theta_amp_Go, mean_theta_phi_Go, tspan, ze, ...
    z_depth, [min_amp_theta max_amp_theta])
title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_Go, mean_alpha_phi_Go, tspan, ze, ...
    z_depth, [min_amp_alpha max_amp_alpha])
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_Go, mean_beta_phi_Go, tspan, ze, ...
    z_depth, [min_amp_beta max_amp_beta])
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_Go, mean_gamma_phi_Go, tspan, ze, ...
    z_depth, [min_amp_gamma max_amp_gamma])
title('\gamma')

if noShowFigures == 1
    % Set CreateFcn callback
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
end

% Save Fig file
if saveFigures==1
    saveas(hFig, fullfile(savePathFigs, 'MeanSess_Go.jpg'));
end

%% NC
if noShowFigures == 1
    hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
else
    hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
end

sgtitle('NC Trials')

subplot(2,2,1)
plot_amp_phase_depth(mean_theta_amp_NC, mean_theta_phi_NC, tspan, ze, ...
    z_depth, [min_amp_theta max_amp_theta])
title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_NC, mean_alpha_phi_NC, tspan, ze, ...
    z_depth, [min_amp_alpha max_amp_alpha])
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_NC, mean_beta_phi_NC, tspan, ze, ...
    z_depth, [min_amp_beta max_amp_beta])
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_NC, mean_gamma_phi_NC, tspan, ze, ...
    z_depth, [min_amp_gamma max_amp_gamma])
title('\gamma')

if noShowFigures == 1
    % Set CreateFcn callback
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
end

% Save Fig file
if saveFigures==1
    saveas(hFig, fullfile(savePathFigs, 'MeanSess_NC.jpg'));
end

%% Error - Correct

if noShowFigures == 1
    hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
else
    hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
end

sgtitle('Error-Correct')

% min_amp_theta = min(mean_theta_amp_NC - mean_theta_amp_Go, [], 'all').*1e6;
% min_amp_alpha = min(mean_alpha_amp_NC - mean_alpha_amp_Go, [], 'all').*1e6;
% min_amp_beta = min(mean_beta_amp_NC - mean_beta_amp_Go, [], 'all').*1e6;
% min_amp_gamma = min(mean_gamma_amp_NC - mean_gamma_amp_Go, [], 'all').*1e6;

max_amp_theta = max(mean_theta_amp_NC - mean_theta_amp_Go, [], 'all');
max_amp_alpha = max(mean_alpha_amp_NC - mean_alpha_amp_Go, [], 'all');
max_amp_beta = max(mean_beta_amp_NC - mean_beta_amp_Go, [], 'all');
max_amp_gamma = max(mean_gamma_amp_NC - mean_gamma_amp_Go, [], 'all');

subplot(2,2,1)
plot_amp_phase_depth(mean_theta_amp_NC - mean_theta_amp_Go, ...
    mean_theta_phi_NC - mean_theta_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_theta max_amp_theta])
% colormap
[map,num,typ,scheme] = brewermap(256,'PiYG');
colormap((map));
title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_NC - mean_alpha_amp_Go, ...
    mean_alpha_phi_NC - mean_alpha_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_alpha max_amp_alpha])
% colormap
colormap((map));
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_NC - mean_beta_amp_Go, ...
    mean_beta_phi_NC - mean_beta_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_beta max_amp_beta])
% colormap
colormap((map));
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_NC - mean_gamma_amp_Go, ...
    mean_gamma_phi_NC - mean_gamma_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_gamma max_amp_gamma])
% colormap
colormap((map));
title('\gamma')

if noShowFigures == 1
    % Set CreateFcn callback
    set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
end

% Save Fig file
if saveFigures==1
    saveas(hFig, fullfile(savePathFigs, 'MeanSess_Difference.jpg'));
end
%%
% close all
