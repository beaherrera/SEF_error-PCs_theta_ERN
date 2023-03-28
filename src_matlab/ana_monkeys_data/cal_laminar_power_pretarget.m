%% calculate the laminar power relative to the target

clear
clc

noShowFigures = 1; % if =1; figures don't pop out
saveFigures = 1; % if =1; save figures (jpeg)

%% path 2 filtered lfps

path2data = fullfile('data', 'monkeys_data','laminar_data'); % IMPORTANT: you
% have to download and save the laminar data in this folder 

%% perpendicular sessions

% sessions number
Sess_EuP1 = 14:19; % Eu, site P1
Sess_XP2P3 = [20:25 26:29]; % X, 20-26 site P2, 26-29 site P3
SessNumb = [Sess_EuP1, Sess_XP2P3];

% load electrodes' co-registration across sessions and monkeys
load(fullfile('data', 'monkeys_data','eleAlignment.mat'))

h = 150; % um, inter-electrode distance
Ne = 16; % number of electrodes
ze = 0:h:(Ne-1)*h; % depths
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

%% frequencies of interest

thetaBand = [5 8]; % Hz
alphaBand = [9 14]; % Hz
betaBand = [15 29]; % Hz
gammaBand = [30 80]; % Hz 

%% create saving directory

savePath = fullfile('data', 'monkeys_data','laminar_power_target');
savePathFigs = fullfile(savePath, 'Figs');
if ~exist(savePathFigs, 'dir') % checks if the folder already exists
    mkdir(savePathFigs);  % creates a folder named 'file'
end

%% pre-allocating memory

target_theta_Go = cell(length(SessNumb), 2);
target_theta_NC = cell(length(SessNumb), 2);
target_alpha_Go = cell(length(SessNumb), 2);
target_alpha_NC = cell(length(SessNumb), 2);
target_beta_Go = cell(length(SessNumb), 2);
target_beta_NC = cell(length(SessNumb), 2);
target_gamma_Go = cell(length(SessNumb), 2);
target_gamma_NC = cell(length(SessNumb), 2);

target_theta_trials = cell(length(SessNumb), 2);
target_alpha_trials = cell(length(SessNumb), 2);
target_beta_trials = cell(length(SessNumb), 2);
target_gamma_trials = cell(length(SessNumb), 2);

%% loop for the sessions

for n = 1:length(SessNumb)

    sprintf('session %d', SessNumb(n))

    %% load filtered data
    load(fullfile(path2data,['Session_' num2str(SessNumb(n)) '_Target_filt.mat']),...
        'ft_lfpGo_theta_Target', 'ft_lfpNC_theta_Target', ...
        "ft_lfpNC_alpha_Target", "ft_lfpNC_beta_Target", ...
        "ft_lfpNC_gamma_Target", "ft_lfpGo_gamma_Target", ...
        "ft_lfpGo_beta_Target", "ft_lfpGo_alpha_Target")

    %% calculate the Hilbert Transforms

    lfpGo_theta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_theta_Target.trial, ...
        'UniformOutput',false);
    lfpGo_alpha_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_alpha_Target.trial, ...
        'UniformOutput',false);
    lfpGo_beta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_beta_Target.trial, ...
        'UniformOutput',false);
    lfpGo_gamma_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpGo_gamma_Target.trial, ...
        'UniformOutput',false);

    lfpNC_theta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_theta_Target.trial, ...
        'UniformOutput',false);
    lfpNC_alpha_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_alpha_Target.trial, ...
        'UniformOutput',false);
    lfpNC_beta_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_beta_Target.trial, ...
        'UniformOutput',false);
    lfpNC_gamma_Hilbert = cellfun(@(x) hilbert(x')', ft_lfpNC_gamma_Target.trial, ...
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

    target_theta_trials(n, :) = {lock2event(lfpGo_theta_amp, ft_lfpGo_theta_Target.time),...
        lock2event(lfpNC_theta_amp, ft_lfpNC_theta_Target.time)};
    target_alpha_trials(n, :) = {lock2event(lfpGo_alpha_amp, ft_lfpGo_alpha_Target.time),...
        lock2event(lfpNC_alpha_amp, ft_lfpNC_alpha_Target.time)};
    target_beta_trials(n, :) = {lock2event(lfpGo_beta_amp, ft_lfpGo_beta_Target.time),...
        lock2event(lfpNC_beta_amp, ft_lfpNC_beta_Target.time)};
    target_gamma_trials(n, :) = {lock2event(lfpGo_gamma_amp, ft_lfpGo_gamma_Target.time),...
        lock2event(lfpNC_gamma_amp, ft_lfpNC_gamma_Target.time)};

    target_theta_Go(n, :) = {mean(target_theta_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_theta_phase, ft_lfpGo_theta_Target.time), [], 3)};
    target_theta_NC(n, :) = {mean(target_theta_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_theta_phase, ft_lfpNC_theta_Target.time), [], 3)};

    target_alpha_Go(n, :) = {mean(target_alpha_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_alpha_phase, ft_lfpGo_alpha_Target.time), [], 3)};
    target_alpha_NC(n, :) = {mean(target_alpha_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_alpha_phase, ft_lfpNC_alpha_Target.time), [], 3)};

    target_beta_Go(n, :) = {mean(target_beta_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_beta_phase, ft_lfpGo_beta_Target.time), [], 3)};
    target_beta_NC(n, :) = {mean(target_beta_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_beta_phase, ft_lfpNC_beta_Target.time), [], 3)};

    target_gamma_Go(n, :) = {mean(target_gamma_trials{n, 1},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpGo_gamma_phase, ft_lfpGo_gamma_Target.time), [], 3)};
    target_gamma_NC(n, :) = {mean(target_gamma_trials{n, 2},...
        3, "omitnan"),...
        circ_mean(lock2event(lfpNC_gamma_phase, ft_lfpNC_gamma_Target.time), [], 3)};

    %% plot and save figure

    min_amp_theta = min(abs([target_theta_Go{n, 1} target_theta_NC{n, 1}]).*1e6, [], 'all');
    min_amp_alpha = min(abs([target_alpha_Go{n, 1} target_alpha_NC{n, 1}]).*1e6, [], 'all');
    min_amp_beta = min(abs([target_beta_Go{n, 1} target_beta_NC{n, 1}]).*1e6, [], 'all');
    min_amp_gamma = min(abs([target_gamma_Go{n, 1} target_gamma_NC{n, 1}]).*1e6, [], 'all');

    max_amp_theta = max(abs([target_theta_Go{n, 1} target_theta_NC{n, 1}]).*1e6, [], 'all');
    max_amp_alpha = max(abs([target_alpha_Go{n, 1} target_alpha_NC{n, 1}]).*1e6, [], 'all');
    max_amp_beta = max(abs([target_beta_Go{n, 1} target_beta_NC{n, 1}]).*1e6, [], 'all');
    max_amp_gamma = max(abs([target_gamma_Go{n, 1} target_gamma_NC{n, 1}]).*1e6, [], 'all');

    %% Go
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
    end

    sgtitle('Go Trials')

    tspan = -500:1000;

    subplot(2,2,1)
    plot_amp_phase_depth(target_theta_Go{n, 1}.*1e6, target_theta_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')

    subplot(2,2,2)
    plot_amp_phase_depth(target_alpha_Go{n, 1}.*1e6, target_alpha_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')

    subplot(2,2,3)
    plot_amp_phase_depth(target_beta_Go{n, 1}.*1e6, target_beta_Go{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')

    subplot(2,2,4)
    plot_amp_phase_depth(target_gamma_Go{n, 1}.*1e6, target_gamma_Go{n, 2}.*1e6, tspan, ze, ...
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
    plot_amp_phase_depth(target_theta_NC{n, 1}.*1e6, target_theta_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')

    subplot(2,2,2)
    plot_amp_phase_depth(target_alpha_NC{n, 1}.*1e6, target_alpha_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')

    subplot(2,2,3)
    plot_amp_phase_depth(target_beta_NC{n, 1}.*1e6, target_beta_NC{n, 2}.*1e6, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')

    subplot(2,2,4)
    plot_amp_phase_depth(target_gamma_NC{n, 1}.*1e6, target_gamma_NC{n, 2}.*1e6, tspan, ze, ...
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

%% mean across sessions

mean_theta_amp_Go = mean(cat(3, target_theta_Go{:,1}), 3, "omitnan");
mean_theta_phi_Go = circ_mean(cat(3, target_theta_Go{:,2}), [], 3);

mean_alpha_amp_Go = mean(cat(3, target_alpha_Go{:,1}), 3, "omitnan");
mean_alpha_phi_Go = circ_mean(cat(3, target_alpha_Go{:,2}), [], 3);

mean_beta_amp_Go = mean(cat(3, target_beta_Go{:,1}), 3, "omitnan");
mean_beta_phi_Go = circ_mean(cat(3, target_beta_Go{:,2}), [], 3);

mean_gamma_amp_Go = mean(cat(3, target_gamma_Go{:,1}), 3, "omitnan");
mean_gamma_phi_Go = circ_mean(cat(3, target_gamma_Go{:,2}), [], 3);

mean_theta_amp_NC = mean(cat(3, target_theta_NC{:,1}), 3, "omitnan");
mean_theta_phi_NC = circ_mean(cat(3, target_theta_NC{:,2}), [], 3);

mean_alpha_amp_NC = mean(cat(3, target_alpha_NC{:,1}), 3, "omitnan");
mean_alpha_phi_NC = circ_mean(cat(3, target_alpha_NC{:,2}), [], 3);

mean_beta_amp_NC = mean(cat(3, target_beta_NC{:,1}), 3, "omitnan");
mean_beta_phi_NC = circ_mean(cat(3, target_beta_NC{:,2}), [], 3);

mean_gamma_amp_NC = mean(cat(3, target_gamma_NC{:,1}), 3, "omitnan");
mean_gamma_phi_NC = circ_mean(cat(3, target_gamma_NC{:,2}), [], 3);

%% store data

save(fullfile(savePath, 'sessions_amp_phase_theta.mat'), 'target_theta_NC', ...
    'target_theta_Go', 'mean_theta_amp_Go', 'mean_theta_amp_NC', ...
    'mean_theta_phi_Go', 'mean_theta_phi_NC', 'target_theta_trials')
save(fullfile(savePath, 'sessions_amp_phase_alpha.mat'), 'target_alpha_NC', ...
    'target_alpha_Go', 'mean_alpha_amp_Go', 'mean_alpha_amp_NC', ...
    'mean_alpha_phi_Go', 'mean_alpha_phi_NC', 'target_alpha_trials')
save(fullfile(savePath, 'sessions_amp_phase_beta.mat'), 'target_beta_NC', ...
    'target_beta_Go', 'mean_beta_amp_Go', 'mean_beta_amp_NC', ...
    'mean_beta_phi_Go', 'mean_beta_phi_NC', 'target_beta_trials')
save(fullfile(savePath, 'sessions_amp_phase_gamma.mat'), 'target_gamma_NC', ...
    'target_gamma_Go', 'mean_gamma_amp_Go', 'mean_gamma_amp_NC', ...
    'mean_gamma_phi_Go', 'mean_gamma_phi_NC', 'target_gamma_trials')

%% plot and save figure
min_amp_theta = min(abs([mean_theta_amp_Go mean_theta_amp_NC]).*1e6, [], 'all');
min_amp_alpha = min(abs([mean_alpha_amp_Go mean_alpha_amp_NC]).*1e6, [], 'all');
min_amp_beta = min(abs([mean_beta_amp_Go mean_beta_amp_NC]).*1e6, [], 'all');
min_amp_gamma = min(abs([mean_gamma_amp_Go mean_gamma_amp_NC]).*1e6, [], 'all');

max_amp_theta = max(abs([mean_theta_amp_Go mean_theta_amp_NC]).*1e6, [], 'all');
max_amp_alpha = max(abs([mean_alpha_amp_Go mean_alpha_amp_NC]).*1e6, [], 'all');
max_amp_beta = max(abs([mean_beta_amp_Go mean_beta_amp_NC]).*1e6, [], 'all');
max_amp_gamma = max(abs([mean_gamma_amp_Go mean_gamma_amp_NC]).*1e6, [], 'all');

%% Go
if noShowFigures == 1
    hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
else
    hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
end

sgtitle('Go Trials')

tspan = -500:1000;

subplot(2,2,1)
plot_amp_phase_depth(mean_theta_amp_Go.*1e6, mean_theta_phi_Go, tspan, ze, ...
    z_depth, [min_amp_theta max_amp_theta])
hold on;

title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_Go.*1e6, mean_alpha_phi_Go, tspan, ze, ...
    z_depth, [min_amp_alpha max_amp_alpha])
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_Go.*1e6, mean_beta_phi_Go, tspan, ze, ...
    z_depth, [min_amp_beta max_amp_beta])
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_Go.*1e6, mean_gamma_phi_Go, tspan, ze, ...
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
plot_amp_phase_depth(mean_theta_amp_NC.*1e6, mean_theta_phi_NC, tspan, ze, ...
    z_depth, [min_amp_theta max_amp_theta])
title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_NC.*1e6, mean_alpha_phi_NC, tspan, ze, ...
    z_depth, [min_amp_alpha max_amp_alpha])
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_NC.*1e6, mean_beta_phi_NC, tspan, ze, ...
    z_depth, [min_amp_beta max_amp_beta])
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_NC.*1e6, mean_gamma_phi_NC, tspan, ze, ...
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

max_amp_theta = max(mean_theta_amp_NC - mean_theta_amp_Go, [], 'all').*1e6;
max_amp_alpha = max(mean_alpha_amp_NC - mean_alpha_amp_Go, [], 'all').*1e6;
max_amp_beta = max(mean_beta_amp_NC - mean_beta_amp_Go, [], 'all').*1e6;
max_amp_gamma = max(mean_gamma_amp_NC - mean_gamma_amp_Go, [], 'all').*1e6;

subplot(2,2,1)
plot_amp_phase_depth(mean_theta_amp_NC.*1e6 - mean_theta_amp_Go.*1e6, ...
    mean_theta_phi_NC.*1e6 - mean_theta_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_theta max_amp_theta])
% colormap
[map,num,typ,scheme] = brewermap(256,'PiYG');
colormap((map));
title('\theta')

subplot(2,2,2)
plot_amp_phase_depth(mean_alpha_amp_NC.*1e6 - mean_alpha_amp_Go.*1e6, ...
    mean_alpha_phi_NC.*1e6 - mean_alpha_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_alpha max_amp_alpha])
% colormap
colormap((map));
title('\alpha')

subplot(2,2,3)
plot_amp_phase_depth(mean_beta_amp_NC.*1e6 - mean_beta_amp_Go.*1e6, ...
    mean_beta_phi_NC.*1e6 - mean_beta_phi_Go, tspan, ze, ...
    z_depth, [-max_amp_beta max_amp_beta])
% colormap
colormap((map));
title('\beta')

subplot(2,2,4)
plot_amp_phase_depth(mean_gamma_amp_NC.*1e6 - mean_gamma_amp_Go.*1e6, ...
    mean_gamma_phi_NC.*1e6 - mean_gamma_phi_Go, tspan, ze, ...
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
close all
