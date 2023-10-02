%% calculate the laminar power relative to saccade
% -- no baseline correction is performed to create the figs --

clear
clc

noShowFigures = 0; % if =1; figures don't pop out
saveFigures = 0; % if =1; save figures (jpeg)

%% path 2 filtered lfps

path2data = fullfile('data', 'monkeys_data', 'laminar_data');
% IMPORTANT: you have to download and save the laminar data in this folder

path2tfdata_target = fullfile('data', 'monkeys_data', ...
    'laminar_power_target'); % path 2
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

savePath = fullfile('data', 'monkeys_data', 'laminar_power_saccade');
if ~exist(savePath, 'dir') % checks if the folder already exists
    mkdir(savePath);  % creates a folder named 'file'
end

savePathFigs = fullfile(savePath, 'Figs');
if ~exist(savePathFigs, 'dir') % checks if the folder already exists
    mkdir(savePathFigs);  % creates a folder named 'file'
end

%% pre-allocating memory

tspan = -500:1000;

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
    lfpNC_theta_amp = cellfun(@abs, lfpNC_theta_Hilbert, ...
        'UniformOutput',false);
    
    lfpGo_alpha_amp = cellfun(@abs, lfpGo_alpha_Hilbert, ...
        'UniformOutput',false);
    lfpNC_alpha_amp = cellfun(@abs, lfpNC_alpha_Hilbert, ...
        'UniformOutput',false);
    
    lfpGo_beta_amp = cellfun(@abs, lfpGo_beta_Hilbert, ...
        'UniformOutput',false);
    lfpNC_beta_amp = cellfun(@abs, lfpNC_beta_Hilbert, ...
        'UniformOutput',false);
    
    lfpGo_gamma_amp = cellfun(@abs, lfpGo_gamma_Hilbert, ...
        'UniformOutput',false);
    lfpNC_gamma_amp = cellfun(@abs, lfpNC_gamma_Hilbert, ...
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
    
    theta_Go = mean((theta_trials{n, 1}.*1e6).^2,...
        3, "omitnan");
    theta_NC = mean((theta_trials{n, 2}.*1e6).^2,...
        3, "omitnan");
    
    alpha_Go = mean((alpha_trials{n, 1}.*1e6).^2,...
        3, "omitnan");
    alpha_NC = mean((alpha_trials{n, 2}.*1e6).^2,...
        3, "omitnan");
    
    beta_Go = mean((beta_trials{n, 1}.*1e6).^2,...
        3, "omitnan");
    beta_NC = mean((beta_trials{n, 2}.*1e6).^2,...
        3, "omitnan");
    
    gamma_Go = mean((gamma_trials{n, 1}.*1e6).^2,...
        3, "omitnan");
    gamma_NC = mean((gamma_trials{n, 2}.*1e6).^2,...
        3, "omitnan");
    
    %% plot and save figure
    
    min_amp_theta = min(abs([theta_Go theta_NC]), [], 'all');
    min_amp_alpha = min(abs([alpha_Go alpha_NC]), [], 'all');
    min_amp_beta = min(abs([beta_Go beta_NC]), [], 'all');
    min_amp_gamma = min(abs([gamma_Go gamma_NC]), [], 'all');
    
    max_amp_theta = max(abs([theta_Go theta_NC]), [], 'all');
    max_amp_alpha = max(abs([alpha_Go alpha_NC]), [], 'all');
    max_amp_beta = max(abs([beta_Go beta_NC]), [], 'all');
    max_amp_gamma = max(abs([gamma_Go gamma_NC]), [], 'all');
    
    %% Go
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
    end
    
    sgtitle('Go Trials')
    
    subplot(2,2,1)
    plot_power_depth(theta_Go, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')
    
    subplot(2,2,2)
    plot_power_depth(alpha_Go, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')
    
    subplot(2,2,3)
    plot_power_depth(beta_Go, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')
    
    subplot(2,2,4)
    plot_power_depth(gamma_Go, tspan, ze, ...
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
    plot_power_depth(theta_NC, tspan, ze, ...
        z_depth, [min_amp_theta max_amp_theta])
    title('\theta')
    
    subplot(2,2,2)
    plot_power_depth(alpha_NC, tspan, ze, ...
        z_depth, [min_amp_alpha max_amp_alpha])
    title('\alpha')
    
    subplot(2,2,3)
    plot_power_depth(beta_NC, tspan, ze, ...
        z_depth, [min_amp_beta max_amp_beta])
    title('\beta')
    
    subplot(2,2,4)
    plot_power_depth(gamma_NC, tspan, ze, ...
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
    
    %     close all
    
end

%% store data

save(fullfile(savePath, 'sessions_theta.mat'), 'theta_trials')
save(fullfile(savePath, 'sessions_alpha.mat'), 'alpha_trials')
save(fullfile(savePath, 'sessions_beta.mat'), 'beta_trials')
save(fullfile(savePath, 'sessions_gamma.mat'), 'gamma_trials')

%%

close all

%% ================== sub-functions ==================
function plot_power_depth(amp, tspan, ze, z_depth, caxis_lim)
%PLOT_POWER_DEPTH Plot Amp & Phi across depth and time
%   Plots the laminar maps of the amplitude envelope with the phase
%   superimposed.
% Inputs:
%

if any(isnan(amp(:,1)))
    
    amp = amp(~isnan(amp(:,1)), :);
    ze = ze(~isnan(amp(:,1)));
    z_depth = z_depth([1:8 10:17]);
    z_depth = z_depth(~isnan(amp(:,1)));
    if any(ismember([75 -75], z_depth))
        z_depth = sort([z_depth 0], 'descend');
    end
end

q = length(tspan); % number of interpolating points
[zeq, tspanq] = meshgrid(linspace(ze(1),ze(end), q), tspan);

Amp_interp = interp2(tspan, ze, amp, tspanq, zeq, 'spline');

imagesc(tspan, zeq(1,:), Amp_interp');
caxis(caxis_lim)
c = colorbar;
colormap(jet);
c.Label.FontSize = 9;
ax = gca; % current axes
ax.FontSize = 9;
ylabel('z (mm)','FontSize',9);
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
set(gca,'linewidth',1.5,'fontsize',9,'fontweight','bold','TickDir','out')

end
