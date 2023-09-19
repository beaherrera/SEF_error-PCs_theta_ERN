%% calculate the laminar power relative to the target

clear
clc

%% path 2 filtered lfps

path2data = fullfile('data', 'monkeys_data', 'laminar_data');
% IMPORTANT: you have to download and save the laminar data in this folder

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

%% pre-allocating memory

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
    
    target_theta_trials(n, :) = {lock2event(lfpGo_theta_amp, ft_lfpGo_theta_Target.time),...
        lock2event(lfpNC_theta_amp, ft_lfpNC_theta_Target.time)};
    target_alpha_trials(n, :) = {lock2event(lfpGo_alpha_amp, ft_lfpGo_alpha_Target.time),...
        lock2event(lfpNC_alpha_amp, ft_lfpNC_alpha_Target.time)};
    target_beta_trials(n, :) = {lock2event(lfpGo_beta_amp, ft_lfpGo_beta_Target.time),...
        lock2event(lfpNC_beta_amp, ft_lfpNC_beta_Target.time)};
    target_gamma_trials(n, :) = {lock2event(lfpGo_gamma_amp, ft_lfpGo_gamma_Target.time),...
        lock2event(lfpNC_gamma_amp, ft_lfpNC_gamma_Target.time)};
    
end

%% store data

save(fullfile(savePath, 'sessions_amp_phase_theta.mat'), 'target_theta_trials')
save(fullfile(savePath, 'sessions_amp_phase_alpha.mat'), 'target_alpha_trials')
save(fullfile(savePath, 'sessions_amp_phase_beta.mat'), 'target_beta_trials')
save(fullfile(savePath, 'sessions_amp_phase_gamma.mat'), 'target_gamma_trials')
