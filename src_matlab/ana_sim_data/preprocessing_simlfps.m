%% preprocessing of sim LFPs

clear
clc

%% path to data
% Simulated data for individual neurons was not included in the repository
% because of its large size.
save_folder = fullfile("data/sim_data/processed_lfps");
if ~exist(save_folder, 'dir') % checks if the folder already exists
    mkdir(save_folder);  % creates a folder named 'file'
end

path_simData_Go_L3 = fullfile('data', 'sim_data', 'Gotrials_625L3PCs');
path_simData_Go_L5 = fullfile('data', 'sim_data', 'Gotrials_1000L5PCs');
path_simData_NC_L3 = fullfile('data', 'sim_data', 'NCtrials_625L3PCs');
path_simData_NC_L5 = fullfile('data', 'sim_data', 'NCtrials_1000L5PCs');

file_name_L3 = 'Dend_r3.5_Apic_r2';
file_name_L5 = 'Dend_r2_Apic_r0.5';

%% simulation details

trials_number = 1:20;
num_trials = length(trials_number); % number of simulated trials
num_neuronsL3 = 625;
num_neuronsL5 = 1000; % number of simulated neurons
num_ele = 16; % number of electrodes in the simulated linear probe

warmup_period = 500e-3; % sec, warm up period. Discarted in analyses
simulation_length = 10000e-3; % sec, simulation length
dt = (2^(-4))*1e-3;  %[s] time increment
Fs = 1/(dt); % [Hz] sampling frequency
ts_all = 0:dt:(warmup_period + simulation_length);  %[s] time span,
% whole simulation
ts = 0:dt:simulation_length;  %[s] time span neglecting the warmup period

Fs_lfp = 1e3; % targeted LFP sampling rate
ts_out = 0:(Fs_lfp^-1):simulation_length;

%% load events time | ms
trials_number = 1:20;

load(fullfile('data', 'sim_data', 'events_timing_Go.mat'))
saccade_times_Go = double(saccade_times(1,trials_number))';
target_times_Go = double(target_times(1,trials_number))';

load(fullfile('data', 'sim_data', 'events_timing_NC.mat'))
saccade_times_NC = double(saccade_times(1,trials_number))';
target_times_NC = double(target_times(1,trials_number))';

%% get band-pass filters coefficients
Hd_theta = FilterW1_W2Hz(4,8);
Hd_alpha = FilterW1_W2Hz(9,14);
Hd_beta = FilterW1_W2Hz(15,29);
Hd_gamma = FilterW1_W2Hz(30,80);

%% 100 Hz low-pass filter: a two-pass fourth order Butterworth filter.
N = 4; % filter order
fsample = Fs_lfp; % lfp sampling rate
type = 'but'; % Butterworth filter
dir='twopass'; % two-pass filter
instabilityfix = 'reduce';
filterFreq = 100; % cutoff frequency

filterFunBut = @(x) ft_preproc_lowpassfilter(x, fsample, filterFreq, N,...
    type, dir, instabilityfix); % filter function

%% calculate the avg LFP across trials relative to the target and saccade onset

LFP_Go_trials = cell(num_ele, num_trials); % LFP Go trials
LFP_NC_trials = cell(num_ele, num_trials); % LFP NC trials

LFP_Go_L3_trials = cell(num_ele, num_trials); % LFP Go trials
LFP_NC_L3_trials = cell(num_ele, num_trials); % LFP NC trials

LFP_Go_L5_trials = cell(num_ele, num_trials); % LFP Go trials
LFP_NC_L5_trials = cell(num_ele, num_trials); % LFP NC trials

% loop for the trials
for r = trials_number
    
    loadPath2 = fullfile(path_simData_Go_L3, ['SimData_r' num2str(r) '_PS'...
        num2str(num_neuronsL3) '_' file_name_L3 '.mat']);
    load(loadPath2, 'LFP', 'ze') % loading extracellular potentials
    Ve_Go_L3 = LFP;
    ze_Go_L3 = ze;
    clearvars LFP ze
    
    loadPath2 = fullfile(path_simData_Go_L5, ['SimData_r' num2str(r) '_PS'...
        num2str(num_neuronsL5) '_' file_name_L5 '.mat']);
    load(loadPath2, 'LFP', 'ze') % loading extracellular potentials
    Ve_Go_L5 = LFP;
    ze_Go_L5 = ze;
    clearvars LFP ze
    
    loadPath2 = fullfile(path_simData_NC_L3, ['SimData_r' num2str(r) '_PS'...
        num2str(num_neuronsL3) '_' file_name_L3 '.mat']);
    load(loadPath2, 'LFP', 'ze') % loading extracellular potentials
    Ve_NC_L3 = LFP;
    ze_NC_L3 = ze;
    clearvars LFP ze
    
    loadPath2 = fullfile(path_simData_NC_L5, ['SimData_r' num2str(r) '_PS'...
        num2str(num_neuronsL5) '_' file_name_L5 '.mat']);
    load(loadPath2, 'LFP', 'ze') % loading extracellular potentials
    Ve_NC_L5 = LFP;
    ze_NC_L5 = ze;
    clearvars LFP ze
    
    if ~all((ze_Go_L3 == ze_Go_L5) & (ze_NC_L3 == ze_NC_L5))
        error("Channels coordinates are different in L3 and L5 sim.")
    end
    
    % -- remove warm-up period
    Ve_Go_L3 = Ve_Go_L3(:, ts_all>=warmup_period);
    Ve_Go_L5 = Ve_Go_L5(:, ts_all>=warmup_period);
    Ve_NC_L3 = Ve_NC_L3(:, ts_all>=warmup_period);
    Ve_NC_L5 = Ve_NC_L5(:, ts_all>=warmup_period);
    
    % -- get total LFP
    Ve_Go = scaling_L3.*Ve_Go_L3 + scaling_L5.*Ve_Go_L5;
    Ve_NC = scaling_L3.*Ve_NC_L3 + scaling_L5.*Ve_NC_L5;
    
    % -- downsampling the LFPs to Fs_lfp
    [VeD_Go_L3, ~] = process_resample('Compute', double(Ve_Go_L3), ts, Fs_lfp);
    [VeD_NC_L3, ~] = process_resample('Compute', double(Ve_NC_L3), ts, Fs_lfp);
    [VeD_Go_L5, ~] = process_resample('Compute', double(Ve_Go_L5), ts, Fs_lfp);
    [VeD_NC_L5, ~] = process_resample('Compute', double(Ve_NC_L5), ts, Fs_lfp);
    
    [VeD_Go, ~] = process_resample('Compute', double(Ve_Go), ts, Fs_lfp);
    [VeD_NC, ~] = process_resample('Compute', double(Ve_NC), ts, Fs_lfp);
    
    % -- store LFPs in a cell array
    LFP_Go_trials(:, (trials_number==r)) = num2cell(VeD_Go', 1)';
    LFP_NC_trials(:, (trials_number==r)) = num2cell(VeD_NC', 1)';
    LFP_Go_L3_trials(:, (trials_number==r)) = num2cell(VeD_Go_L3', 1)';
    LFP_NC_L3_trials(:, (trials_number==r)) = num2cell(VeD_NC_L3', 1)';
    LFP_Go_L5_trials(:, (trials_number==r)) = num2cell(VeD_Go_L5', 1)';
    LFP_NC_L5_trials(:, (trials_number==r)) = num2cell(VeD_NC_L5', 1)';
    
end

%% concatenate data as a continuous recording for filtering

cont_lfps = @(lfp_cell) arrayfun(@(elei) ...
    reshape(cell2mat(lfp_cell(elei, :)),[],1), ...
    1:size(lfp_cell, 1), 'UniformOutput', false);

LFP_Go_cont = cont_lfps(LFP_Go_trials);
LFP_NC_cont = cont_lfps(LFP_NC_trials);

LFP_Go_L3_cont = cont_lfps(LFP_Go_L3_trials);
LFP_NC_L3_cont = cont_lfps(LFP_NC_L3_trials);

LFP_Go_L5_cont = cont_lfps(LFP_Go_L5_trials);
LFP_NC_L5_cont = cont_lfps(LFP_NC_L5_trials);

%% filter at freq bands of interest
% theta
theta_filt = @(lfp_cell) cellfun(@(x) filtfilt(Hd_theta.Numerator,1,x), ...
    lfp_cell, 'UniformOutput', false);

LFP_Go_cont_theta = theta_filt(LFP_Go_cont);
LFP_NC_cont_theta = theta_filt(LFP_NC_cont);

LFP_Go_L3_cont_theta = theta_filt(LFP_Go_L3_cont);
LFP_NC_L3_cont_theta = theta_filt(LFP_NC_L3_cont);

LFP_Go_L5_cont_theta = theta_filt(LFP_Go_L5_cont);
LFP_NC_L5_cont_theta = theta_filt(LFP_NC_L5_cont);

% alpha
alpha_filt = @(lfp_cell) cellfun(@(x) filtfilt(Hd_alpha.Numerator,1,x), ...
    lfp_cell, 'UniformOutput', false);

LFP_Go_cont_alpha = alpha_filt(LFP_Go_cont);
LFP_NC_cont_alpha = alpha_filt(LFP_NC_cont);

LFP_Go_L3_cont_alpha = alpha_filt(LFP_Go_L3_cont);
LFP_NC_L3_cont_alpha = alpha_filt(LFP_NC_L3_cont);

LFP_Go_L5_cont_alpha = alpha_filt(LFP_Go_L5_cont);
LFP_NC_L5_cont_alpha = alpha_filt(LFP_NC_L5_cont);

% beta
beta_filt = @(lfp_cell) cellfun(@(x) filtfilt(Hd_beta.Numerator,1,x), ...
    lfp_cell, 'UniformOutput', false);

LFP_Go_cont_beta = beta_filt(LFP_Go_cont);
LFP_NC_cont_beta = beta_filt(LFP_NC_cont);

LFP_Go_L3_cont_beta = beta_filt(LFP_Go_L3_cont);
LFP_NC_L3_cont_beta = beta_filt(LFP_NC_L3_cont);

LFP_Go_L5_cont_beta = beta_filt(LFP_Go_L5_cont);
LFP_NC_L5_cont_beta = beta_filt(LFP_NC_L5_cont);

% gamma
gamma_filt = @(lfp_cell) cellfun(@(x) filtfilt(Hd_gamma.Numerator,1,x), ...
    lfp_cell, 'UniformOutput', false);

LFP_Go_cont_gamma = gamma_filt(LFP_Go_cont);
LFP_NC_cont_gamma = gamma_filt(LFP_NC_cont);

LFP_Go_L3_cont_gamma = gamma_filt(LFP_Go_L3_cont);
LFP_NC_L3_cont_gamma = gamma_filt(LFP_NC_L3_cont);

LFP_Go_L5_cont_gamma = gamma_filt(LFP_Go_L5_cont);
LFP_NC_L5_cont_gamma = gamma_filt(LFP_NC_L5_cont);

%% reshpae lfps into its original cell array shape

reshape_lfps = @(lfp_cell) arrayfun(@(elei) ...
    num2cell(reshape(lfp_cell{1, elei}, size(LFP_Go_trials{1,1}, 1), ...
    size(LFP_Go_trials, 2)), 1), ...
    (1:size(lfp_cell, 2))', 'UniformOutput', false);

% -- theta
LFP_Go_theta = reshape_lfps(LFP_Go_cont_theta);
LFP_Go_theta = cat(1, LFP_Go_theta{:});
LFP_NC_theta = reshape_lfps(LFP_NC_cont_theta);
LFP_NC_theta = cat(1, LFP_NC_theta{:});

LFP_Go_L3_theta = reshape_lfps(LFP_Go_L3_cont_theta);
LFP_Go_L3_theta = cat(1, LFP_Go_L3_theta{:});
LFP_NC_L3_theta = reshape_lfps(LFP_NC_L3_cont_theta);
LFP_NC_L3_theta = cat(1, LFP_NC_L3_theta{:});

LFP_Go_L5_theta = reshape_lfps(LFP_Go_L5_cont_theta);
LFP_Go_L5_theta = cat(1, LFP_Go_L5_theta{:});
LFP_NC_L5_theta = reshape_lfps(LFP_NC_L5_cont_theta);
LFP_NC_L5_theta = cat(1, LFP_NC_L5_theta{:});

% -- alpha
LFP_Go_alpha = reshape_lfps(LFP_Go_cont_alpha);
LFP_Go_alpha = cat(1, LFP_Go_alpha{:});
LFP_NC_alpha = reshape_lfps(LFP_NC_cont_alpha);
LFP_NC_alpha = cat(1, LFP_NC_alpha{:});

LFP_Go_L3_alpha = reshape_lfps(LFP_Go_L3_cont_alpha);
LFP_Go_L3_alpha = cat(1, LFP_Go_L3_alpha{:});
LFP_NC_L3_alpha = reshape_lfps(LFP_NC_L3_cont_alpha);
LFP_NC_L3_alpha = cat(1, LFP_NC_L3_alpha{:});

LFP_Go_L5_alpha = reshape_lfps(LFP_Go_L5_cont_alpha);
LFP_Go_L5_alpha = cat(1, LFP_Go_L5_alpha{:});
LFP_NC_L5_alpha = reshape_lfps(LFP_NC_L5_cont_alpha);
LFP_NC_L5_alpha = cat(1, LFP_NC_L5_alpha{:});

% -- beta
LFP_Go_beta = reshape_lfps(LFP_Go_cont_beta);
LFP_Go_beta = cat(1, LFP_Go_beta{:});
LFP_NC_beta = reshape_lfps(LFP_NC_cont_beta);
LFP_NC_beta = cat(1, LFP_NC_beta{:});

LFP_Go_L3_beta = reshape_lfps(LFP_Go_L3_cont_beta);
LFP_Go_L3_beta = cat(1, LFP_Go_L3_beta{:});
LFP_NC_L3_beta = reshape_lfps(LFP_NC_L3_cont_beta);
LFP_NC_L3_beta = cat(1, LFP_NC_L3_beta{:});

LFP_Go_L5_beta = reshape_lfps(LFP_Go_L5_cont_beta);
LFP_Go_L5_beta = cat(1, LFP_Go_L5_beta{:});
LFP_NC_L5_beta = reshape_lfps(LFP_NC_L5_cont_beta);
LFP_NC_L5_beta = cat(1, LFP_NC_L5_beta{:});

% -- gamma
LFP_Go_gamma = reshape_lfps(LFP_Go_cont_gamma);
LFP_Go_gamma = cat(1, LFP_Go_gamma{:});
LFP_NC_gamma = reshape_lfps(LFP_NC_cont_gamma);
LFP_NC_gamma = cat(1, LFP_NC_gamma{:});

LFP_Go_L3_gamma = reshape_lfps(LFP_Go_L3_cont_gamma);
LFP_Go_L3_gamma = cat(1, LFP_Go_L3_gamma{:});
LFP_NC_L3_gamma = reshape_lfps(LFP_NC_L3_cont_gamma);
LFP_NC_L3_gamma = cat(1, LFP_NC_L3_gamma{:});

LFP_Go_L5_gamma = reshape_lfps(LFP_Go_L5_cont_gamma);
LFP_Go_L5_gamma = cat(1, LFP_Go_L5_gamma{:});
LFP_NC_L5_gamma = reshape_lfps(LFP_NC_L5_cont_gamma);
LFP_NC_L5_gamma = cat(1, LFP_NC_L5_gamma{:});

%% convert sim trials to ft lfp structure

% -- raw sim LFPs
ft_lfpGo = convSimLFPs2ftRawDataStr(LFP_Go_trials, Fs_lfp, Fs);
ft_lfpNC = convSimLFPs2ftRawDataStr(LFP_NC_trials, Fs_lfp, Fs);

ft_lfpGo_L3 = convSimLFPs2ftRawDataStr(LFP_Go_L3_trials, Fs_lfp, Fs);
ft_lfpNC_L3 = convSimLFPs2ftRawDataStr(LFP_NC_L3_trials, Fs_lfp, Fs);

ft_lfpGo_L5 = convSimLFPs2ftRawDataStr(LFP_Go_L5_trials, Fs_lfp, Fs);
ft_lfpNC_L5 = convSimLFPs2ftRawDataStr(LFP_NC_L5_trials, Fs_lfp, Fs);

% ---- low-pass filter the lfp at 100Hz
ft_lfpGo.trial = cellfun(filterFunBut, ft_lfpGo.trial, ...
    'UniformOutput', false); % filtered data
ft_lfpNC.trial = cellfun(filterFunBut, ft_lfpNC.trial, ...
    'UniformOutput', false); % filtered data

ft_lfpGo_L3.trial = cellfun(filterFunBut, ft_lfpGo_L3.trial, ...
    'UniformOutput', false); % filtered data
ft_lfpNC_L3.trial = cellfun(filterFunBut, ft_lfpNC_L3.trial, ...
    'UniformOutput', false); % filtered data

ft_lfpGo_L5.trial = cellfun(filterFunBut, ft_lfpGo_L5.trial, ...
    'UniformOutput', false); % filtered data
ft_lfpNC_L5.trial = cellfun(filterFunBut, ft_lfpNC_L5.trial, ...
    'UniformOutput', false); % filtered data

% -- filtered LFPs
% --- theta
ft_lfpGo_theta = convSimLFPs2ftRawDataStr(LFP_Go_theta, Fs_lfp, Fs);
ft_lfpNC_theta = convSimLFPs2ftRawDataStr(LFP_NC_theta, Fs_lfp, Fs);

ft_lfpGo_L3_theta = convSimLFPs2ftRawDataStr(LFP_Go_L3_theta, Fs_lfp, Fs);
ft_lfpNC_L3_theta = convSimLFPs2ftRawDataStr(LFP_NC_L3_theta, Fs_lfp, Fs);

ft_lfpGo_L5_theta = convSimLFPs2ftRawDataStr(LFP_Go_L5_theta, Fs_lfp, Fs);
ft_lfpNC_L5_theta = convSimLFPs2ftRawDataStr(LFP_NC_L5_theta, Fs_lfp, Fs);

% --- alpha
ft_lfpGo_alpha = convSimLFPs2ftRawDataStr(LFP_Go_alpha, Fs_lfp, Fs);
ft_lfpNC_alpha = convSimLFPs2ftRawDataStr(LFP_NC_alpha, Fs_lfp, Fs);

ft_lfpGo_L3_alpha = convSimLFPs2ftRawDataStr(LFP_Go_L3_alpha, Fs_lfp, Fs);
ft_lfpNC_L3_alpha = convSimLFPs2ftRawDataStr(LFP_NC_L3_alpha, Fs_lfp, Fs);

ft_lfpGo_L5_alpha = convSimLFPs2ftRawDataStr(LFP_Go_L5_alpha, Fs_lfp, Fs);
ft_lfpNC_L5_alpha = convSimLFPs2ftRawDataStr(LFP_NC_L5_alpha, Fs_lfp, Fs);

% --- beta
ft_lfpGo_beta = convSimLFPs2ftRawDataStr(LFP_Go_beta, Fs_lfp, Fs);
ft_lfpNC_beta = convSimLFPs2ftRawDataStr(LFP_NC_beta, Fs_lfp, Fs);

ft_lfpGo_L3_beta = convSimLFPs2ftRawDataStr(LFP_Go_L3_beta, Fs_lfp, Fs);
ft_lfpNC_L3_beta = convSimLFPs2ftRawDataStr(LFP_NC_L3_beta, Fs_lfp, Fs);

ft_lfpGo_L5_beta = convSimLFPs2ftRawDataStr(LFP_Go_L5_beta, Fs_lfp, Fs);
ft_lfpNC_L5_beta = convSimLFPs2ftRawDataStr(LFP_NC_L5_beta, Fs_lfp, Fs);

% --- gamma
ft_lfpGo_gamma = convSimLFPs2ftRawDataStr(LFP_Go_gamma, Fs_lfp, Fs);
ft_lfpNC_gamma = convSimLFPs2ftRawDataStr(LFP_NC_gamma, Fs_lfp, Fs);

ft_lfpGo_L3_gamma = convSimLFPs2ftRawDataStr(LFP_Go_L3_gamma, Fs_lfp, Fs);
ft_lfpNC_L3_gamma = convSimLFPs2ftRawDataStr(LFP_NC_L3_gamma, Fs_lfp, Fs);

ft_lfpGo_L5_gamma = convSimLFPs2ftRawDataStr(LFP_Go_L5_gamma, Fs_lfp, Fs);
ft_lfpNC_L5_gamma = convSimLFPs2ftRawDataStr(LFP_NC_L5_gamma, Fs_lfp, Fs);

%% redefine trials relative to target and saccade
%% target

% -- biginning of the trials in samples
trialStartSampGo = ft_lfpGo.sampleinfo(:,1); % Go trials
trialStartSampNC = ft_lfpNC.sampleinfo(:,1); % Non-canceled trials

trialEndSampGo = ft_lfpGo.sampleinfo(:,2); % Go trials
trialEndSampNC = ft_lfpNC.sampleinfo(:,2); % Non-canceled trials

% -- target onset samples
targetCenterSampGo = target_times_Go + ...
    ft_lfpGo.sampleinfo(:,1); % Go trials
targetCenterSampNC = target_times_NC + ...
    ft_lfpNC.sampleinfo(:,1); % Non-canceled trials

% -- creating the trial definition structure
% Go trials
trlbeginGo = trialStartSampGo;
trlendGo   = trialEndSampGo;
offsetGo   = trialStartSampGo - targetCenterSampGo + 1;
trlGo = [trlbeginGo trlendGo offsetGo];
% Non-canceled trials
trlbeginNC = trialStartSampNC;
trlendNC   = trialEndSampNC;
offsetNC   = trialStartSampNC - targetCenterSampNC + 1;
trlNC = [trlbeginNC trlendNC offsetNC];

% -- redefining trials relative to target onset
% Go trials
cfg_Go = [];
cfg_Go.trl = trlGo;
ft_lfpGo_theta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_theta);
ft_lfpGo_L3_theta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_theta);
ft_lfpGo_L5_theta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_theta);

ft_lfpGo_alpha_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_alpha);
ft_lfpGo_L3_alpha_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_alpha);
ft_lfpGo_L5_alpha_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_alpha);

ft_lfpGo_beta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_beta);
ft_lfpGo_L3_beta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_beta);
ft_lfpGo_L5_beta_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_beta);

ft_lfpGo_gamma_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_gamma);
ft_lfpGo_L3_gamma_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_gamma);
ft_lfpGo_L5_gamma_Target = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_gamma);

% Non-canceled trials
cfg_NC = [];
cfg_NC.trl = trlNC;
ft_lfpNC_theta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_theta);
ft_lfpNC_L3_theta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_theta);
ft_lfpNC_L5_theta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_theta);

ft_lfpNC_alpha_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_alpha);
ft_lfpNC_L3_alpha_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_alpha);
ft_lfpNC_L5_alpha_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_alpha);

ft_lfpNC_beta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_beta);
ft_lfpNC_L3_beta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_beta);
ft_lfpNC_L5_beta_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_beta);

ft_lfpNC_gamma_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_gamma);
ft_lfpNC_L3_gamma_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_gamma);
ft_lfpNC_L5_gamma_Target = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_gamma);

%% saccade

% -- biginning of the trials in samples
trialStartSampGo = ft_lfpGo.sampleinfo(:,1); % Go trials
trialStartSampNC = ft_lfpNC.sampleinfo(:,1); % Non-canceled trials

trialEndSampGo = ft_lfpGo.sampleinfo(:,2); % Go trials
trialEndSampNC = ft_lfpNC.sampleinfo(:,2); % Non-canceled trials

% -- saccade onset samples
saccadeCenterSampGo = target_times_Go + saccade_times_Go + ...
    ft_lfpGo.sampleinfo(:,1); % Go trials
saccadeCenterSampNC = target_times_NC + saccade_times_NC + ...
    ft_lfpNC.sampleinfo(:,1); % Non-canceled trials

% -- creating the trial definition structure
% Go trials
trlbeginGo = trialStartSampGo;
trlendGo   = trialEndSampGo;
offsetGo   = trialStartSampGo - saccadeCenterSampGo + 1;
trlGo = [trlbeginGo trlendGo offsetGo];
% Non-canceled trials
trlbeginNC = trialStartSampNC;
trlendNC   = trialEndSampNC;
offsetNC   = trialStartSampNC - saccadeCenterSampNC + 1;
trlNC = [trlbeginNC trlendNC offsetNC];

% -- redefining trials relative to saccade onset
% Go trials
cfg_Go = [];
cfg_Go.trl = trlGo;
ft_lfpGo_theta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_theta);
ft_lfpGo_L3_theta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_theta);
ft_lfpGo_L5_theta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_theta);

ft_lfpGo_alpha_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_alpha);
ft_lfpGo_L3_alpha_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_alpha);
ft_lfpGo_L5_alpha_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_alpha);

ft_lfpGo_beta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_beta);
ft_lfpGo_L3_beta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_beta);
ft_lfpGo_L5_beta_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_beta);

ft_lfpGo_gamma_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_gamma);
ft_lfpGo_L3_gamma_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L3_gamma);
ft_lfpGo_L5_gamma_Saccade = ft_redefinetrial(cfg_Go,ft_lfpGo_L5_gamma);

% Non-canceled trials
cfg_NC = [];
cfg_NC.trl = trlNC;
ft_lfpNC_theta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_theta);
ft_lfpNC_L3_theta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_theta);
ft_lfpNC_L5_theta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_theta);

ft_lfpNC_alpha_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_alpha);
ft_lfpNC_L3_alpha_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_alpha);
ft_lfpNC_L5_alpha_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_alpha);

ft_lfpNC_beta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_beta);
ft_lfpNC_L3_beta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_beta);
ft_lfpNC_L5_beta_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_beta);

ft_lfpNC_gamma_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_gamma);
ft_lfpNC_L3_gamma_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L3_gamma);
ft_lfpNC_L5_gamma_Saccade = ft_redefinetrial(cfg_NC,ft_lfpNC_L5_gamma);

%% save data

save(fullfile(save_folder, 'sim_lfps.mat'), 'ft_lfpGo', ...
    'ft_lfpNC', 'ft_lfpGo_L3', 'ft_lfpNC_L3', 'ft_lfpGo_L5', ...
    'ft_lfpNC_L5')

save(fullfile(save_folder, 'filtered_sim_lfps.mat'), 'ft_lfpGo_theta', ...
    'ft_lfpNC_theta', 'ft_lfpGo_L3_theta', ...
    'ft_lfpNC_L3_theta', 'ft_lfpGo_L5_theta', ...
    'ft_lfpNC_L5_theta', 'ft_lfpGo_alpha', ...
    'ft_lfpNC_alpha', 'ft_lfpGo_L3_alpha', ...
    'ft_lfpNC_L3_alpha', 'ft_lfpGo_L5_alpha', ...
    "ft_lfpNC_L5_alpha", 'ft_lfpGo_beta', ...
    'ft_lfpNC_beta', 'ft_lfpGo_L3_beta', ...
    'ft_lfpNC_L3_beta', 'ft_lfpGo_L5_beta', ...
    'ft_lfpNC_L5_beta', 'ft_lfpGo_gamma', ...
    'ft_lfpNC_gamma', 'ft_lfpGo_L3_gamma', ...
    'ft_lfpNC_L3_gamma', 'ft_lfpGo_L5_gamma', ...
    'ft_lfpNC_L5_gamma')

save(fullfile(save_folder, 'filtered_sim_lfps_Target.mat'), 'ft_lfpGo_theta_Target', ...
    'ft_lfpNC_theta_Target', 'ft_lfpGo_L3_theta_Target', ...
    'ft_lfpNC_L3_theta_Target', 'ft_lfpGo_L5_theta_Target', ...
    'ft_lfpNC_L5_theta_Target', 'ft_lfpGo_alpha_Target', ...
    'ft_lfpNC_alpha_Target', 'ft_lfpGo_L3_alpha_Target', ...
    'ft_lfpNC_L3_alpha_Target', 'ft_lfpGo_L5_alpha_Target', ...
    "ft_lfpNC_L5_alpha_Target", 'ft_lfpGo_beta_Target', ...
    'ft_lfpNC_beta_Target', 'ft_lfpGo_L3_beta_Target', ...
    'ft_lfpNC_L3_beta_Target', 'ft_lfpGo_L5_beta_Target', ...
    'ft_lfpNC_L5_beta_Target', 'ft_lfpGo_gamma_Target', ...
    'ft_lfpNC_gamma_Target', 'ft_lfpGo_L3_gamma_Target', ...
    'ft_lfpNC_L3_gamma_Target', 'ft_lfpGo_L5_gamma_Target', ...
    'ft_lfpNC_L5_gamma_Target')

save(fullfile(save_folder, 'filtered_sim_lfps_Saccade.mat'), 'ft_lfpGo_theta_Saccade', ...
    'ft_lfpNC_theta_Saccade', 'ft_lfpGo_L3_theta_Saccade', ...
    'ft_lfpNC_L3_theta_Saccade', 'ft_lfpGo_L5_theta_Saccade', ...
    'ft_lfpNC_L5_theta_Saccade', 'ft_lfpGo_alpha_Saccade', ...
    'ft_lfpNC_alpha_Saccade', 'ft_lfpGo_L3_alpha_Saccade', ...
    'ft_lfpNC_L3_alpha_Saccade', 'ft_lfpGo_L5_alpha_Saccade', ...
    "ft_lfpNC_L5_alpha_Saccade", 'ft_lfpGo_beta_Saccade', ...
    'ft_lfpNC_beta_Saccade', 'ft_lfpGo_L3_beta_Saccade', ...
    'ft_lfpNC_L3_beta_Saccade', 'ft_lfpGo_L5_beta_Saccade', ...
    'ft_lfpNC_L5_beta_Saccade', 'ft_lfpGo_gamma_Saccade', ...
    'ft_lfpNC_gamma_Saccade', 'ft_lfpGo_L3_gamma_Saccade', ...
    'ft_lfpNC_L3_gamma_Saccade', 'ft_lfpGo_L5_gamma_Saccade', ...
    'ft_lfpNC_L5_gamma_Saccade')
