%% calculate the ISI dist for the sim L3 neurons and its density maps

clear
clc

%% path to data
% Simulated data for individual neurons was not included in the repository
% because of its large size.
path_simData_Go = ['D:\Theta_paper_sim\results_L3PCsPopMky\' ...
    'Go_trial_Sept9_mpi\neurons#10_clustered_synp\' ...
    'StimDend#4_StimApic#4'];

path_simData_NC = ['D:\Theta_paper_sim\results_L3PCsPopMky\' ...
    'NC_trial_Sept9_mpi\neurons#10_clustered_synp\' ...
    'StimDend#4_StimApic#4'];

file_name_NC = 'Dend_r3.5_Apic_r2';
file_name_Go = 'Dend_r3.5_Apic_r2';

%% load events time

load(fullfile('data', 'sim_data', 'events_timing_Go.mat'))
saccade_times_Go = saccade_times;
target_times_Go = target_times;

load(fullfile('data', 'sim_data', 'events_timing_NC.mat'))
saccade_times_NC = saccade_times;
target_times_NC = target_times;

load("data\sim_data\tone_times.mat")

%% simulation details
target_times_Go = double(target_times_Go(1,:)).*1e-3; % sec, target times
saccade_times_Go = double(saccade_times_Go(1,:)).*1e-3; % sec,
% saccade times relative to target
target_times_NC = double(target_times_NC(1,:)).*1e-3; % sec, target times
saccade_times_NC = double(saccade_times_NC(1,:)).*1e-3; % sec,
% saccade times relative to target
tone_times_Go = double(tone_times_Go(1,:)).*1e-3; % sec, tone times
% relative to saccade
tone_times_NC = double(tone_times_NC(1,:)).*1e-3; % sec, tone times
% relative to saccade

num_trials = 116; % number of simulated trials
num_neurons = 10; % number of simulated neurons
warmup_period = 500e-3; % sec
simulation_length = 10000e-3; % simulation length
dt = (2^(-4))*1e-3;  %[s] time increment
Fs = 1/(dt); % [Hz] sampling frequency

%% isi bins
max_isi = 0.1; % sec, max value for isi
bins = 0:0.002:max_isi;
bins_middle = bins(1:end-1) + diff(bins)./2;
bins_middle = bins_middle.*1e3;

%% ========== creat field-trip spike structure ==========
%% define trials relative to the event of interest
%% --- relative to target ---
%% ---- Go trials ----
clearvars spikes_Go_target
spikes_Go_target.label = cell(1, num_neurons); % name of the units
spikes_Go_target.timestamp = cell(1, num_neurons); % spike times in recording samples
spikes_Go_target.trial = cell(1, num_neurons);
spikes_Go_target.time = cell(1, num_neurons);
spikes_Go_target.trialtime = zeros(num_trials, 2);
spikes_Go_target.trialtime(:,1) = spikes_Go_target.trialtime(:,1) - target_times_Go';
spikes_Go_target.trialtime(:,2) = simulation_length - target_times_Go';
spikes_Go_target.timestampdimord = '{chan}_spike'; % dimension of the spiking data

for n=0:(num_neurons-1)
    
    spikes_Go_target.label{n+1} = ['neuron_' num2str(n)];
    
    for t=1:num_trials
        
        load(fullfile(path_simData_Go, ['NeuronsData_r' num2str(t) '_n#' ...
            '' num2str(n) '_' file_name_Go '.mat']), ...
            'soma_spkTimes')
        soma_spkTimesT = soma_spkTimes - target_times_Go(1,t);
        soma_spkTimesT_tone = soma_spkTimesT(soma_spkTimesT < tone_times_Go(1,t));
        spikes_Go_target.timestamp{n+1} = cat(2, spikes_Go_target.timestamp{n+1}, ...
            (soma_spkTimesT_tone).*Fs+1);
        spikes_Go_target.time{n+1} = cat(2, spikes_Go_target.time{n+1}, ...
            (soma_spkTimesT_tone));
        spikes_Go_target.trial{n+1} = cat(2, spikes_Go_target.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimesT_tone)]));
    end
end

%% ---- NC trials ----
clearvars spikes_NC_target soma_spkTimesT soma_spkTimesT_tone
spikes_NC_target.label = cell(1, num_neurons); % name of the units
spikes_NC_target.timestamp = cell(1, num_neurons); % spike times in recording samples
spikes_NC_target.trial = cell(1, num_neurons);
spikes_NC_target.time = cell(1, num_neurons);
spikes_NC_target.trialtime = zeros(num_trials, 2);
spikes_NC_target.trialtime(:,1) = spikes_NC_target.trialtime(:,1) - target_times_NC';
spikes_NC_target.trialtime(:,2) = simulation_length - target_times_NC';
spikes_NC_target.timestampdimord = '{chan}_spike'; % dimension of the spiking data

for n=0:(num_neurons-1)
    
    spikes_NC_target.label{n+1} = ['neuron_' num2str(n)];
    
    for t=1:num_trials
        
        load(fullfile(path_simData_NC, ['NeuronsData_r' num2str(t) '_n#' ...
            '' num2str(n) '_' file_name_NC '.mat']), ...
            'soma_spkTimes')
        soma_spkTimesT = soma_spkTimes - target_times_NC(1,t);
        soma_spkTimesT_tone = soma_spkTimesT(soma_spkTimesT < tone_times_NC(1,t));
        spikes_NC_target.timestamp{n+1} = cat(2, spikes_NC_target.timestamp{n+1}, ...
            (soma_spkTimesT_tone).*Fs+1);
        spikes_NC_target.time{n+1} = cat(2, spikes_NC_target.time{n+1}, ...
            (soma_spkTimesT_tone));
        spikes_NC_target.trial{n+1} = cat(2, spikes_NC_target.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimesT_tone)]));
    end
end

%% --- relative to saccade ---
%% ---- Go trials ----
clearvars spikes_Go_saccade soma_spkTimesT soma_spkTimesT_tone
spikes_Go_saccade.label = cell(1, num_neurons); % name of the units
spikes_Go_saccade.timestamp = cell(1, num_neurons); % spike times in recording samples
spikes_Go_saccade.trial = cell(1, num_neurons);
spikes_Go_saccade.time = cell(1, num_neurons);
spikes_Go_saccade.trialtime = zeros(num_trials, 2);
spikes_Go_saccade.trialtime(:,1) = spikes_Go_saccade.trialtime(:,1) - (target_times_Go ...
    + saccade_times_Go)';
spikes_Go_saccade.trialtime(:,2) = simulation_length - (target_times_Go ...
    + saccade_times_Go)';
spikes_Go_saccade.timestampdimord = '{chan}_spike'; % dimension of the spiking data

for n=0:(num_neurons-1)
    
    spikes_Go_saccade.label{n+1} = ['neuron_' num2str(n)];
    
    for t=1:num_trials
        
        load(fullfile(path_simData_Go, ['NeuronsData_r' num2str(t) '_n#' ...
            '' num2str(n) '_' file_name_Go '.mat']), ...
            'soma_spkTimes')
        soma_spkTimesT = soma_spkTimes - (target_times_Go(1,t) ...
            + saccade_times_Go(1,t));
        soma_spkTimesT_tone = soma_spkTimesT(soma_spkTimesT < tone_times_Go(1,t));
        spikes_Go_saccade.timestamp{n+1} = cat(2, spikes_Go_saccade.timestamp{n+1}, ...
            (soma_spkTimesT_tone).*Fs+1);
        spikes_Go_saccade.time{n+1} = cat(2, spikes_Go_saccade.time{n+1}, ...
            (soma_spkTimesT_tone));
        spikes_Go_saccade.trial{n+1} = cat(2, spikes_Go_saccade.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimesT_tone)]));
    end
end

%% ---- NC trials ----
clearvars spikes_NC_saccade soma_spkTimesT soma_spkTimesT_tone
spikes_NC_saccade.label = cell(1, num_neurons); % name of the units
spikes_NC_saccade.timestamp = cell(1, num_neurons); % spike times in recording samples
spikes_NC_saccade.trial = cell(1, num_neurons);
spikes_NC_saccade.time = cell(1, num_neurons);
spikes_NC_saccade.trialtime = zeros(num_trials, 2);
spikes_NC_saccade.trialtime(:,1) = spikes_NC_saccade.trialtime(:,1) - (target_times_NC ...
    + saccade_times_NC)';
spikes_NC_saccade.trialtime(:,2) = simulation_length - (target_times_NC ...
    + saccade_times_NC)';
spikes_NC_saccade.timestampdimord = '{chan}_spike'; % dimension of the spiking data

for n=0:(num_neurons-1)
    
    spikes_NC_saccade.label{n+1} = ['neuron_' num2str(n)];
    
    for t=1:num_trials
        
        load(fullfile(path_simData_NC, ['NeuronsData_r' num2str(t) '_n#' ...
            '' num2str(n) '_' file_name_NC '.mat']), ...
            'soma_spkTimes')
        soma_spkTimesT = soma_spkTimes - (target_times_NC(1,t) ...
            + saccade_times_NC(1,t));
        soma_spkTimesT_tone = soma_spkTimesT(soma_spkTimesT < tone_times_NC(1,t));
        spikes_NC_saccade.timestamp{n+1} = cat(2, spikes_NC_saccade.timestamp{n+1}, ...
            (soma_spkTimesT_tone).*Fs+1);
        spikes_NC_saccade.time{n+1} = cat(2, spikes_NC_saccade.time{n+1}, ...
            (soma_spkTimesT_tone));
        spikes_NC_saccade.trial{n+1} = cat(2, spikes_NC_saccade.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimesT_tone)]));
    end
end

%% isi-pretarget
cfg = [];
cfg.bins = bins;
latensy_str = 'prestim'; % latency
cfg.latency = latensy_str;
isih_preTarget_Go = ft_spike_isi(cfg, spikes_Go_target);
isih_preTarge_NC = ft_spike_isi(cfg, spikes_NC_target);

%%
for ii=1:num_neurons
    figure('Visible', 'off');
    cfg.interpolate = 5;
    cfg.spikechannel = spikes_NC_saccade.label{ii};
    cfg.window = 'gausswin';
    cfg.winlen = 0.005;
    cfg.scatter = 'no';
    cfg.colormap = jet(300);
    
    subplot(1,2,1)
    [~, hdl_Go_preTarget.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_preTarget_Go);
    
    subplot(1,2,2)
    [~, hdl_NC_preTarget.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_preTarge_NC);
    
end

% close all

%% get the mean dens dist across neurons

dens_Go_preTarget = structfun(@(x) x.density.CData, hdl_Go_preTarget,...
    "UniformOutput",false);
dens_Go_preTarget = struct2cell(dens_Go_preTarget);
dens_Go_preTarget = cat(3,dens_Go_preTarget{:})./size(spikes_Go_target.trialtime, 1);

dens_NC_preTarget = structfun(@(x) x.density.CData, hdl_NC_preTarget,...
    "UniformOutput",false);
dens_NC_preTarget = struct2cell(dens_NC_preTarget);
dens_NC_preTarget = cat(3,dens_NC_preTarget{:})./size(spikes_Go_target.trialtime, 1);

isih_preTarget = mean(cat(3, isih_preTarget_Go.avg, ...
    isih_preTarge_NC.avg), 3, ...
    'omitnan')./size(spikes_Go_target.trialtime, 1);

isidens_preTarget = mean(mean(cat(4, dens_Go_preTarget, dens_NC_preTarget), ...
    4,'omitnan'), 3,'omitnan');

%% isi-postsaccade
cfg = [];
cfg.bins = bins;
latensy_str = 'poststim'; % latency
cfg.latency = latensy_str;
isih_postSaccade_Go = ft_spike_isi(cfg, spikes_Go_saccade);
isih_postSaccade_NC = ft_spike_isi(cfg, spikes_NC_saccade);

%%
for ii=1:num_neurons
    figure('Visible', 'off');
    cfg.interpolate = 5;
    cfg.spikechannel = spikes_NC_saccade.label{ii};
    cfg.window = 'gausswin';
    cfg.winlen = 0.005;
    cfg.scatter = 'no';
    cfg.colormap = jet(300);
    
    subplot(1,2,1)
    [~, hdl_Go_postSaccade.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_postSaccade_Go);
    
    subplot(1,2,2)
    [~, hdl_NC_postSaccade.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_postSaccade_NC);
    
end

%% get mean isi dens across neurons

dens_Go_postSaccade = structfun(@(x) x.density.CData, hdl_Go_postSaccade,...
    "UniformOutput",false);
dens_Go_postSaccade = struct2cell(dens_Go_postSaccade);
dens_Go_postSaccade = cat(3,dens_Go_postSaccade{:})./size(spikes_Go_target.trialtime, 1);

dens_NC_postSaccade = structfun(@(x) x.density.CData, hdl_NC_postSaccade,...
    "UniformOutput",false);
dens_NC_postSaccade = struct2cell(dens_NC_postSaccade);
dens_NC_postSaccade = cat(3,dens_NC_postSaccade{:})./size(spikes_Go_target.trialtime, 1);

isidens_Go_postSaccade = mean(dens_Go_postSaccade, 3,'omitnan');
isidens_NC_postSaccade = mean(dens_NC_postSaccade, 3,'omitnan');

close all

%% save ISI dists

save(fullfile('data', 'sim_data','L3_isi_dist_sim.mat'), ...
    'isih_preTarget', 'isidens_preTarget', ...
    'isih_postSaccade_Go', 'isih_postSaccade_NC', ...
    'bins_middle', 'isidens_Go_postSaccade', 'isidens_NC_postSaccade')

save(fullfile('data', 'sim_data','L3_spk_tone_ft_struct_sim.mat'), ...
    'spikes_Go_target', 'spikes_NC_target', ...
    'spikes_Go_saccade', 'spikes_NC_saccade')
