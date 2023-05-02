%% calculate the spk rate for L3 sim neurons

clear
clc

%% path to data

path_simData_Go = ['D:\Theta_paper_sim\results_L3PCsPopMky\' ...
    'Go_trial_Sept9_mpi\neurons#10_clustered_synp\' ...
    'StimDend#4_StimApic#4'];

path_simData_NC = ['D:\Theta_paper_sim\results_L3PCsPopMky\' ...
    'NC_trial_Sept9_mpi\neurons#10_clustered_synp\' ...
    'StimDend#4_StimApic#4'];

file_name = 'Dend_r3.5_Apic_r2';

%% load events time

load(fullfile('data', 'sim_data', 'events_timing_Go.mat'))
saccade_times_Go = saccade_times;
target_times_Go = target_times;

load(fullfile('data', 'sim_data', 'events_timing_NC.mat'))
saccade_times_NC = saccade_times;
target_times_NC = target_times;

%% simulation details
target_times_Go = double(target_times_Go(1,:)).*1e-3; % sec, target times
saccade_times_Go = double(saccade_times_Go(1,:)).*1e-3; % sec, 
% saccade times relative to target
target_times_NC = double(target_times_NC(1,:)).*1e-3; % sec, target times
saccade_times_NC = double(saccade_times_NC(1,:)).*1e-3; % sec, 
% saccade times relative to target
num_trials = 116; % number of simulated trials
num_neurons = 10; % number of simulated neurons
warmup_period = 500e-3; % sec
simulation_length = 10000e-3; % simulation length
dt = (2^(-4))*1e-3;  %[s] time increment
Fs = 1/(dt); % [Hz] sampling frequency

%% parameters for firing rate calculation
binsize =  10e-3; % [s] bin size
gauss_SD = 10; %/binsize; % 0.02 seconds (20ms) SD
flt = true; % if True: use Gaussian filter in smoothdata;
% otherwise: convol with Gaussian

%% ========== creat field-trip spike structure ==========
%% define trials relative to saccade
%% ---- Go trials ----
clearvars spikes_Go_saccade
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
            '' num2str(n) '_' file_name '.mat']), ...
            'soma_spkTimes')
        spikes_Go_saccade.timestamp{n+1} = cat(2, spikes_Go_saccade.timestamp{n+1}, ...
            (soma_spkTimes - (target_times_Go(1,t) + saccade_times_Go(1,t))).*Fs+1);
        spikes_Go_saccade.time{n+1} = cat(2, spikes_Go_saccade.time{n+1}, ...
            (soma_spkTimes - (target_times_Go(1,t) + saccade_times_Go(1,t))));
        spikes_Go_saccade.trial{n+1} = cat(2, spikes_Go_saccade.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimes)]));
    end
end

%% ---- NC trials ----
clearvars spikes_NC_saccade
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
            '' num2str(n) '_' file_name '.mat']), ...
            'soma_spkTimes')
        spikes_NC_saccade.timestamp{n+1} = cat(2, spikes_NC_saccade.timestamp{n+1}, ...
            (soma_spkTimes - (target_times_NC(1,t) + saccade_times_NC(1,t))).*Fs+1);
        spikes_NC_saccade.time{n+1} = cat(2, spikes_NC_saccade.time{n+1}, ...
            (soma_spkTimes - (target_times_NC(1,t) + saccade_times_NC(1,t))));
        spikes_NC_saccade.trial{n+1} = cat(2, spikes_NC_saccade.trial{n+1}, ...
            repmat(t, [1, length(soma_spkTimes)]));
    end
end

%% calculate firing rate relative to saccade

clearvars rate_Go_saccade rate_NC_saccade
ts_edges_saccade = (-0.5):binsize:2;
ts_saccade = ts_edges_saccade + mean(diff(ts_edges_saccade), 'omitnan');
ts_saccade = ts_saccade(1:end-1);
rate_Go_saccade = cal_spkrate_ftspkstr(spikes_Go_saccade, ts_edges_saccade, ...
    gauss_SD + 1, flt);
rate_NC_saccade = cal_spkrate_ftspkstr(spikes_NC_saccade, ts_edges_saccade, ...
    gauss_SD + 1, flt);

%% save data

save(fullfile('data', 'sim_data', 'L3_spk_rate_sim.mat'), 'rate_NC_saccade', ...
    'rate_Go_saccade', 'ts_saccade')
save(fullfile('data', 'sim_data', 'L3_spk_ft_struct_sim.mat'), 'spikes_NC_saccade', ...
    'spikes_Go_saccade')
