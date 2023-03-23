%% process sim data for Fig 5

clear
clc

POPULATION_SIZE = 100;

% specify pop of neurons to be analyzed
population_type_stim = 'L3PCs_off'; %
% options: 
%   L3PCs_off or L5PCs_off -> no synchronized stimulus in either of the
%       populations.
%   L3PCs_on or L5PCs_on -> synchronized stimulus at 1s.
runNumb = 1;

%% load simulated data: extracellular potentials, membrane potentials, etc...

switch population_type_stim
    case 'L3PCs_off'
        data_folder = ['D:\theta-paper\results_L3PCsPopMky\' ...
            'Go_trial_Feb23_2023_rhythmicity\' ...
            'neurons#100_distributed_synp\' ...
            'StimDend#4_StimApic#4'];
        fileName = 'Dend_r2_Apic_r1';
    case 'L5PCs_off'
        data_folder = ['D:\Theta_paper_sim\results_L5PCsPopMky\' ...
            'Go_trial_Feb23_2023_rhythmicity\' ...
            'neurons#100_distributed_synp\StimDend#4_StimOblq#4_StimApic#4'];
        fileName = 'Dend_r5_Oblq_4_Apic_r1';
    case 'L3PCs_on'
        data_folder = ['D:\theta-paper\results_L3PCsPopMky\' ...
            'Go_trial_Feb23_2023_rhythmicity_stim\' ...
            'neurons#100_distributed_synp\' ...
            'StimDend#4_StimApic#4'];
        fileName = 'Dend_r2_Apic_r1';
    case 'L5PCs_on'
        data_folder = ['D:\Theta_paper_sim\results_L5PCsPopMky\' ...
            'Go_trial_Feb23_2023_rhythmicity\' ...
            'neurons#100_distributed_synp\StimDend#4_StimOblq#4_StimApic#4'];
        fileName = 'Dend_r5_Oblq_4_Apic_r1';
    otherwise
        error('')
end

% - path and file name of the simulated data
loadPath2 = fullfile(data_folder, ['SimData_r' num2str(runNumb) '_PS' ...
    num2str(POPULATION_SIZE) ...
    '_' fileName '.mat']);

load(loadPath2, 'LFP', 'ze') % load extracellular potentials

%% simulation details

warmup_period = 500e-3; % sec
simulation_length = 10000e-3; % sec, simulation length
dt = (2^(-4))*1e-3;  %[s] time increment
Fs = 1/(dt); % [Hz] sampling frequency
ts_all = 0:dt:(warmup_period + simulation_length);  %[s] time span,
% whole simulation
ts = 0:dt:simulation_length;  %[s] time span neglecting the warmup period

Fs_lfp = 1e3; % targeted LFP sampling rate

%% pre-allocating memory for APs and Ca spks times cell-array
APs_times = cell(POPULATION_SIZE,1);
Caspks_times = cell(POPULATION_SIZE,1);

%% loading data and performing pre-processing

for ii = 1:POPULATION_SIZE % loop for the cells

    neuron_num = ii - 1;

    sprintf('Neuron #%d', ii)

    % load somatic and dendritic membrane potentials
    loadPath1 = fullfile(data_folder,...
        ['NeuronsData_r' num2str(runNumb) '_n#'...
        num2str(neuron_num) '_'...
        fileName '.mat']);

    load(loadPath1, 'Vs', 'v_mbp', 'soma_spkTimes', 'dend_spkTimes')

    %% saving spikes times

    APs_times{ii} = soma_spkTimes;
    Caspks_times{ii} = dend_spkTimes;

    %% remove warm up period
    Vs = Vs(ts_all>=warmup_period);
    v_mbp = v_mbp(ts_all>=warmup_period);

    %% amplitude/phase analysis on the membrane potential

    % -- filter at theta and alpha bands
    Vs_theta = eegfilt(Vs, Fs, 4, 8);

    v_mbp_theta = eegfilt(v_mbp, Fs, 4, 8);

    % -- compute Hilbert Transform
    Vs_theta_hilbert = hilbert(Vs_theta);
    Vmbp_theta_hilbert = hilbert(v_mbp_theta);

    % -- compute instantaneous phase and amplitude envelope
    Phi_VsTheta = angle(Vs_theta_hilbert); % time series of the phases of
    % the signal
    Amp_VsTheta = abs(Vs_theta_hilbert); % time series of the amplitude
    % envelope of the signal

    Phi_VmbpTheta = angle(Vmbp_theta_hilbert); % time series of the phases
    % of the signal
    Amp_VmbpTheta = abs(Vmbp_theta_hilbert); % time series of the amplitude
    % envelope of the signal

    %% store data

    if ii==1
        Phi_VsThetaR = Phi_VsTheta;
        Amp_VsThetaR = Amp_VsTheta;
        Phi_VdThetaR = Phi_VmbpTheta;
        Amp_VdThetaR = Amp_VmbpTheta;

        VsR = Vs;
        VdR = v_mbp;
    else
        Phi_VsThetaR = cat(1, Phi_VsThetaR, Phi_VsTheta);
        Amp_VsThetaR = cat(1, Amp_VsThetaR, Amp_VsTheta);
        Phi_VdThetaR = cat(1, Phi_VdThetaR, Phi_VmbpTheta);
        Amp_VdThetaR = cat(1, Amp_VdThetaR, Amp_VmbpTheta);

        VsR = cat(1, VsR, Vs);
        VdR = cat(1, VdR, v_mbp);
    end

end

%% save Phase/Amp Data and membrane potentials of all cells

saveName = fullfile(data_folder, ['PhaseAmpData_MemPot_' fileName '.mat']);
save(saveName, 'Phi_VsThetaR', 'Amp_VsThetaR', 'Phi_VdThetaR', 'Amp_VdThetaR')

saveName = fullfile(data_folder, ['MemPotentials_' fileName '.mat']);
save(saveName, 'VsR', 'VdR')

%% save the spike times -> APs and Ca spikes

saveName = fullfile(data_folder, ['SpikeTimes_' fileName '.mat']);
save(saveName, 'APs_times', 'Caspks_times')

%% process the extracellular potentials

% -- downsampling the LFPs to 1kHz
Fs_new = 1e3; % 1kHz new sampling frequency
[VeD, ts_out] = process_resample('Compute', double(LFP(:, ts_all>=warmup_period)), ...
    ts, Fs_new);

%% create FieldTrip str for the LFP data
% relative to the start of the trial
% data is low-pass filtered at 100 Hz using a two-pass fourth order Butterworth filter
% to prevent spike artifacts in the LFP (Voloh and Womelsdorf 2018 Cerebral Cortex)

ft_lfp = convert2ftLFPDataStr(VeD, Fs_new, Fs, 1);

%% save the LFP

saveName = fullfile(data_folder, ['LFP1kHzAndMUA_' fileName '.mat']);
save(saveName, 'VeD','ts_out','Fs_new', 'ft_lfp')

%% get power | membrane potential

[p_Vs, f_Vs] = pspectrum(VsR',Fs,'FrequencyLimits',[2 30],...
    'FrequencyResolution', 1);
[p_Vd, f_Vd] = pspectrum(VdR',Fs,'FrequencyLimits',[2 30],...
    'FrequencyResolution', 1);

%% Fig 5 - membrane potential plots

font = 10;

figure('Units', 'inches','Position',[0 0 4.5 5.5]);
tiledlayout(3, 1,'TileSpacing','Compact','Padding','Compact');

nexttile;
plot(ts, VsR(1,:), '-', 'LineWidth',1.5, 'Color',[0 128 255]./255)
[~,locs] = findpeaks(VsR(1,:),ts,"MinPeakProminence",60);
hold on;
plot(locs, 50.*ones(length(locs)),'.k')
hold on;
plot(ts, VdR(1,:), '-', 'LineWidth',1.5, 'Color',[255 102 102]./255)
box off
xlim([0 2])
ylim([-100 60])
xlabel('Time (s)')
ylabel('mV')
title('Membrane Potential')
% legend({'Soma', 'Dendrites', 'Location', 'best'})
set(gca,'linewidth',1,'fontsize',font,'fontweight','bold')

nexttile;
hold on;
plot(f_Vs,(p_Vs), ...
    'LineWidth', 1.5)
hold off;
xticks([1:9, 10:5:60])
xticklabels({})
xlim([2 30])
box off
title('Somatic Membrane Potential')
ylabel({'Power (mV^2)'})
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

nexttile;
hold on;
plot(f_Vd,(p_Vd), ...
    'LineWidth', 1.5)
hold off;
xticks([1:9, 10:5:60])
xticksLbs = {''};
for jj=[2:9, 10:5:60]
    if jj ==5 || jj>=10
        xticksLbs = cat(1,xticksLbs,{sprintf('%.0f',jj)});
    else
        xticksLbs = cat(1,xticksLbs,{''});
    end
end
xticklabels(xticksLbs)
xlim([2 30])
title('Dendritic Membrane Potential')
xlabel('Frequency (Hz)')
box off
ylabel({'Power (mV^2)'})
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% Fig 5 - plot theta phase

font = 10;

figure('Units', 'inches','Position',[0 0 4.5 4]);
tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');

nexttile;
hold on;
plot(ts, Phi_VsThetaR', '-', 'LineWidth',1.5)
box off
xlim([0 2])
xticklabels({})
% xlabel('Time (s)')
ylabel('Soma')
title('Membrane Potential Theta Phase (rad)')
% legend({'Soma', 'Dendrites', 'Location', 'best'})
set(gca,'linewidth',1,'fontsize',font,'fontweight','bold')

nexttile;
hold on;
plot(ts, Phi_VdThetaR', '-', 'LineWidth',1.5)
box off
xlim([0 2])
xlabel('Time (s)')
ylabel('Dendrites')
% legend({'Soma', 'Dendrites', 'Location', 'best'})
set(gca,'linewidth',1,'fontsize',font,'fontweight','bold')

%% LFP power spectrum

[p_LFP, f_LFP] = pspectrum(ft_lfp.trial{1,1}'.*1e3,Fs_lfp,'FrequencyLimits',[2 30],...
    'FrequencyResolution', 1);

%% plot

font = 10;

figure('Units', 'inches','Position',[0.05 0.05 4 1.5]);
plot(f_LFP,(p_LFP), ...
    'LineWidth', 1.25)
xticks([1:9, 10:5:60])
xticksLbs = {''};
for jj=[2:9, 10:5:60]
    if jj ==5 || jj>=10
        xticksLbs = cat(1,xticksLbs,{sprintf('%.0f',jj)});
    else
        xticksLbs = cat(1,xticksLbs,{''});
    end
end
xticklabels(xticksLbs)
xlim([4 30])
box off
title('LFP')
xlabel('Frequency (Hz)')
ylabel({'Power (\muV^2)'})
set(gca,'linewidth',1.5,'fontsize',9,'fontweight','bold')

%% calculate the lfp laminar theta power

% baseline correction
LFP_theta = eegfilt(ft_lfp.trial{1,1} - mean(ft_lfp.trial{1,1}(:,1:200),2), ...
    Fs_lfp, 4, 8);

% get hilbert transform
LFP_theta_hilbert = hilbert(LFP_theta')';
% calculate the power
power_VeDtheta = abs(LFP_theta_hilbert).^2;

%% plot

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
imagesc(ts_out, 1:16, power_VeDtheta);
xlim([0 2])
hold off
c = colorbar;
c.Label.String = '\muV^2';
ylabel('Electrodes')
title('LFP theta power')
colormap(jet);
xlabel('Time (s)')
set(gca,'linewidth',2,'fontsize',10,'fontweight','bold')
