%% process sim data for Fig 6

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
        data_folder = fullfile('data', 'sim_data', 'L3PCs_off');
        fileName = 'Dend_r2_Apic_r1';
    case 'L5PCs_off'
        data_folder = fullfile('data', 'sim_data', 'L5PCs_off');
        fileName = 'Dend_r5_Oblq_4_Apic_r1';
    case 'L3PCs_on'
        data_folder = fullfile('data', 'sim_data', 'L3PCs_on');
        fileName = 'Dend_r2_Apic_r1';
    case 'L5PCs_on'
        data_folder = fullfile('data', 'sim_data', 'L5PCs_on');
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

%% loading data and performing pre-processing

for ii = 1:POPULATION_SIZE % loop for the cells
    
    neuron_num = ii - 1;
    
    sprintf('Neuron #%d', ii)
    
    % load somatic and dendritic membrane potentials
    loadPath1 = fullfile(data_folder,...
        ['NeuronsData_r' num2str(runNumb) '_n#'...
        num2str(neuron_num) '_'...
        fileName '.mat']);
    
    load(loadPath1, 'Vs', 'v_mbp')
    
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

%% fit a/f^b curve to the power spectrum of membrane potentials
% if population_type_stim = 'L3PCs_off' or 'L3PCs_on'

if strcmp(population_type_stim, 'L3PCs_off') || strcmp(population_type_stim, 'L3PCs_on')
    % somatic membrane potential
    p_Vs_mean = mean(p_Vs, 2, "omitnan")';
    p_Vs_std = std(p_Vs, 0, 2, "omitnan")';
    
    % dendritic membrane potential
    p_Vd_mean = mean(p_Vd, 2, "omitnan")';
    p_Vd_std = std(p_Vd, 0, 2, "omitnan")';
    
    % fit type
    ft=fittype(@(a, b, x) a./(x.^b), 'coefficients',{'a', 'b'},'independent',{'x'});
    
    % fitting
    [fitobject_pVs, gof_pVs] = fit(f_Vs, ...
        p_Vs_mean', ft);
    [fitobject_pVd, gof_pVd] = fit(f_Vd, ...
        p_Vd_mean', ft);
    
end

%% Fig 5 & suppl fig - membrane potential plots

font = 10;

figure('Units', 'inches','Position',[0 0 4.5 5.5]);
tiledlayout(2, 1,'TileSpacing','Compact','Padding','Compact');

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
% soma
x2 = [f_Vs', fliplr(f_Vs')];
inBetween1 = [p_Vs_mean,...
    fliplr(p_Vs_mean-(1.96.*p_Vs_std))];
fill(x2, inBetween1, [0 128 255]./255, 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
inBetween2 = [p_Vs_mean+(1.96.*p_Vs_std),...
    fliplr(p_Vs_mean)];
fill(x2, inBetween2, [0 128 255]./255, 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
p1 = plot(f_Vs',p_Vs_mean, ...
    'LineWidth', 1.5, 'Color', [0 128 255]./255);
if strcmp(population_type_stim, 'L3PCs_off') || strcmp(population_type_stim, 'L3PCs_on')
    p2 = plot(f_Vs',feval(fitobject_pVs, f_Vs'), '--', ...
        'LineWidth', 2, 'Color', [0 204 0]./255);
end
% dendrites
hold on;
x2 = [f_Vs', fliplr(f_Vs')];
inBetween1 = [p_Vd_mean,...
    fliplr(p_Vd_mean-(1.96.*p_Vd_std))];
fill(x2, inBetween1, [255 102 102]./255, 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
inBetween2 = [p_Vd_mean+(1.96.*p_Vd_std),...
    fliplr(p_Vd_mean)];
fill(x2, inBetween2, [255 102 102]./255, 'FaceAlpha', 0.2, ...
    'EdgeColor', 'none');
p3 = plot(f_Vd,p_Vd_mean, ...
    'LineWidth', 1.5, 'Color', [255 102 102]./255);
if strcmp(population_type_stim, 'L3PCs_off') || strcmp(population_type_stim, 'L3PCs_on')
    p4 = plot(f_Vd',feval(fitobject_pVd, f_Vd'), '--', ...
        'LineWidth', 2, 'Color', [0 204 0]./255);
end
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
if strcmp(population_type_stim, 'L3PCs_off') || strcmp(population_type_stim, 'L3PCs_on')
    legend([p1 p2 p3 p4], {'Soma', sprintf(['a/f^b fit: a = %.2f, b = %.2f ' ...
        '| R^2 = %.2f'], round(fitobject_pVs.a,2), ...
        round(fitobject_pVs.b,2),round(gof_pVs.rsquare,2)), ...
        'Dendrites', sprintf(['a/f^b fit: a = %.2f, b = %.2f ' ...
        '| R^2 = %.2f'], round(fitobject_pVd.a,2), ...
        round(fitobject_pVd.b,2),round(gof_pVd.rsquare,2))}, 'box', 'off')
else
    legend([p1 p3], {'Soma', 'Dendrites'}, 'box', 'off')
end
xticklabels(xticksLbs)
xlim([2 20])
xlabel('Frequency (Hz)')
box off
ylabel({'Power (mV^2)'})
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% calculate the LFP laminar theta power

% baseline correction and filtering of the LFP
LFP_theta = eegfilt(ft_lfp.trial{1,1} - mean(ft_lfp.trial{1,1}(:,1:200),2), ...
    Fs_lfp, 4, 8);

% get hilbert transform
LFP_theta_hilbert = hilbert(LFP_theta')';
% calculate the power
power_VeDtheta = abs(LFP_theta_hilbert).^2;

%% Fig 5 - plot theta phase and LFP laminar theta power

font = 10;

figure('Units', 'inches','Position',[0 0 4.5 5]);
tiledlayout(3, 1,'TileSpacing','Compact','Padding','Compact');

t1 = nexttile;
imagesc(ts(ts>=0 & ts<=2) - 1, 1:100, Phi_VsThetaR(:,(ts>=0 & ts<=2)));
xticklabels({})
colormap(t1, "hsv")
caxis([-3.1416 3.1416])
c = colorbar;
c.Label.String = 'Theta Phase';
c.LineWidth = 1.5;
c.Ticks = [-3.1416 0 3.1416];
c.TickLabels = {'-\pi', '0', '\pi'};
c.FontSize = font;
c.FontWeight = 'bold';
ylabel('Neurons')
title('Somatic Membrane Potential')
set(gca,'linewidth',1.75,'fontsize',font,'fontweight','bold')

t2 = nexttile;
imagesc(ts(ts>=0 & ts<=2) - 1, 1:100, Phi_VdThetaR(:,(ts>=0 & ts<=2)));
xticklabels({})
ylabel('Neurons')
colormap(t2, "hsv")
caxis([-3.1416 3.1416])
c = colorbar;
c.Label.String = 'Theta Phase';
c.LineWidth = 1.5;
c.Ticks = [-3.1416 0 3.1416];
c.TickLabels = {'-\pi', '0', '\pi'};
c.FontSize = font;
c.FontWeight = 'bold';
ylabel('Neurons')
title('Dendritic Membrane Potential')
set(gca,'linewidth',1.75,'fontsize',font,'fontweight','bold')

t3 = nexttile;
imagesc(ts_out(ts_out>=0 & ts_out<=2) - 1, 1:16, power_interp');
c = colorbar;
c.Label.String = '\muV^2';
c.FontSize = font;
c.FontWeight = 'bold';
c.LineWidth = 1.5;
ylabel('Electrodes')
title('LFP theta power')
colormap(t3, jet);
xlabel('Time from Stimulus Onset (s)')
set(gca,'linewidth',1.75,'fontsize',font,'fontweight','bold')
