%% make Fig 2D - L3 spk rate

clear
clc

rng('default')

%% load spk rate - experimental data

load(fullfile("data/monkeys_data/idx_L3neurons_sac.mat"), 'idx_4'); 

path2data = fullfile('data', 'monkeys_data','spkrate_errorneurons');
load(fullfile(path2data, 'fr_error_cells_saccade.mat'), 'mean_fr_L3_NC', ...
    'mean_fr_L3_Go', 'ts_sess')
cells_L3_mean_fr_NC_saccade = mean_fr_L3_NC(idx_4, :);
cells_L3_mean_fr_Go_saccade = mean_fr_L3_Go(idx_4, :);
L3_mean_fr_NC_saccade = mean(cells_L3_mean_fr_NC_saccade, 1, "omitnan");
L3_mean_fr_Go_saccade = mean(cells_L3_mean_fr_Go_saccade, 1, "omitnan");
L3_sem_fr_NC_saccade = std(cells_L3_mean_fr_NC_saccade, 0, 1, ...
    "omitnan")./sqrt(length(idx_4));
L3_sem_fr_Go_saccade = std(cells_L3_mean_fr_Go_saccade, 0, 1, ...
    "omitnan")./sqrt(length(idx_4));
ts_saccade_exp = ts_sess{1};

%% load spk rate - sim data

load(fullfile('data', 'sim_data','L3_spk_rate_sim.mat'), 'rate_NC_saccade', ...
    'rate_Go_saccade', 'ts_saccade')

mean_fr_Go_saccade = mean(cell2mat(rate_Go_saccade.firingRatet'), 1, ...
    "omitnan");
mean_fr_NC_saccade = mean(cell2mat(rate_NC_saccade.firingRatet'), 1, ...
    "omitnan");
sem_fr_Go_saccade = std(cell2mat(rate_Go_saccade.firingRatet'), 0, 1, ...
    "omitnan")./sqrt(length(rate_Go_saccade.firingRatet));
sem_fr_NC_saccade = std(cell2mat(rate_NC_saccade.firingRatet'), 0, 1, ...
    "omitnan")./sqrt(length(rate_NC_saccade.firingRatet));

ref_saccade = mean([mean(mean_fr_Go_saccade(ts_saccade<-0.2), "omitnan") ...
    mean(mean_fr_NC_saccade(ts_saccade<-0.2), "omitnan")]);

%% get mean spk rate features

[pks_Go_saccade_exp,locs_Go_saccade_exp,widths_Go_saccade_exp,proms_Go_saccade_exp] = findpeaks( ...
    L3_mean_fr_Go_saccade - ...
    L3_mean_fr_Go_saccade(1), ...
    ts_saccade_exp,'MinPeakProminence', 10, 'WidthReference', 'halfheight');
[pks_NC_saccade_exp,locs_NC_saccade_exp,widths_NC_saccade_exp,proms_NC_saccade_exp] = findpeaks( ...
    L3_mean_fr_NC_saccade - ...
    L3_mean_fr_NC_saccade(2), ...
    ts_saccade_exp,'MinPeakProminence', 10, 'WidthReference', 'halfheight');

[pks_Go_saccade_sim,locs_Go_saccade_sim,widths_Go_saccade_sim,proms_Go_saccade_sim] = findpeaks( ...
    mean_fr_Go_saccade - ref_saccade, ...
    ts_saccade,'MinPeakProminence', 10, 'WidthReference', 'halfheight');
[pks_NC_saccade_sim,locs_NC_saccade_sim,widths_NC_saccade_sim,proms_NC_saccade_sim] = findpeaks( ...
    mean_fr_NC_saccade - ref_saccade, ...
    ts_saccade,'MinPeakProminence', 10, 'WidthReference', 'halfheight');

%% plot 

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
hold on;
plot(ts_saccade.*1e3, mean_fr_NC_saccade - ref_saccade, ':r', ...
    'linewidth',3);
plot(ts_saccade.*1e3, mean_fr_Go_saccade - ref_saccade, '-r', ...
    'linewidth', 1.25);

plot(ts_saccade.*1e3, L3_mean_fr_NC_saccade - L3_mean_fr_Go_saccade(2), ':k', ...
    'linewidth',3);
plot(ts_saccade.*1e3, L3_mean_fr_Go_saccade - L3_mean_fr_Go_saccade(1), '-k', ...
    'linewidth', 1.25);
xlabel('Time from saccade (ms)')
ylabel('Relative firing rate (spks/s)')
xlim([-0.5 0.5].*1e3)
ylim([-1 45])
xticks([-400 -200 0 200 400])
legend({'Error - Sim', 'Correct - Sim', 'Error - Exp', 'Correct - Exp'}, ...
    'Location','best', 'Box','off')
set(gca, 'box', 'off','linewidth',2,'fontsize',10,'fontweight','bold')

%% calculate peak features across cells
%% exp data

idx_ts_0_exp = find(ts_saccade>=0, 1, 'first');
idx_ts_500_exp = find(ts_saccade>=0.5, 1, 'first');

cells_ref_L3_saccade_exp = min([cells_L3_mean_fr_Go_saccade(:,1) ...
    cells_L3_mean_fr_NC_saccade(:,1)], [], 2);
[cells_peak_Go_saccade_exp, cells_peak_idx_Go_saccade_exp] = max( ...
    cells_L3_mean_fr_Go_saccade(:, idx_ts_0_exp:idx_ts_500_exp) ...
    - cells_ref_L3_saccade_exp, [], 2);
[cells_peak_NC_saccade_exp, cells_peak_idx_NC_saccade_exp] = max( ...
    cells_L3_mean_fr_NC_saccade(:, idx_ts_0_exp:idx_ts_500_exp) ...
    - cells_ref_L3_saccade_exp, [], 2);

% peak location and half width
[~,locs_cells_Go_saccade_exp,widths_cells_Go_saccade_exp_all,~] = cellfun(@(x) ...
    findpeaks(x, ts_saccade(1:idx_ts_500_exp), ...
    'MinPeakProminence', 3, 'WidthReference', 'halfheight'), ...
    num2cell(cells_L3_mean_fr_Go_saccade(:, 1:idx_ts_500_exp) ...
    - cells_ref_L3_saccade_exp, 2), 'UniformOutput',false);

locs_idx = cellfun(@(x) find(x>=0, 1, 'first'), locs_cells_Go_saccade_exp);
widths_cells_Go_saccade_exp = cellfun(@(x, idx) x(idx), ...
    widths_cells_Go_saccade_exp_all, num2cell(locs_idx,2));

[~,~,widths_cells_NC_saccade_exp,~] = cellfun(@(x) ...
    findpeaks(x, ts_saccade(1:idx_ts_500_exp), ...
    'MinPeakProminence', 10, 'WidthReference', 'halfheight'), ...
    num2cell(cells_L3_mean_fr_NC_saccade(:, 1:idx_ts_500_exp) ...
    - cells_ref_L3_saccade_exp, 2));

%% sim data
idx_ts_500 = find(ts_saccade>=0.5, 1, 'first');
cells_L3_mean_fr_Go_saccade_sim = cell2mat(rate_Go_saccade.firingRatet');
cells_L3_mean_fr_NC_saccade_sim = cell2mat(rate_NC_saccade.firingRatet');

% get peak | global max
[cells_peak_Go_saccade_sim, cells_peak_idx_Go_saccade_sim] = max( ...
    cells_L3_mean_fr_Go_saccade_sim ...
    - ref_saccade, [], 2); 
[cells_peak_NC_saccade_sim, cells_peak_idx_NC_saccade_sim] = max( ...
    cells_L3_mean_fr_NC_saccade_sim ...
    - ref_saccade, [], 2);

% get latency and half width of the peak
[~,locs_cells_Go_saccade_sim,widths_cells_Go_saccade_sim_all,~] = cellfun(@(x) ...
    findpeaks(x, ts_saccade(1:idx_ts_500_exp), ...
    'MinPeakProminence', 3, 'WidthReference', 'halfheight'), ...
    num2cell(cells_L3_mean_fr_Go_saccade_sim(:, 1:idx_ts_500_exp) ...
    - ref_saccade, 2), 'UniformOutput',false);

locs_idx_sim = cellfun(@(x) find(x>=0, 1, 'first'), locs_cells_Go_saccade_sim);
widths_cells_Go_saccade_sim = cellfun(@(x, idx) x(idx), ...
    widths_cells_Go_saccade_sim_all, num2cell(locs_idx_sim,2));

[~,~,widths_cells_NC_saccade_sim,~] = cellfun(@(x) ...
    findpeaks(x, ts_saccade(1:idx_ts_500_exp), ...
    'MinPeakProminence', 10, 'WidthReference', 'halfheight'), ...
    num2cell(cells_L3_mean_fr_NC_saccade_sim(:, 1:idx_ts_500_exp) ...
    - ref_saccade, 2));

%% get dist of values across cells 

font = 8;

figure('Units', 'inches','Position',[0.05 0.05 3 4]);

% amplitude
subplot(3,1,1)
hold on;
b1 = barh([mean(cells_peak_Go_saccade_exp, 'omitnan') ...
    mean(cells_peak_NC_saccade_exp, 'omitnan'); 0 0]);
b1(1).FaceColor = [1 1 1].*224/255;
b1(2).FaceColor = [1 1 1].*96/255;
b1(1).LineWidth = 1.25;
b1(2).LineWidth = 3;
b1(1).LineStyle = '-';
b1(2).LineStyle = ':';
ylim([0.5 1.5])

errorbar(mean(cells_peak_Go_saccade_exp, ...
    'omitnan'), 1 - 0.15, std(cells_peak_Go_saccade_exp, ...
    "omitnan"), 'horizontal', ...
    '.k','linewidth', 1.25)
errorbar(mean(cells_peak_NC_saccade_exp, 'omitnan'), 1 + 0.15, ...
    std(cells_peak_NC_saccade_exp, ...
    "omitnan"), 'horizontal', '.k', 'linewidth', 3)

plot(cells_peak_Go_saccade_sim, 1 - 0.15, ...
    '.r','MarkerSize', 10)
hold on;
plot(cells_peak_NC_saccade_sim, 1 + 0.15, ...
    '.r','MarkerSize', 10)

yticklabels({})
% title('Amplitude')
xlabel('Amplitude (spks/s)')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

% latency
subplot(3,1,2)
hold on;
b1 = barh([mean(ts_saccade(cells_peak_idx_Go_saccade_exp+idx_ts_0_exp-1), 'omitnan').*1e3 ...
    mean(ts_saccade(cells_peak_idx_NC_saccade_exp+idx_ts_0_exp-1), 'omitnan').*1e3; 0 0]);
b1(1).FaceColor = [1 1 1].*224/255;
b1(2).FaceColor = [1 1 1].*96/255;
b1(1).LineWidth = 1;
b1(2).LineWidth = 3;
b1(1).LineStyle = '-';
b1(2).LineStyle = ':';
ylim([0.5 1.5])

errorbar(mean(ts_saccade(cells_peak_idx_Go_saccade_exp+idx_ts_0_exp-1), ...
    'omitnan').*1e3, 1 - 0.15, std(ts_saccade(cells_peak_idx_Go_saccade_exp+idx_ts_0_exp-1).*1e3, ...
    "omitnan"), 'horizontal', ...
    '.k','linewidth', 1)
errorbar(mean(ts_saccade(cells_peak_idx_NC_saccade_exp+idx_ts_0_exp-1), ...
    'omitnan').*1e3, 1 + 0.15, ...
    std(ts_saccade(cells_peak_idx_NC_saccade_exp+idx_ts_0_exp-1).*1e3, ...
    "omitnan"), 'horizontal', '.k', 'linewidth', 3)

plot(ts_saccade(cells_peak_idx_Go_saccade_sim)'.*1e3, ...
    1 - 0.15, ...
    '.r','MarkerSize', 10)
hold on;
plot(ts_saccade(cells_peak_idx_NC_saccade_sim)'.*1e3, ...
    1 + 0.15, ...
    '.r','MarkerSize', 10)

yticklabels({})
% title('Peak Latency')
xlabel('Peak Latency (ms)')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

% half-width
subplot(3,1,3)
hold on;
b1 = barh([mean(widths_cells_Go_saccade_exp, 'omitnan').*1e3 ...
    mean(widths_cells_NC_saccade_exp, 'omitnan').*1e3; 0 0]);
b1(1).FaceColor = [1 1 1].*224/255;
b1(2).FaceColor = [1 1 1].*96/255;
b1(1).LineWidth = 1;
b1(2).LineWidth = 3;
b1(1).LineStyle = '-';
b1(2).LineStyle = ':';
ylim([0.5 1.5])

errorbar(mean(widths_cells_Go_saccade_exp, ...
    'omitnan').*1e3, 1 - 0.15, std(widths_cells_Go_saccade_exp.*1e3, ...
    "omitnan"), 'horizontal', ...
    '.k','linewidth', 1)
errorbar(mean(widths_cells_NC_saccade_exp, ...
    'omitnan').*1e3, 1 + 0.15, ...
    std(widths_cells_NC_saccade_exp.*1e3, ...
    "omitnan"), 'horizontal', '.k', 'linewidth', 3)

plot(widths_cells_Go_saccade_sim.*1e3, ...
    1 - 0.15, ...
    '.r','MarkerSize', 10)
hold on;
plot(widths_cells_NC_saccade_sim.*1e3, ...
    1 + 0.15, ...
    '.r','MarkerSize', 10)

yticklabels({})
% title('Peak Latency')
xlabel('Peak Half-width (ms)')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

%% MC test p-values | NC greater than Go | 100000 permutations
p_peak_amp_Go_simvsdata = mc_perm_test(cells_peak_Go_saccade_sim, ...
    cells_peak_Go_saccade_exp);
p_peak_time_Go_simvsdata = mc_perm_test(ts_saccade(cells_peak_idx_Go_saccade_sim)', ...
    ts_saccade(cells_peak_idx_Go_saccade_exp+idx_ts_0_exp-1)');
p_width_Go_amp_simvsdata = mc_perm_test(widths_cells_Go_saccade_sim, ...
    widths_cells_NC_saccade_exp);

p_peak_amp_NC_simvsdata = mc_perm_test(cells_peak_NC_saccade_sim, ...
    cells_peak_NC_saccade_exp);
p_peak_time_NC_simvsdata = mc_perm_test(ts_saccade(cells_peak_idx_NC_saccade_sim)', ...
    ts_saccade(cells_peak_idx_NC_saccade_exp+idx_ts_0_exp-1)');
p_width_NC_amp_simvsdata = mc_perm_test(widths_cells_NC_saccade_sim, ...
    widths_cells_NC_saccade_exp);

