%% create Fig 1 B - spiking rate plot

clear
clc

%% load spk rate data

path2data = fullfile('data', 'monkeys_data','spkrate_errorneurons');
load(fullfile(path2data, 'fr_error_cells_saccade.mat'), 'mean_fr_L*_NC', ...
    'mean_fr_L*_Go', 'ts_sess')

ts_saccade = ts_sess{1}; 

load(fullfile("data/monkeys_data/idx_L3neurons_sac.mat"), 'idx_4'); 
load(fullfile("data/monkeys_data/idx_L5neurons_sac.mat"), 'idx'); 

cells_L3_mean_fr_NC_saccade = mean_fr_L3_NC(idx_4, :);
cells_L3_mean_fr_Go_saccade = mean_fr_L3_Go(idx_4, :);
cells_L5_mean_fr_NC_saccade = mean_fr_L5_NC(idx, :);
cells_L5_mean_fr_Go_saccade = mean_fr_L5_Go(idx, :);

%% L3
font = 10;

cell_num_L3 = 3;

figure('Units', 'inches','Position',[0.05 0.05 4 2]);
hold on;
% old color: 'b'
p2 = plot(ts_saccade.*1e3, cells_L3_mean_fr_NC_saccade(cell_num_L3,:) ...
    , ':', 'Color', 'k', ...
    'linewidth',3);

% old color: 'k'
hold on;
p1 = plot(ts_saccade.*1e3, cells_L3_mean_fr_Go_saccade(cell_num_L3,:) ...
    , 'Color', 'k', ...
    'linewidth',1.25);

xlim([-.5 .5].*1e3)
ylim([10 70])
box off
legend([p1 p2], {'Correct', 'Error'}, 'Box','off', 'Location','northwest')
xlabel('Time from saccade (ms)')
ylabel({'Relative firing', 'rate (spks/s)'})
set(gca, 'YTickLabel',[],'YTick',[],'YColor','None', 'LineWidth', 2, ...
    'FontSize', 10, 'FontWeight', 'bold')

%% L5

cell_num_L5 = 4;

figure('Units', 'inches','Position',[0.05 0.05 4 2]);
hold on;
% old color: 'b'
p2 = plot(ts_saccade.*1e3, cells_L5_mean_fr_NC_saccade(cell_num_L5,:) ...
    , ':', 'Color', 'k', ...
    'linewidth',3);

% old color: 'k'
hold on;
p1 = plot(ts_saccade.*1e3, cells_L5_mean_fr_Go_saccade(cell_num_L5,:) ...
    , 'Color', 'k', ...
    'linewidth',1.25);
hold on; 
plot([1 1].*-500, [15 35], '-k', 'LineWidth',2)

xlim([-.5 .5].*1e3)
ylim([5 65])
box off
% legend([p1 p2], {'Correct', 'Error'}, 'Box','off', 'Location','northwest')
xlabel('Time from saccade (ms)')
ylabel({'Relative firing', 'rate (spks/s)'})
set(gca, 'YTickLabel',[],'YTick',[],'YColor','None', 'LineWidth', 2, ...
    'FontSize', 10, 'FontWeight', 'bold')
