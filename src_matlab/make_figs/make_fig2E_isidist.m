%% make Fig 2E - ISI distribution plots

clear
clc

rng('default')

%% load ISI distribution - experimental data

max_isi = 0.1; % sec, max value for isi

bins1 = 0:0.002:max_isi;
bins_middle1 = bins1(1:end-1) + diff(bins1)./2;
bins_middle1 = bins_middle1.*1e3; % sec -> ms

load(fullfile("data/monkeys_data/idx_L3neurons_sac.mat"), 'idx_4'); 

path_expData = fullfile('data', 'monkeys_data','isi_pretarget_errorBS'); %
% isidens is normalized by the number of trials
load(fullfile(path_expData, 'mean_isi_dist.mat'), 'isih_L3', ...
    'isih_L3_mean', 'isidens_L3')
isih_L3_preTarget_exp = isih_L3;
isih_L3_mean_preTarget_exp = isih_L3_mean; 
clearvars isih_L3 isih_L3_mean    

path_expData = fullfile('data', 'monkeys_data','isi_postsaccade_errorBS'); %
% isidens is normalized by the number of trials
load(fullfile(path_expData, 'mean_isi_dist.mat'), 'isih_L3_Go_mean', ...
    'isih_L3_Go', 'isih_L3_NC_mean', ...
    'isih_L3_NC', 'isidens_L3_Go', 'isidens_L3_NC')

mean_isi_preTarget_exp = mean(isih_L3_preTarget_exp(idx_4, :), 1,'omitnan');
sem_isi_preTarget_exp = std(isih_L3_preTarget_exp(idx_4, :), 0, 1, ...
    'omitnan')./sqrt(length(idx_4));

isidens_preTarget_exp = mean(isidens_L3(:,:, idx_4), 3, 'omitnan');

mean_isi_Go_exp = mean(isih_L3_Go(idx_4, :), 1,'omitnan');
sem_isi_Go_exp = std(isih_L3_Go(idx_4, :), 0, 1, ...
    'omitnan')./sqrt(length(idx_4));

isidens_Go_exp = mean(isidens_L3_Go(:,:, idx_4), 3, 'omitnan');

mean_isi_NC_exp = mean(isih_L3_NC(idx_4, :), 1,'omitnan');
sem_isi_NC_exp = std(isih_L3_NC(idx_4, :), 0, 1, ...
    'omitnan')./sqrt(length(idx_4));

isidens_NC_exp = mean(isidens_L3_NC(:,:, idx_4), 3, 'omitnan');

%% load ISI distribution - sim data

load(fullfile('data', 'sim_data','L3_isi_dist_sim.mat'), ...
    'isih_preTarget', 'isidens_preTarget', ...
    'isih_postSaccade_Go', 'isih_postSaccade_NC', ...
    'bins_middle', 'isidens_Go_postSaccade', 'isidens_NC_postSaccade')

load(fullfile('data', 'sim_data','L3_spk_tone_ft_struct_sim.mat'), ...
    'spikes_Go_target', 'spikes_NC_target', ...
    'spikes_Go_saccade', 'spikes_NC_saccade')

mean_isi_preTarget = mean(isih_preTarget, 1, 'omitnan');
sem_isi_preTarget = std(isih_preTarget, 0, 1, 'omitnan')./sqrt( ...
    length(spikes_Go_target.label));

mean_isi_Go = mean(isih_postSaccade_Go.avg./size(spikes_Go_target.trialtime, ...
    1), 1, 'omitnan');
sem_isi_Go = std(isih_postSaccade_Go.avg./size(spikes_Go_target.trialtime, 1), ...
    0, 1, 'omitnan')./sqrt(length(spikes_Go_target.label));

mean_isi_NC = mean(isih_postSaccade_NC.avg./size(spikes_Go_target.trialtime, 1), ...
    1, 'omitnan');
sem_isi_NC = std(isih_postSaccade_NC.avg./size(spikes_Go_target.trialtime, 1), ...
    0, 1, 'omitnan')./sqrt(length(spikes_Go_target.label));

%% fit exponential to mean ISI dist across cells
fitType = 'exp1';

%% experimental data 

[fitobject_preTarget1_exp, gof_pretarget_exp] = fit(bins_middle1(bins_middle1>=7.5)', ...
    mean(isih_L3_preTarget_exp(idx_4, (bins_middle1>=7.5)), ...
    1,'omitnan')', fitType);
fitobject_preTarget2_exp = fit(bins_middle1(bins_middle1>=7.5)', ...
    isih_L3_mean_preTarget_exp(bins_middle1>=7.5)', fitType);

[~, idx_max_Go_exp] = max(mean(isih_L3_Go(idx_4,:)));
[~, idx_max_NC_exp] = max(mean(isih_L3_NC(idx_4,:)));
[fitobject_Go_exp, gof_Go_exp] = fit(bins_middle1(idx_max_Go_exp:end)', ...
    mean(isih_L3_Go(idx_4,idx_max_Go_exp:end))', fitType);
[fitobject_NC_exp, gof_NC_exp] = fit(bins_middle1(idx_max_NC_exp:end)', ...
    mean(isih_L3_NC(idx_4,idx_max_NC_exp:end))', fitType);

%% simulated data

[~, idx_max_Go] = max(mean_isi_Go);
[~, idx_max_NC] = max(mean_isi_NC);

[fitobject_preTarget1, gof_pretarget] = fit(bins_middle(bins_middle>=7.5)', ...
    mean_isi_preTarget(bins_middle>=7.5)', fitType);

[fitobject_Go, gof_Go] = fit(bins_middle(idx_max_Go:end)', ...
    mean_isi_Go(idx_max_Go:end)', fitType);
[fitobject_NC, gof_NC] = fit(bins_middle(idx_max_NC:end)', ...
    mean_isi_NC(idx_max_NC:end)', fitType);

%% plot ISI(n) Vs ISI(n+1) 
font = 16;

%% ----------- Exp Data -----------

fig1 = figure('Units', 'inches','Position',[0 0 5 5]);
plot_isi_dist(fig1, isidens_preTarget_exp, mean_isi_preTarget_exp, bins_middle1, ...
    [], [], font)

fig2 = figure('Units', 'inches','Position',[0 0 5 5]);
minCMap = 0;
maxCMap = 0.0265; %4.7126;
plot_isi_dist(fig2, isidens_Go_exp, mean_isi_Go_exp, bins_middle1, minCMap, maxCMap, ...
    font)

fig3 = figure('Units', 'inches','Position',[0 0 5 5]);
plot_isi_dist(fig3, isidens_NC_exp, mean_isi_NC_exp, bins_middle1, minCMap, maxCMap, ...
    font)

%% ----------- Sim Data -----------

fig1 = figure('Units', 'inches','Position',[0 0 5 5]);
plot_isi_dist(fig1, isidens_preTarget, mean_isi_preTarget, bins_middle, ...
    [], [], font)

fig2 = figure('Units', 'inches','Position',[0 0 5 5]);
minCMap = 0;
maxCMap = 0.0574;
plot_isi_dist(fig2, isidens_Go_postSaccade, mean_isi_Go, bins_middle, ...
    minCMap, maxCMap, font)

fig3 = figure('Units', 'inches','Position',[0 0 5 5]);
plot_isi_dist(fig3, isidens_NC_postSaccade, mean_isi_NC, bins_middle, ...
    minCMap, maxCMap, font)

%% compare decay and intersect between cells pre vs post target (trial type)
%% --- fit exp to all cells
%% ---- exp data
fitobject_preTarget_allcells_exp = cellfun(@(x, y) myfit_fun(x, y), ...
    num2cell(repmat(bins_middle1(bins_middle1>=7.5), ...
    size(isih_L3_preTarget_exp(idx_4,:), 1), 1), 2), ...
    num2cell(isih_L3_preTarget_exp(idx_4, (bins_middle1>=7.5)),2), ...
    'UniformOutput', false);

fitobject_Go_allcells_exp = cellfun(@(x, y) myfit_fun(x, y), ...
    num2cell(repmat(bins_middle1(bins_middle1>=7.5), ...
    size(isih_L3_Go(idx_4,:), 1), 1), 2), ...
    num2cell(isih_L3_Go(idx_4, (bins_middle1>=7.5)),2), ...
    'UniformOutput', false);

fitobject_NC_allcells_exp = cellfun(@(x, y) myfit_fun(x, y), ...
    num2cell(repmat(bins_middle1(bins_middle1>=7.5), ...
    size(isih_L3_NC(idx_4,:), 1), 1), 2), ...
    num2cell(isih_L3_NC(idx_4, (bins_middle1>=7.5)),2), ...
    'UniformOutput', false);

b_preTarget_exp = cellfun(@(x) x.b, fitobject_preTarget_allcells_exp);
b_Go_exp = cellfun(@(x) x.b, fitobject_Go_allcells_exp);
b_NC_exp = cellfun(@(x) x.b, fitobject_NC_allcells_exp);

a_preTarget_exp = cellfun(@(x) x.a, fitobject_preTarget_allcells_exp);
a_Go_exp = cellfun(@(x) x.a, fitobject_Go_allcells_exp);
a_NC_exp = cellfun(@(x) x.a, fitobject_NC_allcells_exp);

%% 
b_exp = [b_preTarget_exp, b_Go_exp, b_NC_exp];
a_exp = [a_preTarget_exp, a_Go_exp, a_NC_exp];

%% ---- sim data

fitobject_preTarget_allcells = cellfun(@(x, y) myfit_fun_sim(x, y, [2.181, -0.03949]), ...
    num2cell(repmat(bins_middle(bins_middle>=7.5), ...
    size(isih_preTarget, 1), 1), 2), ...
    num2cell(isih_preTarget(:, (bins_middle>=7.5)),2), ...
    'UniformOutput', false);

[~, idx_max_target_cells] = max(isih_preTarget, [], 2);
[~, idx_max_Go_cells] = max(isih_postSaccade_Go.avg, [], 2);
[~, idx_max_NC_cells] = max(isih_postSaccade_NC.avg, [], 2);

bins = num2cell(repmat(bins_middle, ...
    size(isih_postSaccade_Go.avg, 1), 1), 2);

bins_target = cellfun(@(x, i) x(i:end), bins, num2cell(idx_max_target_cells, 2), ...
    'UniformOutput', false);
bins_Go = cellfun(@(x, i) x(i:end), bins, num2cell(idx_max_Go_cells, 2), ...
    'UniformOutput', false);
bins_NC = cellfun(@(x, i) x(i:end), bins, num2cell(idx_max_NC_cells, 2), ...
    'UniformOutput', false);

isi_Go_cells = cellfun(@(x, i) x(i:end), ...
    num2cell(isih_postSaccade_Go.avg./size( ...
    spikes_Go_target.trialtime, 1),2), num2cell(idx_max_Go_cells, 2), ...
    'UniformOutput', false);
isi_NC_cells = cellfun(@(x, i) x(i:end), ...
    num2cell(isih_postSaccade_NC.avg./size( ...
    spikes_NC_target.trialtime, 1),2), num2cell(idx_max_NC_cells, 2), ...
    'UniformOutput', false);

fitobject_Go_allcells = cellfun(@(x, y) myfit_fun_sim(x, y, [8.069, -0.07797]), ...
    bins_Go, isi_Go_cells, 'UniformOutput', false);
fitobject_NC_allcells = cellfun(@(x, y) myfit_fun_sim(x, y, [24.13, -0.1507]), ...
    bins_NC, isi_NC_cells, 'UniformOutput', false);

b_preTarget = cellfun(@(x) x.b, fitobject_preTarget_allcells);
b_Go = cellfun(@(x) x.b, fitobject_Go_allcells);
b_NC = cellfun(@(x) x.b, fitobject_NC_allcells);

y_preTarget = cellfun(@(cfun, point) feval(cfun,point), ...
    fitobject_preTarget_allcells, cellfun(@(x) x(1), bins_target, ...
    'UniformOutput', false));
y_Go = cellfun(@(cfun, point) feval(cfun,point), ...
    fitobject_Go_allcells, cellfun(@(x) x(1), bins_Go, 'UniformOutput', false));
y_NC = cellfun(@(cfun, point) feval(cfun,point), ...
    fitobject_NC_allcells, cellfun(@(x) x(1), bins_NC, 'UniformOutput', false));

%% plot cdf of b and y

[f_b_preTarget_exp,x_b_preTarget_exp] = ecdf(b_preTarget_exp);
[f_b_preTarget,x_b_preTarget] = ecdf(b_preTarget);
[f_b_Go_exp,x_b_Go_exp] = ecdf(b_Go_exp);
[f_b_Go,x_b_Go] = ecdf(b_Go);
[f_b_NC_exp,x_b_NC_exp] = ecdf(b_NC_exp);
[f_b_NC,x_b_NC] = ecdf(b_NC);

[f_a_preTarget_exp,x_a_preTarget_exp] = ecdf(a_preTarget_exp);
[f_y_preTarget,x_y_preTarget] = ecdf(y_preTarget);
[f_a_Go_exp,x_a_Go_exp] = ecdf(a_Go_exp);
[f_y_Go,x_y_Go] = ecdf(y_Go);
[f_a_NC_exp,x_a_NC_exp] = ecdf(a_NC_exp);
[f_y_NC,x_y_NC] = ecdf(y_NC);

%%
font = 10;

figure('Units', 'inches', 'Position', [0 0 6 3]);
t = tiledlayout(1, 4,'TileSpacing','Compact','Padding','Compact');

% exp decay
ii = 1;
nexttile(ii, [1 2])   
hold on;
plot(x_b_preTarget_exp, f_b_preTarget_exp, 'k-', 'linewidth', 1.5)
plot(x_b_preTarget, f_b_preTarget, 'r-', 'linewidth', 1.5)
ylabel('CDF')
% xlabel('Exponential Decay b (1/ms)')
title('pre-target')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

ii = ii + 2;
nexttile(ii, [1 2]);
hold on;
plot(x_b_Go_exp, f_b_Go_exp, 'k-', 'linewidth', 1.25)
plot(x_b_Go, f_b_Go, 'r-', 'linewidth', 1.25)
plot(x_b_NC_exp, f_b_NC_exp, 'k:', 'linewidth', 3)
plot(x_b_NC, f_b_NC, 'r:', 'linewidth', 3)
title('post-saccade')
legend({'Correct - Exp', 'Correct - Sim', 'Error - Exp', 'Error - Sim'}, ...
    'box', 'off', 'location', 'bestoutside')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

xlabel(t, 'Exponential Decay b (1/ms)','fontsize',font,'fontweight','bold')

%%
figure('Units', 'inches', 'Position', [0 0 6 3]);
t = tiledlayout(1, 4,'TileSpacing','Compact','Padding','Compact');

% norm num of spikes
ii = 1;
nexttile(ii, [1 2])   
hold on;
plot(x_a_preTarget_exp, f_a_preTarget_exp, 'k-', 'linewidth', 1.5)
plot(x_y_preTarget, f_y_preTarget, 'r-', 'linewidth', 1.5)
ylabel('CDF')
% title('pre-target')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

ii = ii + 2;
nexttile(ii, [1 2]);
hold on;
plot(x_a_Go_exp, f_a_Go_exp, 'k-', 'linewidth', 1.25)
plot(x_y_Go, f_y_Go, 'r-', 'linewidth', 1.25)
plot(x_a_NC_exp, f_a_NC_exp, 'k:', 'linewidth', 3)
plot(x_y_NC, f_y_NC, 'r:', 'linewidth', 3)
% title('post-saccade')
legend({'Correct - Exp', 'Correct - Sim', 'Error - Exp', 'Error - Sim'}, ...
    'box', 'off', 'location', 'bestoutside')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

xlabel(t, 'Max Norm Spk Count','fontsize',font,'fontweight','bold')

%% sub-functions

function f = myfit_fun(x, y)

f = fit(x', y', 'exp1');

end

function f = myfit_fun_sim(x, y, startpoint)

f = fit(x', y', 'exp1', 'StartPoint', startpoint);

end
