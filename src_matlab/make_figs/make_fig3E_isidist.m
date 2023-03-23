%% make Fig 3E - ISI distributions

clear
clc

rng('default')

%% load ISI distribution - experimental data

max_isi = 0.1; % sec, max value for isi

bins1 = 0:0.002:max_isi;
bins_middle1 = bins1(1:end-1) + diff(bins1)./2;
bins_middle1 = bins_middle1.*1e3; % sec -> ms

load(fullfile("data/monkeys_data/idx_L5neurons_sac.mat"), 'idx'); 

path_expData = fullfile('data', 'monkeys_data','isi_pretarget_errorBS'); %
% isidens is normalized by the number of trials
load(fullfile(path_expData, 'mean_isi_dist.mat'), 'isih_L5', ...
    'isih_L5_mean', 'isidens_L5')
isih_L5_preTarget_exp = isih_L5;
isih_L5_mean_preTarget_exp = isih_L5_mean; 
clearvars isih_L5 isih_L5_mean    

path_expData = fullfile('data', 'monkeys_data','isi_postsaccade_errorBS'); %
% isidens is normalized by the number of trials
load(fullfile(path_expData, 'mean_isi_dist.mat'), 'isih_L5_Go_mean', ...
    'isih_L5_Go', 'isih_L5_NC_mean', ...
    'isih_L5_NC', 'isidens_L5_Go', 'isidens_L5_NC')

mean_isi_preTarget_exp = mean(isih_L5_preTarget_exp(idx, :), 1,'omitnan');
sem_isi_preTarget_exp = std(isih_L5_preTarget_exp(idx, :), 0, 1, ...
    'omitnan')./sqrt(length(idx));

isidens_preTarget_exp = mean(isidens_L5(:,:, idx), 3, 'omitnan');

mean_isi_Go_exp = mean(isih_L5_Go(idx, :), 1,'omitnan');
sem_isi_Go_exp = std(isih_L5_Go(idx, :), 0, 1, ...
    'omitnan')./sqrt(length(idx));

isidens_Go_exp = mean(isidens_L5_Go(:,:, idx), 3, 'omitnan');

mean_isi_NC_exp = mean(isih_L5_NC(idx, :), 1,'omitnan');
sem_isi_NC_exp = std(isih_L5_NC(idx, :), 0, 1, ...
    'omitnan')./sqrt(length(idx));

isidens_NC_exp = mean(isidens_L5_NC(:,:, idx), 3, 'omitnan');

%% load ISI distribution - sim data

load(fullfile('data', 'sim_data','L5_isi_dist_sim.mat'), ...
    'isih_preTarget', 'isidens_preTarget', ...
    'isih_postSaccade_Go', 'isih_postSaccade_NC', ...
    'bins_middle', 'isidens_Go_postSaccade', 'isidens_NC_postSaccade')

load(fullfile('data', 'sim_data','L5_spk_tone_ft_struct_sim.mat'), ...
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
fitType = 'exp2';

%% experimental data 

[fitobject_preTarget1_exp, gof_preTarget1_exp] = fit(bins_middle1(2:end)', ...
    mean_isi_preTarget((2:end))', 'poly1');

[fitobject_Go_exp, gof_Go_exp] = fit(bins_middle1(2:end)', ...
    mean_isi_Go_exp(2:end)', fitType);

[fitobject_NC_exp, gof_NC_exp] = fit(bins_middle1(2:end)', ...
    mean_isi_NC_exp(2:end)', fitType);

%% simulated data

[fitobject_preTarget1, gof_preTarget1] = fit(bins_middle', ...
    mean_isi_preTarget', 'poly1');

[fitobject_Go, gof_Go] = fit(bins_middle', ...
    mean_isi_Go', fitType, ...
    'StartPoint', [-870.6, -0.06919, 870.7, -0.06912]);
[fitobject_NC, gof_NC] = fit(bins_middle', ...
    mean_isi_NC', fitType, ...
    'StartPoint', [-870.6, -0.06919, 870.7, -0.06912]);

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

fitobject_Go_allcells_exp = cellfun(@(x, y) myfit_fun(x, y), ...
    num2cell(repmat(bins_middle1(2:end), ...
    size(isih_L5_Go(idx, :), 1), 1), 2), ...
    num2cell(isih_L5_Go(idx, (2:end)),2), ...
    'UniformOutput', false);

fitobject_NC_allcells_exp = cellfun(@(x, y) myfit_fun(x, y), ...
    num2cell(repmat(bins_middle1(2:end), ...
    size(isih_L5_NC(idx, :), 1), 1), 2), ...
    num2cell(isih_L5_NC(idx, (2:end)),2), ...
    'UniformOutput', false);

b_Go_exp = cellfun(@(x) x.b, fitobject_Go_allcells_exp);
b_NC_exp = cellfun(@(x) x.b, fitobject_NC_allcells_exp);

d_Go_exp = cellfun(@(x) x.d, fitobject_Go_allcells_exp);
d_NC_exp = cellfun(@(x) x.d, fitobject_NC_allcells_exp);

b_exp = [b_Go_exp, b_NC_exp];
d_exp = [d_Go_exp, d_NC_exp];

%% ---- sim data

fitobject_Go_allcells = cellfun(@(x, y) myfit_fun_sim(x, y, ...
    [-870.6, -0.06919, 870.7, -0.06912]), ...
    num2cell(repmat(bins_middle, ...
    size(isih_postSaccade_Go.avg, 1), 1), 2), ...
    num2cell(isih_postSaccade_Go.avg./size( ...
    spikes_Go_target.trialtime, 1),2), ...
    'UniformOutput', false);

fitobject_NC_allcells = cellfun(@(x, y) myfit_fun_sim(x, y, ...
    [-870.6, -0.06919, 870.7, -0.06912]), ...
    num2cell(repmat(bins_middle, ...
    size(isih_postSaccade_NC.avg, 1), 1), 2), ...
    num2cell(isih_postSaccade_NC.avg./size( ...
    spikes_Go_target.trialtime, 1),2), ...
    'UniformOutput', false);

for ii=1:length(fitobject_NC_allcells)

    figure;
    subplot(2,1,1)
    hold on
    plot(bins_middle, isih_postSaccade_Go.avg(ii,:)./size( ...
    spikes_Go_target.trialtime, 1), '-r', ...
        'LineWidth', 1.5);
    plot(bins_middle, ...
        feval(fitobject_Go_allcells{1,1},bins_middle), '-k', ...
        'LineWidth', 1.5);

    subplot(2,1,2)
    hold on
    plot(bins_middle, isih_postSaccade_NC.avg(ii,:)./size( ...
    spikes_NC_target.trialtime, 1), ':k', ...
        'LineWidth', 1.5);
    plot(bins_middle, ...
        feval(fitobject_NC_allcells{1,1},bins_middle), '-k', ...
        'LineWidth', 1.5);
end

b_Go = cellfun(@(x) x.b, fitobject_Go_allcells);
b_NC = cellfun(@(x) x.b, fitobject_NC_allcells);

d_Go = cellfun(@(x) x.d, fitobject_Go_allcells);
d_NC = cellfun(@(x) x.d, fitobject_NC_allcells);

b = [b_Go, b_NC];
d = [d_Go, d_NC];


%% plot cdf of b and y

[f_b_Go_exp,x_b_Go_exp] = ecdf(b_Go_exp);
[f_b_Go,x_b_Go] = ecdf(b_Go);
[f_b_NC_exp,x_b_NC_exp] = ecdf(b_NC_exp);
[f_b_NC,x_b_NC] = ecdf(b_NC);

[f_d_Go_exp,x_d_Go_exp] = ecdf(d_Go_exp);
[f_d_Go,x_d_Go] = ecdf(d_Go);
[f_d_NC_exp,x_d_NC_exp] = ecdf(d_NC_exp);
[f_d_NC,x_d_NC] = ecdf(d_NC);

%%
font = 10;

figure('Units', 'inches', 'Position', [0 0 6 3]);
t = tiledlayout(1, 4,'TileSpacing','Compact','Padding','Compact');

% exp decay
ii = 1;
nexttile(ii, [1 2])
hold on;
plot(x_b_Go_exp, f_b_Go_exp, 'k-', 'linewidth', 1.25)
plot(x_b_Go, f_b_Go, 'r-', 'linewidth', 1.25)
plot(x_b_NC_exp, f_b_NC_exp, 'k:', 'linewidth', 3)
plot(x_b_NC, f_b_NC, 'r:', 'linewidth', 3)
title('b')
% legend({'Correct - Exp', 'Correct - Sim', 'Error - Exp', 'Error - Sim'}, ...
%     'box', 'off', 'location', 'bestoutside')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

xlabel(t, 'Exponential Decay b (1/ms)','fontsize',font,'fontweight','bold')

ii = ii + 2;
nexttile(ii, [1 2]);
hold on;
plot(x_d_Go_exp, f_d_Go_exp, 'k-', 'linewidth', 1.25)
plot(x_d_Go, f_d_Go, 'r-', 'linewidth', 1.25)
plot(x_d_NC_exp, f_d_NC_exp, 'k:', 'linewidth', 3)
plot(x_d_NC, f_d_NC, 'r:', 'linewidth', 3)
title('d')
legend({'Correct - Exp', 'Correct - Sim', 'Error - Exp', 'Error - Sim'}, ...
    'box', 'off', 'location', 'bestoutside')
set(gca, 'box', 'off','linewidth',1,'fontsize',font,'fontweight','bold')

xlabel(t, 'Exp Time Constant (1/ms)','fontsize',font,'fontweight','bold')

%% sub-functions

function f = myfit_fun(x, y)

f = fit(x', y', 'exp2');

end

function f = myfit_fun_sim(x, y, startpoint)

f = fit(x', y', 'exp2', 'StartPoint', startpoint);

end
