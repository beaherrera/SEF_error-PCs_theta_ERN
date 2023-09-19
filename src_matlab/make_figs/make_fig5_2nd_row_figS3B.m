%% analyze simulated LFP evoked by the activity of L3 and L5 pop of pyr error neurons

clear
clc

scaling = 27.48; % scaling factor to correct by the number of L3 and L5  error PCs in SEF

%% path to data

path2data = fullfile("data/sim_data/processed_lfps");

%% load data
load(fullfile(path2data, 'sim_lfps.mat'), 'ft_lfpGo', ...
    'ft_lfpNC', 'ft_lfpGo_L3', 'ft_lfpNC_L3', 'ft_lfpGo_L5', ...
    'ft_lfpNC_L5')

%% load events time | ms
trials_number = 1:20;

load(fullfile('data', 'sim_data', 'events_timing_Go.mat'))
saccade_times_Go = double(saccade_times(1,trials_number))';
target_times_Go = double(target_times(1,trials_number))';

load(fullfile('data', 'sim_data', 'events_timing_NC.mat'))
saccade_times_NC = double(saccade_times(1,trials_number))';
target_times_NC = double(target_times(1,trials_number))';

% convert event times to sec
target_times_secs_Go = double(target_times_Go).*1e-3; % sec, target times
saccade_times_secs_Go = double(saccade_times_Go).*1e-3; % sec, saccade times
target_times_secs_NC = double(target_times_NC).*1e-3; % sec, target times
saccade_times_secs_NC = double(saccade_times_NC).*1e-3; % sec, saccade times

%% maximum pre/post-event time windows

pre_saccade = 500; % ms before saccade
post_saccade = 1000; % ms after saccade

pre_target = 200;

%% depth information
Ne = 16; % number of electrodes in the shank
a = 1e-6; % [mm] position of the first electrode | a = 0 in the simulations
elec_spacing = 0.15; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % position of the electrodes
% along z
ze = ze.*1e3; % mm to um

%% center trials at saccade and calculate the avg lfp
select_window = @(x, time) x(:, (time - pre_saccade):(time + post_saccade));
select_bswindow = @(x, time) x(:,(time - pre_target):time);

cat_lfps = @(lfp_cell) cellfun(@(x, y) select_window(x,y), lfp_cell.trial', ...
    num2cell(target_times_Go + saccade_times_Go), 'UniformOutput', false);
cat_bslfp = @(lfp_cell) cellfun(@(x, y) select_bswindow(x,y), lfp_cell.trial', ...
    num2cell(target_times_Go), 'UniformOutput', false);

LFP_trials_Go = cat_lfps(ft_lfpGo);
LFP_3d_trials_Go = cat(3,LFP_trials_Go{:});

LFP_trials_NC = cat_lfps(ft_lfpNC);
LFP_3d_trials_NC = cat(3,LFP_trials_NC{:});

LFP_trials_Go_bs = cat_bslfp(ft_lfpGo);
LFP_3d_trials_Go_bs = cat(3,LFP_trials_Go_bs{:});

LFP_trials_NC_bs = cat_bslfp(ft_lfpNC);
LFP_3d_trials_NC_bs = cat(3,LFP_trials_NC_bs{:});

LFP_Go_saccade = mean(LFP_3d_trials_Go - mean(LFP_3d_trials_Go_bs, 2, ...
    'omitnan'), 3, 'omitnan');

LFP_NC_saccade = mean(LFP_3d_trials_NC - mean(LFP_3d_trials_NC_bs, 2, ...
    'omitnan'), 3, 'omitnan');

LFP_trials_Go_L3 = cat_lfps(ft_lfpGo_L3);
LFP_3d_trials_Go_L3 = cat(3,LFP_trials_Go_L3{:});

LFP_trials_NC_L3 = cat_lfps(ft_lfpNC_L3);
LFP_3d_trials_NC_L3 = cat(3,LFP_trials_NC_L3{:});

LFP_trials_Go_L3_bs = cat_bslfp(ft_lfpGo_L3);
LFP_3d_trials_Go_L3_bs = cat(3,LFP_trials_Go_L3_bs{:});

LFP_trials_NC_L3_bs = cat_bslfp(ft_lfpNC_L3);
LFP_3d_trials_NC_L3_bs = cat(3,LFP_trials_NC_L3_bs{:});

LFP_Go_L3_saccade = mean(LFP_3d_trials_Go_L3 - mean(LFP_3d_trials_Go_L3_bs, 2, ...
    'omitnan'), 3, 'omitnan');
LFP_NC_L3_saccade = mean(LFP_3d_trials_NC_L3 - mean(LFP_3d_trials_NC_L3_bs, 2, ...
    'omitnan'), 3, 'omitnan');

LFP_trials_Go_L5 = cat_lfps(ft_lfpGo_L5);
LFP_3d_trials_Go_L5 = cat(3,LFP_trials_Go_L5{:});

LFP_trials_NC_L5 = cat_lfps(ft_lfpNC_L5);
LFP_3d_trials_NC_L5 = cat(3,LFP_trials_NC_L5{:});

LFP_trials_Go_L5_bs = cat_bslfp(ft_lfpGo_L5);
LFP_3d_trials_Go_L5_bs = cat(3,LFP_trials_Go_L5_bs{:});

LFP_trials_NC_L5_bs = cat_bslfp(ft_lfpNC_L5);
LFP_3d_trials_NC_L5_bs = cat(3,LFP_trials_NC_L5_bs{:});

LFP_Go_L5_saccade = mean(LFP_3d_trials_Go_L5 - mean(LFP_3d_trials_Go_L5_bs, 2, ...
    'omitnan'), 3, 'omitnan');
LFP_NC_L5_saccade = mean(LFP_3d_trials_NC_L5 - mean(LFP_3d_trials_NC_L5_bs, 2, ...
    'omitnan'), 3, 'omitnan');

%% time
Fs_lfp = 1e3;
ts_saccade = (0:size(LFP_Go_saccade,2)-1) - pre_saccade; % msec

%% plot simulated LFPs
font = 9;

%% figure

LFPmax = 5e-3*max(abs(scaling.*[LFP_NC_saccade.*1e3 LFP_Go_saccade.*1e3]), ...
    [],'all'); % scale for the LFP plot

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
% -- correct
subplot(1, 3, 1)
hold on
for ii = 1:Ne
    plot(ts_saccade, scaling.*LFP_Go_saccade(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([300 3100])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error
subplot(1, 3, 2)
hold on
for ii = 1:Ne
    plot(ts_saccade,scaling.*LFP_NC_saccade(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
plot([1 1].*-250, [0 0.02]./LFPmax + (ze(16)) + 600, '-k', 'LineWidth', 1.5)
text(-190, median([0 0.02]./LFPmax + (ze(16)) + 600), '20\muV', ...
    'FontSize', font,'HorizontalAlignment','left')
ax = gca;
ax.FontSize = font;
ylim([300 3100])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error - correct
subplot(1,3,3)
LFPmax = 5e-3*max(abs(scaling.*LFP_NC_saccade.*1e3 - ...
    scaling.*LFP_Go_saccade.*1e3),[],'all'); % scale for the LFP plot
hold on
for ii = 1:Ne
    plot(ts_saccade,(scaling.*LFP_NC_saccade(Ne-ii+1,:).*1e3 - ...
        scaling.*LFP_Go_saccade(Ne-ii+1,:).*1e3)./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([300 3100])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.5)

%% CSD analysis
%% -- Parameters
el_pos = ze*1e-6; % um -> m
cond = 0.33; %[S/m] gray matter conductance
cond_top = 0; %[S/m] conductance at the top (cylinder)
gauss_sigma = 0.1*1e-3;   %[m] Gaussian filter std
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
diam = 3e-3; % [m] cylinder diameter

%% - solve Pettersen model
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);

[zs_Go_saccade,CSD_cs_Go_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_Go_saccade.*1e-3,Fcs);
[zs_NC_saccade,CSD_cs_NC_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_NC_saccade.*1e-3,Fcs);

[zs_Go_L3_saccade,CSD_cs_Go_L3_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_Go_L3_saccade.*1e-3,Fcs);
[zs_NC_L3_saccade,CSD_cs_NC_L3_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_NC_L3_saccade.*1e-3,Fcs);

[zs_Go_L5_saccade,CSD_cs_Go_L5_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_Go_L5_saccade.*1e-3,Fcs);
[zs_NC_L5_saccade,CSD_cs_NC_L5_saccade] = make_cubic_splines(el_pos, ...
    scaling.*LFP_NC_L5_saccade.*1e-3,Fcs);

if ~isempty(gauss_sigma) && gauss_sigma~=0 %filter iCSD
    
    [zs_Go_saccade,CSD_cs_Go_saccade] = gaussian_filtering(zs_Go_saccade, ...
        CSD_cs_Go_saccade, gauss_sigma, filter_range);
    [zs_NC_saccade,CSD_cs_NC_saccade] = gaussian_filtering(zs_NC_saccade, ...
        CSD_cs_NC_saccade, gauss_sigma, filter_range);
    
    [zs_Go_L3_saccade,CSD_cs_Go_L3_saccade] = gaussian_filtering(zs_Go_L3_saccade, ...
        CSD_cs_Go_L3_saccade, gauss_sigma, filter_range);
    [zs_NC_L3_saccade,CSD_cs_NC_L3_saccade] = gaussian_filtering(zs_NC_L3_saccade, ...
        CSD_cs_NC_L3_saccade, gauss_sigma, filter_range);
    
    [zs_Go_L5_saccade,CSD_cs_Go_L5_saccade] = gaussian_filtering(zs_Go_L5_saccade, ...
        CSD_cs_Go_L5_saccade, gauss_sigma, filter_range);
    [zs_NC_L5_saccade,CSD_cs_NC_L5_saccade] = gaussian_filtering(zs_NC_L5_saccade, ...
        CSD_cs_NC_L5_saccade, gauss_sigma, filter_range);
end

iCSD_Go_saccade = CSD_cs_Go_saccade; % [nA/mm3] current source density
iCSD_NC_saccade = CSD_cs_NC_saccade; % [nA/mm3] current source density
iCSD_Go_L3_saccade = CSD_cs_Go_L3_saccade; % [nA/mm3] current source density
iCSD_NC_L3_saccade = CSD_cs_NC_L3_saccade; % [nA/mm3] current source density
iCSD_Go_L5_saccade = CSD_cs_Go_L5_saccade; % [nA/mm3] current source density
iCSD_NC_L5_saccade = CSD_cs_NC_L5_saccade; % [nA/mm3] current source density

zs_Go_saccade = zs_Go_saccade*1e3; % m -> mm
zs_NC_saccade = zs_NC_saccade*1e3; % m -> mm

%% - saving CSD matrix
% save('data\sim_data\simCSD.mat', 'zs_Go_saccade', 'zs_NC_saccade', ...
%     'iCSD_Go_saccade', 'iCSD_NC_saccade')

%% --- Plotting the CSD
h = 150; % um, inter-electrode distance
ze_plot = 0:h:(Ne-1)*h;
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

%% L3 + L5
max_CSD = max(abs([iCSD_Go_saccade iCSD_NC_saccade]), [], 'all');
font = 12;

figure('Units', 'inches','Position',[0.05 0.05 5.2 3]);
tiledlayout(1, 2,'TileSpacing','Compact','Padding','Compact');

nexttile(1); % subplot(1,3,1)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_Go_saccade);
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

nexttile(2); %subplot(1,3,2)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_saccade);
colormap(flip(jet));
max_CSD = max(abs(iCSD_NC_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
c = colorbar;
c.Label.FontSize = font;
c.Layout.Tile = 'east';

%%
figure('Units', 'inches','Position',[0.05 0.05 2.7 3]);
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_saccade - iCSD_Go_saccade);
colorbar;
[map,~,~,~] = brewermap(256,'PiYG');
colormap(gca, (map));
max_CSD = max(abs(iCSD_NC_saccade - iCSD_Go_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% L3
font = 10;

max_CSD = max(abs([iCSD_Go_L3_saccade iCSD_NC_L3_saccade ...
    iCSD_Go_L5_saccade iCSD_NC_L5_saccade]), [], 'all');

figure('Units', 'inches','Position',[0.05 0.05 5 2]);
tiledlayout(1, 3,'TileSpacing','Compact','Padding','Compact');

nexttile(1); % subplot(1,3,1)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_Go_L3_saccade);
colormap(flip(jet));
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
xticklabels({})
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

nexttile(2); %subplot(1,3,2)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_L3_saccade);
colormap(flip(jet));
max_CSD = max(abs(iCSD_NC_L3_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
% xticklabels({})
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
% c = colorbar;
% c.Label.FontSize = font;

nexttile(3)
% ax_diff = subplot(1,3,3);
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_saccade - iCSD_Go_saccade);
[map,~,~,~] = brewermap(256,'PiYG');
colormap(gca, (map));
max_CSD = max(abs(iCSD_NC_saccade - iCSD_Go_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
xticklabels({})
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% L5

figure('Units', 'inches','Position',[0.05 0.05 5 2]);
tiledlayout(1, 3,'TileSpacing','Compact','Padding','Compact');

nexttile(1); % subplot(1,3,1)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_Go_L5_saccade);
colormap(flip(jet));
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
xticklabels({})
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

nexttile(2); %subplot(1,3,2)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_L5_saccade);
colormap(flip(jet));
max_CSD = max(abs(iCSD_NC_L5_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

nexttile(3)
imagesc(ts_saccade, zs_Go_saccade.*1e3, iCSD_NC_saccade - iCSD_Go_saccade);
[map,~,~,~] = brewermap(256,'PiYG');
colormap(gca, (map));
max_CSD = max(abs(iCSD_NC_saccade - iCSD_Go_saccade), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
xticks([-500 0 500 1000])
xticklabels({})
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
