%% create fig 5 top panels and Fig S3A

clear
clc

%% perpendicular sessions number

Sess_EuP1 = 14:19; % Eu, site P1
Sess_XP2P3 = 20:29; % X, 20-25 site P2, 26-29 site P3

SessNumb = [Sess_EuP1, Sess_XP2P3];

%% load data
load("data\monkeys_data\lfp_csd\avg_lfp_csd_sess.mat")
load("data\monkeys_data\lfp_csd\grand_avg_lfp_csd.mat", ...
    'tspan', 'zs_NC', 'zs_Go', 'grand_avg_LFP_NC', 'grand_avg_LFP_Go', ...
    'grand_avg_iCSD_NC', 'grand_avg_iCSD_Go')

%% grand avg lfp per monkey

avg_LFP_sess_Go = structfun(@(x) x.avg_LFP_Go, avg_LFP, ...
    'UniformOutput',false);
avg_LFP_sess_Go = struct2cell(avg_LFP_sess_Go);
avg_LFP_sess_Go = cat(3, avg_LFP_sess_Go{:});
grand_avg_LFP_Go_Eu = mean(avg_LFP_sess_Go(:,:,ismember(SessNumb, Sess_EuP1)), ...
    3, 'omitnan');
grand_avg_LFP_Go_X = mean(avg_LFP_sess_Go(:,:,ismember(SessNumb, Sess_XP2P3)), ...
    3, 'omitnan');

avg_LFP_sess_NC = structfun(@(x) x.avg_LFP_NC, avg_LFP, ...
    'UniformOutput',false);
avg_LFP_sess_NC = struct2cell(avg_LFP_sess_NC);
avg_LFP_sess_NC = cat(3, avg_LFP_sess_NC{:});
grand_avg_LFP_NC_Eu = mean(avg_LFP_sess_NC(:,:,ismember(SessNumb, Sess_EuP1)), ...
    3, 'omitnan');
grand_avg_LFP_NC_X = mean(avg_LFP_sess_NC(:,:,ismember(SessNumb, Sess_XP2P3)), ...
    3, 'omitnan');

%% calculate ind monkeys CSD
% -- Parameters CSD calculation
Ne = 16; % number of electrodes in the shank
a = 1e-6; % [mm] position of the first electrode
elec_spacing = 0.15; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % position of the electrodes
% along z
el_pos = ze*1e-3;  % [m] electrode positions with respect to the pia surface
cond = 0.33; %[S/m] gray matter conductance 0.4
cond_top = cond; % conductance outside the gray matter (above pia matter)
this_tol = 1e-6;  % tolerance
gauss_sigma = 0.1e-3;   %[m] Gaussian filter std
diam = 3e-3; % [m] cylinder diameter
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent

% -- calculate cubic splines
Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);

% -- get CSD
[zs_Go_Eu,CSD_cs_Go_Eu] = make_cubic_splines(el_pos, ...
    grand_avg_LFP_Go_Eu, Fcs);
[zs_Go_X,CSD_cs_Go_X] = make_cubic_splines(el_pos, ...
    grand_avg_LFP_Go_X, Fcs);
if ~isempty(gauss_sigma) && gauss_sigma~=0 %filter iCSD
    [zs_Go_Eu,CSD_cs_Go_Eu] = gaussian_filtering(zs_Go_Eu, ...
        CSD_cs_Go_Eu, gauss_sigma, filter_range);
    [zs_Go_X,CSD_cs_Go_X] = gaussian_filtering(zs_Go_X, ...
        CSD_cs_Go_X, gauss_sigma, filter_range);
end
grand_avg_iCSD_Go_Eu = CSD_cs_Go_Eu; % [nA/mm3] current source density
zs_Go_Eu = zs_Go_Eu*1e3; % [mm]
grand_avg_iCSD_Go_X = CSD_cs_Go_X; % [nA/mm3] current source density
zs_Go_X = zs_Go_X*1e3; % [mm]

[zs_NC_Eu,CSD_cs_NC_Eu] = make_cubic_splines(el_pos, ...
    grand_avg_LFP_NC_Eu, Fcs);
[zs_NC_X,CSD_cs_NC_X] = make_cubic_splines(el_pos, ...
    grand_avg_LFP_NC_X, Fcs);
if ~isempty(gauss_sigma) && gauss_sigma~=0 %filter iCSD
    [zs_NC_Eu,CSD_cs_NC_Eu] = gaussian_filtering(zs_NC_Eu, ...
        CSD_cs_NC_Eu, gauss_sigma, filter_range);
    [zs_NC_X,CSD_cs_NC_X] = gaussian_filtering(zs_NC_X, ...
        CSD_cs_NC_X, gauss_sigma, filter_range);
end
grand_avg_iCSD_NC_Eu = CSD_cs_NC_Eu; % [nA/mm3] current source density
zs_NC_Eu = zs_NC_Eu*1e3; % [mm]
grand_avg_iCSD_NC_X = CSD_cs_NC_X; % [nA/mm3] current source density
zs_NC_X = zs_NC_X*1e3; % [mm]

%% depth information
zs_Go = zs_Go.*1e3; % mm to um
zs_NC = zs_NC.*1e3; % mm to um

Ne = 16; % number of electrodes in the shank
a = 0; % [mm] position of the first electrode
elec_spacing = 0.15; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % position of the electrodes
% along z
ze = ze.*1e3; % mm to um

%% plot LFP
font = 9;

%% figure
%% avg across monkeys

LFPmax = 5e-3*max(abs([grand_avg_LFP_NC.*1e3 ...
    grand_avg_LFP_Go.*1e3]),[],'all'); % scale for the LFP plot

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
% -- correct
subplot(1, 3, 1)
hold on
for ii = 1:Ne
    plot(tspan, grand_avg_LFP_Go(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Correct')
box 'off'
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error
subplot(1, 3, 2)
hold on
for ii = 1:Ne
    plot(tspan,grand_avg_LFP_NC(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
plot([1 1].*-250, [0 0.05]./LFPmax + (ze(16)) + 600, '-k', 'LineWidth', 1.5)
text(-190, median([0 0.05]./LFPmax + (ze(16)) + 600), '50\muV', ...
    'FontSize', font,'HorizontalAlignment','left')

ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
xlabel('Time from Saccade (s)')
% title('Error')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error - correct
subplot(1,3,3)
LFPmax = 5e-3*max(abs(grand_avg_LFP_NC.*1e3 - ...
    grand_avg_LFP_Go.*1e3),[],'all'); % scale for the LFP plot
hold on
for ii = 1:Ne
    plot(tspan,(grand_avg_LFP_NC(Ne-ii+1,:).*1e3 - ...
        grand_avg_LFP_Go(Ne-ii+1,:).*1e3)./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Error - Correct')
box 'off'
% xlabel('Time from saccade onset (ms)')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.5)

%% Eu

LFPmax = 5e-3*max(abs([grand_avg_LFP_NC_Eu.*1e3 ...
    grand_avg_LFP_Go_Eu.*1e3]),[],'all'); % scale for the LFP plot

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
% -- correct
subplot(1, 3, 1)
hold on
for ii = 1:Ne
    plot(tspan, grand_avg_LFP_Go_Eu(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Correct')
box 'off'
% xlabel('Time from saccade onset (ms)')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error
subplot(1, 3, 2)
hold on
for ii = 1:Ne
    plot(tspan,grand_avg_LFP_NC_Eu(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
plot([1 1].*-250, [0 0.05]./LFPmax + (ze(16)) + 600, '-k', 'LineWidth', 1.5)
text(-190, median([0 0.05]./LFPmax + (ze(16)) + 600), '50\muV', ...
    'FontSize', font,'HorizontalAlignment','left')
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
xlabel('Time from Saccade (s)')
% title('Error')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error - correct
subplot(1,3,3)
LFPmax = 5e-3*max(abs(grand_avg_LFP_NC_Eu.*1e3 - ...
    grand_avg_LFP_Go_Eu.*1e3),[],'all'); % scale for the LFP plot
hold on
for ii = 1:Ne
    plot(tspan,(grand_avg_LFP_NC_Eu(Ne-ii+1,:).*1e3 - ...
        grand_avg_LFP_Go_Eu(Ne-ii+1,:).*1e3)./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Error - Correct')
box 'off'
% xlabel('Time from saccade onset (ms)')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.5)

%% X

LFPmax = 5e-3*max(abs([grand_avg_LFP_NC_X.*1e3 ...
    grand_avg_LFP_Go_X.*1e3]),[],'all'); % scale for the LFP plot

figure('Units', 'inches','Position',[0.05 0.05 4 2.5]);
% -- correct
subplot(1, 3, 1)
hold on
for ii = 1:Ne
    plot(tspan, grand_avg_LFP_Go_X(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Correct')
box 'off'
% xlabel('Time from saccade onset (ms)')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error
subplot(1, 3, 2)
hold on
for ii = 1:Ne
    plot(tspan,grand_avg_LFP_NC_X(Ne-ii+1,:).*1e3./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
plot([1 1].*-250, [0 0.02]./LFPmax + (ze(16)) + 600, '-k', 'LineWidth', 1.5)
text(-190, median([0 0.02]./LFPmax + (ze(16)) + 600), '20\muV', ...
    'FontSize', font,'HorizontalAlignment','left')
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
box 'off'
xlabel('Time from Saccade (s)')
% title('Error')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.25)

% -- error - correct
subplot(1,3,3)
LFPmax = 5e-3*max(abs(grand_avg_LFP_NC_X.*1e3 - ...
    grand_avg_LFP_Go_X.*1e3),[],'all'); % scale for the LFP plot
hold on
for ii = 1:Ne
    plot(tspan,(grand_avg_LFP_NC_X(Ne-ii+1,:).*1e3 - ...
        grand_avg_LFP_Go_X(Ne-ii+1,:).*1e3)./LFPmax + ...
        (ze(ii)) + 500, ...
        'color', 'k', 'clipping','on', 'linewidth', 1)
end
ax = gca;
ax.FontSize = font;
ylim([200 3200])
xlim([-500 500])
xticks([-500 0 500 1000])
% title('Error - Correct')
box 'off'
% xlabel('Time from saccade onset (ms)')
set(ax, 'YTickLabel',[],'YTick',[],'YColor','None',...
    'FontSize',font, 'FontWeight', 'bold','LineWidth',1.5)

%% plot CSD
font = 10;
h = 150; % um, inter-electrode distance
ze_plot = 0:h:(Ne-1)*h;
z_depth  = [1125:-150:0, 0, -75:-150:-1125]; % mm depth

max_CSD = max(abs([grand_avg_iCSD_Go grand_avg_iCSD_NC]), [], 'all');

figure('Units', 'inches','Position',[0.05 0.05 10 3]);
subplot(1,3,1)
imagesc(tspan, zs_Go, grand_avg_iCSD_Go);
colormap(flip(jet));
c = colorbar;
% c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
yticks(sort([ze 1125]))
ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
% xlabel('Time relative to saccade (ms)', 'FontSize',font);
% ylabel('Cortical Depth (\mum)','FontSize',font);
% title('Correct')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

subplot(1,3,2)
imagesc(tspan, zs_Go, grand_avg_iCSD_NC);
colormap(flip(jet));
c = colorbar;
% c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
caxis([bar_min bar_max]);
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
% xlabel('Time from Saccade (s)', 'FontSize',font);
% ylabel('z (\mum)','FontSize',font);
% title('Error')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

ax_diff = subplot(1,3,3);
imagesc(tspan, zs_Go, grand_avg_iCSD_NC - grand_avg_iCSD_Go);
c = colorbar;
% c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
[map,~,~,~] = brewermap(256,'PiYG');
colormap(ax_diff, (map));
max_CSD = max(abs(grand_avg_iCSD_NC - grand_avg_iCSD_Go), [], 'all');
bar_min = -max_CSD;
bar_max = max_CSD;
caxis([bar_min bar_max]);
xlim([-500 500])
yticks(sort([ze 1125]))
ytick_labels = cellfun(@(x) ' ', num2cell(z_depth), 'UniformOutput', false);
yticklabels(ytick_labels)
ax = gca;
ax.FontSize = font;
ax.YDir='reverse';
% xlabel('Time relative to saccade (ms)', 'FontSize',font);
% ylabel('z (\mum)','FontSize',font);
% title('Error - Correct')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% Eu

font = 12;

max_CSD = max(abs([grand_avg_iCSD_Go_Eu grand_avg_iCSD_NC_Eu]), [], 'all');

figure('Units', 'inches','Position',[0.05 0.05 5 3]);
tiledlayout(1, 2,'TileSpacing','Compact','Padding','Compact');

nexttile(1);
imagesc(tspan, zs_Go, grand_avg_iCSD_Go_Eu);
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

nexttile(2);
imagesc(tspan, zs_Go, grand_avg_iCSD_NC_Eu);
colormap(flip(jet));
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
c = colorbar;
c.Label.FontSize = font;
c.Layout.Tile = 'east';

%%
figure('Units', 'inches','Position',[0.05 0.05 2.7 3]);
imagesc(tspan, zs_Go, grand_avg_iCSD_NC_Eu - grand_avg_iCSD_Go_Eu);
c = colorbar;
c.Label.FontSize = font;
[map,~,~,~] = brewermap(256,'PiYG');
colormap(gca, (map));
max_CSD = max(abs(grand_avg_iCSD_NC_Eu - grand_avg_iCSD_Go_Eu), [], 'all');
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

%% X

max_CSD = max(abs([grand_avg_iCSD_Go_X grand_avg_iCSD_NC_X]), [], 'all');

figure('Units', 'inches','Position',[0.05 0.05 10 3]);
subplot(1,3,1)
imagesc(tspan, zs_Go, grand_avg_iCSD_Go_X);
colormap(flip(jet));
c = colorbar;
% c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
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
% xlabel('Time relative to saccade (ms)', 'FontSize',font);
% ylabel('Cortical Depth (\mum)','FontSize',font);
% title('Correct')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

subplot(1,3,2)
imagesc(tspan, zs_Go, grand_avg_iCSD_NC_X);
colormap(flip(jet));
c = colorbar;
c.Label.FontSize = font;
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
% xlabel('Time from Saccade (s)', 'FontSize',font);
% ylabel('z (\mum)','FontSize',font);
% title('Error')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

ax_diff = subplot(1,3,3);
imagesc(tspan, zs_Go, grand_avg_iCSD_NC_X - grand_avg_iCSD_Go_X);
c = colorbar;
% c.Label.String = 'nA/mm^{3}';
c.Label.FontSize = font;
[map,~,~,~] = brewermap(256,'PiYG');
colormap(ax_diff, (map));
max_CSD = max(abs(grand_avg_iCSD_NC_X - grand_avg_iCSD_Go_X), [], 'all');
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
% xlabel('Time relative to saccade (ms)', 'FontSize',font);
% ylabel('z (\mum)','FontSize',font);
% title('Error - Correct')
set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')

%% plot all sess CSDs

rows = 6;
cols = 6;

idx = 1:36;
idx = reshape(idx, [6 6])';
idx = reshape(idx(1:2:6,:)',1,[]);
idx = idx(1:end-2);

sessiCSD_Go = struct2cell(structfun(@(str) str.iCSD_Go, iCSD, ...
    'UniformOutput', false));
sessiCSD_NC = struct2cell(structfun(@(str) str.iCSD_NC, iCSD, ...
    'UniformOutput', false));

max_CSD = max(abs([cell2mat(sessiCSD_Go) cell2mat(sessiCSD_NC)]), [], 'all');

%%

figure;
tiledlayout(rows,cols,'TileSpacing','Compact','Padding','Compact');

jj = 1;

for ii = idx
    
    max_CSD = max(abs([sessiCSD_Go{jj, 1} sessiCSD_NC{jj, 1}]), [], 'all');
    
    nexttile(ii);
    imagesc(tspan, zs_Go, sessiCSD_Go{jj, 1}./max_CSD);
    colormap(flip(jet));
    bar_min = -1;
    bar_max = 1;
    caxis([bar_min bar_max]);
    xlim([-500 500])
    xticks([-500 0 500 1000])
    xticklabels({})
    yticks(sort([ze 1125]))
    yticklabels({})
    ax = gca;
    ax.FontSize = font;
    ax.YDir='reverse';
    set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
    
    nexttile(cols + ii);
    imagesc(tspan, zs_Go, sessiCSD_NC{jj, 1}./max_CSD);
    colormap(flip(jet));
    caxis([bar_min bar_max]);
    xlim([-500 500])
    xticks([-500 0 500 1000])
    if ii + 1 < 25
        xticklabels({})
    end
    yticks(sort([ze 1125]))
    if ii==25
        ytick_labels = cellfun(@num2str, num2cell(z_depth), 'UniformOutput', false);
        ytick_labels(~ismember(z_depth, [1125 0 -1125])) = {' '};
        yticklabels(ytick_labels)
    else
        yticklabels({})
    end
    ax = gca;
    ax.FontSize = font;
    ax.YDir='reverse';
    % xlabel('Time from Saccade (s)', 'FontSize',font);
    % ylabel('z (\mum)','FontSize',font);
    % title('Error')
    set(gca,'linewidth',1.5,'fontsize',font,'fontweight','bold')
    
    jj = jj + 1;
    
end

