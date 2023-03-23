%% calculate the spk rate of pyramidal error neurons relative to saccade

clear
clc

noShowFigures = 1; % if =1; figures don't pop out
saveFigures = 0; % if =1; figures are save in jpeg

%% path 2 spiking data

loadPath = fullfile('data', 'monkeys_data','laminar_data'); % IMPORTANT: you
% have to download and save the laminar data in this folder 

%% perpendicular sessions

% sessions number
Sess_EuP1 = 14:19; % Euler, site P1
Sess_XP2P3 = [20:25 26:29]; % Xena, 21-22 site P2, 26-29 site P3
SessNumb = [Sess_EuP1, Sess_XP2P3];

% load electrodes' co-registration across sessions and monkeys
load(fullfile('data','eleAlignment.mat'))
load(fullfile('data','neuronsInfoAmirSteven.mat'))

%% selected neurons based on their spike rate

selL5unitsBasedonSpking = cell(1, length(SessNumb));
selL5unitsBasedonSpking(ismember(SessNumb,[15 17 18 19])) = {[1 2],...
    1, [1 3], [1 2]}; % for sessions with
% error neurons in layer 5. Sessions 15, 17, 18 and 19 in that order.

selL3unitsBasedonSpking = cell(1, length(SessNumb));
selL3unitsBasedonSpking(ismember(SessNumb,14:19)) = {[1 3 4], [1 2], 1,...
    [1 2], 1:3, 1:3};

%% perpendicular sessions

% sessions number
Sess_EuP1 = 14:19; % Euler, site P1
Sess_XP2P3 = [20:25 26:29]; % Xena, 21-22 site P2, 26-29 site P3
SessNumb = [Sess_EuP1, Sess_XP2P3];

% load electrodes' co-registration across sessions and monkeys
load(fullfile('data','eleAlignment.mat'))
load(fullfile('data','neuronsInfoAmirSteven.mat'))

%% electrodes layer assignment

num_layers = 4;
L12 = [1 4]; % mine 1-5
L3 = [5 8]; % mine 6-8
L5 = [9 12];
L6 = [13 16];

%% parameters for the analysis
binsize =  10e-3; % [ms] bin size
gauss_SD = 10; % (10ms) standard deviation for smoothing
flt = true; % if True: use Gaussian filter in smoothdata;
% otherwise: convol with Gaussian

%% create saving directory

savePath = fullfile('data', 'monkeys_data','spkrate_errorneurons');

if ~exist(savePath, 'dir') % checks if the folder already exists
    mkdir(savePath);  % creates a folder named 'file'
end

savePathFig = fullfile(savePath,'Figs');

if ~exist(savePathFig, 'dir') % checks if the folder already exists
    mkdir(savePathFig);  % creates a folder named 'file'
end

%% pre-allocating memory
mean_fr_L5 = cell(length(SessNumb), 2);
mean_fr_L3 = cell(length(SessNumb), 2);
ts_sess = cell(length(SessNumb), 1);

%% analyzing each session individually

for n=1:length(SessNumb)

    sprintf('session %d', SessNumb(n))

    %% --- load data
    load(fullfile(loadPath,['Session_' num2str(SessNumb(n)) '_Spks_Saccade.mat']))

    ierrorCells = find(neuronsInfo.(['Session_'...
        num2str(SessNumb(n))]).ErrorNeurons==1);
    [~,iBScells,~] = intersect(neuronsInfo.(['Session_'...
        num2str(SessNumb(n))]).names,...
        neuronsInfo.(['Session_' num2str(SessNumb(n))]).BSunits); % index
    % of BS cells in original cell names array
    iBSerrorcells = intersect(ierrorCells,iBScells); % BS & error cells
    % index in original cell array of names
    BSerrorcells_depths = neuronsInfo.(['Session_' ...
        num2str(SessNumb(n))]).depth(iBSerrorcells,2); % depth associated with
    % BS error cells
    [BSerrorCells, i_BSerrorcells, ~] = intersect(neuronsInfo.(['Session_'...
        num2str(SessNumb(n))]).names(iBSerrorcells), ...
        spikesGo_Saccade.label); % BS and error cells within
    % spikesGo_Target.label
    depths = BSerrorcells_depths(i_BSerrorcells);

    %% --- select BS error cells before calculating the isi dist
    cfg = [];
    cfg.spikechannel = BSerrorCells;
    spksGo_Target = ft_spike_select(cfg, spikesGo_Saccade);
    spksNC_Target = ft_spike_select(cfg, spikesNC_Saccade);

    %% --- calculate the spike rate per trial for each neuron
    ts_edges = (-0.5):binsize:2; % time window relative to saccade for
    % calculating the spike rate
    ts = ts_edges + mean(diff(ts_edges), 'omitnan');
    ts = ts(1:end-1);
    ts_sess{n} = ts;
    flt = true; % use Gaussian filter in smoothdata or convol with Gaussian
    rate_dataGo = cal_spkrate_ftspkstr(spksGo_Target, ts_edges, gauss_SD + 1, flt);
    rate_dataNC = cal_spkrate_ftspkstr(spksNC_Target, ts_edges, gauss_SD + 1, flt);

    %% --- group neurons per depth and calculate mean firing rate

    firingRateDepth_Go = arrayfun(@(e) cell2mat(rate_dataGo.firingRatet(depths==e)'),...
        1:size(eleAlignment,1), "UniformOutput", false)';
    firingRateDepth_NC = arrayfun(@(e) cell2mat(rate_dataNC.firingRatet(depths==e)'),...
        1:size(eleAlignment,1), "UniformOutput", false)';

    depth_empty_Go = (cellfun(@(x) isempty(x), firingRateDepth_Go)==1);
    depth_empty_NC = (cellfun(@(x) isempty(x), firingRateDepth_NC)==1);

    mfiringRateDepth_Go = cellfun(@(x) mean(x, 1,'omitnan'),...
        firingRateDepth_Go, "UniformOutput", false);
    mfiringRateDepth_Go(depth_empty_Go) = {nan(1,length(ts))};
    mfiringRateDepth_Go = cell2mat(mfiringRateDepth_Go);

    mfiringRateDepth_NC = cellfun(@(x) mean(x, 1,'omitnan'),...
        firingRateDepth_NC, "UniformOutput", false);
    mfiringRateDepth_NC(depth_empty_NC) = {nan(1,length(ts))};
    mfiringRateDepth_NC = cell2mat(mfiringRateDepth_NC);

    %% --- group neurons per depth and calculate mean psth

    psthDepth_Go = arrayfun(@(e) cell2mat(rate_dataGo.pstht(depths==e)'),...
        1:size(eleAlignment,1), "UniformOutput", false)';
    psthDepth_NC = arrayfun(@(e) cell2mat(rate_dataNC.pstht(depths==e)'),...
        1:size(eleAlignment,1), "UniformOutput", false)';

    depth_empty_Go = (cellfun(@(x) isempty(x), psthDepth_Go)==1);
    depth_empty_NC = (cellfun(@(x) isempty(x), psthDepth_NC)==1);

    mpsthDepth_Go = cellfun(@(x) mean(x, 1,'omitnan'),...
        psthDepth_Go, "UniformOutput", false);
    mpsthDepth_Go(depth_empty_Go) = {nan(1,length(ts))};
    mpsthDepth_Go = cell2mat(mpsthDepth_Go);

    mpsthDepth_NC = cellfun(@(x) mean(x, 1,'omitnan'),...
        psthDepth_NC, "UniformOutput", false);
    mpsthDepth_NC(depth_empty_NC) = {nan(1,length(ts))};
    mpsthDepth_NC = cell2mat(mpsthDepth_NC);

    %% --- plot psth and firing rate across depth
    %% -------------- firing rate
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 12 4]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 12 4]);
    end

    cols = 3;
    rows = 1;

    % Ne = 16; % number of electrodes
    x  = [1125:-150:0, -75:-150:-1125]; % mm depth

    minC = min([mfiringRateDepth_Go; mfiringRateDepth_NC],[],'all');
    maxC = max([mfiringRateDepth_Go; mfiringRateDepth_NC],[],'all');
    tmp_diff = mfiringRateDepth_NC - mfiringRateDepth_Go;

    subplot(rows,cols,1)
    imagesc(ts, x, interp2(mfiringRateDepth_Go,5));
    %     contourf(ts, x, tmp_Go,numLines1,'linecolor',colorLines)
    hold on;
    plot(ts, zeros(1,length(ts)),'--k')
    c = colorbar;
    c.Label.FontSize=12;
    c.Label.String = 'Spk/s';
    if minC ~= maxC && ~isnan(minC) && ~isnan(maxC)
        caxis([minC maxC]);
    end
    %     caxis([minC maxC])
    ylabel('Depth (\mum)')
    xlabel('Time from saccade')
    title('Go trials')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    subplot(rows,cols,2)
    imagesc(ts, x, interp2(mfiringRateDepth_NC,5));
    %     contourf(ts, x, tmp_NC,numLines1,'linecolor',colorLines)
    hold on;
    plot((ts), zeros(1,length(ts)),'--k')
    c = colorbar;
    c.Label.String = 'Spk/s';
    if minC ~= maxC && ~isnan(minC) && ~isnan(maxC)
        caxis([minC maxC]);
    end
    %     caxis([minC maxC])
    c.Label.FontSize=12;
    xlabel('Time from saccade')
    title('NC trials')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    subplot(rows,cols,3)
    imagesc((ts), x, interp2(tmp_diff,5));
    %     contourf(ts, x, tmp_diff,numLines2,'linecolor',colorLines)
    hold on;
    plot((ts), zeros(1,length(ts)),'--k')
    c = colorbar;
    c.Label.String = '\Delta(NC-Go) (spks/s)';
    %     c.Label.FontSize=12;
    %     caxis([-1 1].*valueM)
    hold on;
    xlabel('Time from saccade')
    title('Difference')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    if noShowFigures == 1
        % Set CreateFcn callback
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end

    % Save Fig file
    if saveFigures==1
        saveas(hFig, fullfile(savePathFig,...
            ['frDepth_Session' num2str(SessNumb(n)) '.jpg']));
    end

    if noShowFigures == 1; close all; end

    %% -------------- psth
    if noShowFigures == 1
        hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 12 4]);
    else
        hFig = figure('Units', 'inches','Position',[0.05 0.05 12 4]);
    end

    minC = min([mpsthDepth_Go; mpsthDepth_NC],[],'all');
    maxC = min([mpsthDepth_Go; mpsthDepth_NC],[],'all');
    tmp_diff = mpsthDepth_NC - mpsthDepth_Go;

    subplot(rows,cols,1)
    imagesc(ts, x, interp2(mpsthDepth_Go,5));
    %     contourf(ts, x, tmp_Go,numLines1,'linecolor',colorLines)
    hold on;
    plot(ts, zeros(1,length(ts)),'--k');
    c = colorbar;
    c.Label.FontSize=12;
    c.Label.String = 'Spk/s';
    if minC ~= maxC && ~isnan(minC) && ~isnan(maxC)
        caxis([minC maxC]);
    end
    ylabel('Depth (\mum)')
    xlabel('Time from saccade')
    title('Go trials')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    subplot(rows,cols,2)
    imagesc(ts, x, interp2(mpsthDepth_NC,5));
    %     contourf(ts, x, tmp_NC,numLines1,'linecolor',colorLines)
    hold on;
    plot((ts), zeros(1,length(ts)),'--k')
    c = colorbar;
    c.Label.String = 'Spk/s';
    if minC ~= maxC && ~isnan(minC) && ~isnan(maxC)
        caxis([minC maxC]);
    end
    %     caxis([minC maxC]);
    c.Label.FontSize=12;
    xlabel('Time from saccade')
    title('NC trials')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    subplot(rows,cols,3)
    imagesc((ts), x, interp2(tmp_diff,5));
    %     contourf(ts, x, tmp_diff,numLines2,'linecolor',colorLines)
    hold on;
    plot((ts), zeros(1,length(ts)),'--k')
    c = colorbar;
    c.Label.String = '\Delta(NC-Go) (Spk/s)';
    %     c.Label.FontSize=12;
    %     caxis([-1 1].*valueM)
    hold on;
    xlabel('Time from saccade')
    title('Difference')
    set(gca,'LineWidth',2, 'FontSize', 12,'fontweight','bold','YDir','normal')

    if noShowFigures == 1
        % Set CreateFcn callback
        set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
    end

    % Save Fig file
    if saveFigures==1
        saveas(hFig, fullfile(savePathFig,...
            ['psthDepth_Session' num2str(SessNumb(n)) '.jpg']));
    end

    if noShowFigures == 1; close all; end

    %% --- get firing rate of neurons in layer 5
    % get mean firing rate across trials of all neurons in layer 5 and
    % compute the mean across neurons
    firingRateL5_Go = cell2mat(rate_dataGo.firingRatet(:,depths>=L5(1)&depths<=L5(2))');
    firingRateL5_NC = cell2mat(rate_dataNC.firingRatet(:,depths>=L5(1)&depths<=L5(2))');

    mfiringRateL5_Go = mean(firingRateL5_Go, 1,'omitnan');
    mfiringRateL5_NC = mean(firingRateL5_NC, 1,'omitnan');

    mean_fr_L5{n,1} = firingRateL5_Go(selL5unitsBasedonSpking{n},:);
    mean_fr_L5{n,2} = firingRateL5_NC(selL5unitsBasedonSpking{n},:);

    % get the firing rate per trial of all layer 5 neurons
    firingRateL5_trials_Go = rate_dataGo.firingRatei(:,depths>=L5(1)&depths<=L5(2));
    firingRateL5_trials_NC = rate_dataNC.firingRatei(:,depths>=L5(1)&depths<=L5(2));
    firingRateL5_trials3D_Go = cell2mat(vertcat(firingRateL5_trials_Go'));
    firingRateL5_trials3D_NC = cell2mat(vertcat(firingRateL5_trials_NC'));

    % calculate the SEM across trails for each neuron
    % -- go trials
    firingRateL5_sem_Go = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL5_trials_Go, 'UniformOutput', false)';
    firingRateL5_sem_Go = cell2mat(firingRateL5_sem_Go);
%     firingRateL5_all_sem_Go = std(firingRateL5_trials3D_Go,'omitnan')...
%         ./sqrt(size(firingRateL5_trials3D_Go,1));
    % -- NC trials
    firingRateL5_sem_NC = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL5_trials_NC, 'UniformOutput', false)';
    firingRateL5_sem_NC = cell2mat(firingRateL5_sem_NC);
%     firingRateL5_all_sem_NC = std(firingRateL5_trials3D_NC,'omitnan')...
%         ./sqrt(size(firingRateL5_trials3D_NC,1));

    % calculate the SEM of across the mean firing rate of the neurons
    firingRateL5_all_sem_Go = std(firingRateL5_Go, 0, 1,...
        'omitnan')./sqrt(size(firingRateL5_Go,1));
    firingRateL5_all_sem_NC = std(firingRateL5_NC, 0, 1,...
        'omitnan')./sqrt(size(firingRateL5_NC,1));
    
    [peakfiringRateL5_Go,ipeakfiringRateL5_Go] = max(firingRateL5_Go, [], 2);
    [peakfiringRateL5_NC,ipeakfiringRateL5_NC] = max(firingRateL5_NC, [], 2);

    %% ----------- figure: mean rate L5 neurons across trials
    if ~isempty(firingRateL5_Go) || ~isempty(firingRateL5_NC)

        if noShowFigures == 1
            hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 8]);
        else
            hFig = figure('Units', 'inches','Position',[0.05 0.05 10 8]);
        end

        % filling the area between these curve
        %     colors = {'r', 'g', 'b', 'c', 'm', 'k',...
        %         [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
        %         [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330],...
        %         [0.6350 0.0780 0.1840]};

        cols = 3;
        rows = 2; %round((size(firingRateL5_Go,1) + 1)/cols);


        %     hold on;
        for ii=1:size(firingRateL5_Go,1)
            subplot(rows,cols,ii)
            maxVn = max([firingRateL5_Go(ii,:) +(1.96.*firingRateL5_sem_Go(ii,:));...
                firingRateL5_NC(ii,:) +(1.96.*firingRateL5_sem_NC(ii,:))],[],'all')+15;

            % -- Go trace
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL5_Go(ii,:),...
                fliplr(firingRateL5_Go(ii,:)-(1.96.*firingRateL5_sem_Go(ii,:)))];
            fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL5_Go(ii,:)+(1.96.*firingRateL5_sem_Go(ii,:)),...
                fliplr(firingRateL5_Go(ii,:))];
            fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p1 = plot(ts, firingRateL5_Go(ii,:), 'Color', 'r', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL5_Go(ii)), peakfiringRateL5_Go(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 8,...
                num2str(ts(ipeakfiringRateL5_Go(ii))'),...
                'HorizontalAlignment','center', 'Color', 'r') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            % -- NC trace
            hold on;
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL5_NC(ii,:),...
                fliplr(firingRateL5_NC(ii,:)-(1.96.*firingRateL5_sem_NC(ii,:)))];
            fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL5_NC(ii,:)+(1.96.*firingRateL5_sem_NC(ii,:)),...
                fliplr(firingRateL5_NC(ii,:))];
            fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p2 = plot(ts, firingRateL5_NC(ii,:), 'Color', 'b', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL5_NC(ii)), peakfiringRateL5_NC(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 14,...
                num2str(ts(ipeakfiringRateL5_NC(ii))'),...
                'HorizontalAlignment','center', 'Color', 'b') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            title(sprintf('Neuron %d', ii))
            legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
            ylim([0 maxVn])
            xlim([min(ts), max(ts)])
            xlabel('Time from saccade')
            set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
        end

        % -- mean across neurons
        maxV = max([mfiringRateL5_Go +(1.96.*firingRateL5_all_sem_Go);...
            mfiringRateL5_NC +(1.96.*firingRateL5_all_sem_NC)],[],'all')+15;

        subplot(rows,cols,size(firingRateL5_Go,1)+1)
        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL5_Go,...
            fliplr(mfiringRateL5_Go-(1.96.*firingRateL5_all_sem_Go))];
        fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL5_Go+(1.96.*firingRateL5_all_sem_Go),...
            fliplr(mfiringRateL5_Go)];
        fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p1 = plot(ts, mfiringRateL5_Go, 'Color', 'r', 'linewidth',1.5);
        hold on;
        [M_Go, i_Go] = max(mfiringRateL5_Go);
        plot(ts(i_Go).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_Go), M_Go+3, num2str(ts(i_Go)),...
            'Color','k','FontWeight','bold')

        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL5_NC,...
            fliplr(mfiringRateL5_NC-(1.96.*firingRateL5_all_sem_NC))];
        fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL5_NC+(1.96.*firingRateL5_all_sem_NC),...
            fliplr(mfiringRateL5_NC)];
        fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p2 = plot(ts, mfiringRateL5_NC, 'Color', 'b', 'linewidth',1.5);
        hold on;
        [M_NC, i_NC] = max(mfiringRateL5_NC);
        plot(ts(i_NC).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_NC), M_NC+3, num2str(ts(i_NC)),...
            'Color','k','FontWeight','bold')

        ylim([0 maxV])
        xlim([min(ts), max(ts)])
        legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
        title('Mean across neurons')
        xlabel('Time from saccade')
        set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')

        if noShowFigures == 1
            % Set CreateFcn callback
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
        end

        % Save Fig file
        if saveFigures==1; saveas(hFig, fullfile(savePathFig, ...
                ['frL5units_Session_' num2str(SessNumb(n)) '.jpg'])); end

        if noShowFigures == 1; close all; end

    end
    %% --- get firing rate of neurons in layer 3
    % get mean firing rate across trials of all neurons in layer 3 and
    % compute the mean across neurons
    firingRateL3_Go = cell2mat(rate_dataGo.firingRatet(:,depths>=L3(1)&depths<=L3(2))');
    firingRateL3_NC = cell2mat(rate_dataNC.firingRatet(:,depths>=L3(1)&depths<=L3(2))');
    mfiringRateL3_Go = mean(firingRateL3_Go, 1,'omitnan');
    mfiringRateL3_NC = mean(firingRateL3_NC, 1,'omitnan');

    mean_fr_L3{n,1} = firingRateL3_Go(selL3unitsBasedonSpking{n},:);
    mean_fr_L3{n,2} = firingRateL3_NC(selL3unitsBasedonSpking{n},:);

    % get the firing rate per trial of all layer 3 neurons
    firingRateL3_trials_Go = rate_dataGo.firingRatei(:,depths>=L3(1)&depths<=L3(2));
    firingRateL3_trials_NC = rate_dataNC.firingRatei(:,depths>=L3(1)&depths<=L3(2));
    firingRateL3_trials3D_Go = cell2mat(vertcat(firingRateL3_trials_Go'));
    firingRateL3_trials3D_NC = cell2mat(vertcat(firingRateL3_trials_NC'));

    % calculate the SEM across trails for each neuron
    % -- go trials
    firingRateL3_sem_Go = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL3_trials_Go, 'UniformOutput', false)';
    firingRateL3_sem_Go = cell2mat(firingRateL3_sem_Go);
%     firingRateL3_all_sem_Go = std(firingRateL3_trials3D_Go,'omitnan')...
%         ./sqrt(size(firingRateL3_trials3D_Go,1));
    % -- NC trials
    firingRateL3_sem_NC = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL3_trials_NC, 'UniformOutput', false)';
    firingRateL3_sem_NC = cell2mat(firingRateL3_sem_NC);
%     firingRateL3_all_sem_NC = std(firingRateL3_trials3D_NC,'omitnan')...
%         ./sqrt(size(firingRateL3_trials3D_NC,1));

    % calculate the SEM of across the mean firing rate of the neurons
    firingRateL3_all_sem_Go = std(firingRateL3_Go, 0, 1,...
        'omitnan')./sqrt(size(firingRateL3_Go,1));
    firingRateL3_all_sem_NC = std(firingRateL3_NC, 0, 1,...
        'omitnan')./sqrt(size(firingRateL3_NC,1));

    [peakfiringRateL3_Go,ipeakfiringRateL3_Go] = max(firingRateL3_Go, [], 2);
    [peakfiringRateL3_NC,ipeakfiringRateL3_NC] = max(firingRateL3_NC, [], 2);

    %% ----------- figure: mean rate L3 neurons across trials
    if ~isempty(firingRateL3_Go) || ~isempty(firingRateL3_NC)

        if noShowFigures == 1
            hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 8]);
        else
            hFig = figure('Units', 'inches','Position',[0.05 0.05 10 8]);
        end

        % filling the area between these curve
        %     colors = {'r', 'g', 'b', 'c', 'm', 'k',...
        %         [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
        %         [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330],...
        %         [0.6350 0.0780 0.1840]};

        cols = 3;
        rows = 2; %round((size(firingRateL3_Go,1) + 1)/cols);


        %     hold on;
        for ii=1:size(firingRateL3_Go,1)
            subplot(rows,cols,ii)
            maxVn = max([firingRateL3_Go(ii,:) +(1.96.*firingRateL3_sem_Go(ii,:));...
                firingRateL3_NC(ii,:) +(1.96.*firingRateL3_sem_NC(ii,:))],[],'all')+15;

            % -- Go trace
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL3_Go(ii,:),...
                fliplr(firingRateL3_Go(ii,:)-(1.96.*firingRateL3_sem_Go(ii,:)))];
            fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL3_Go(ii,:)+(1.96.*firingRateL3_sem_Go(ii,:)),...
                fliplr(firingRateL3_Go(ii,:))];
            fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p1 = plot(ts, firingRateL3_Go(ii,:), 'Color', 'r', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL3_Go(ii)), peakfiringRateL3_Go(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 8,...
                num2str(ts(ipeakfiringRateL3_Go(ii))'),...
                'HorizontalAlignment','center', 'Color', 'r') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            % -- NC trace
            hold on;
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL3_NC(ii,:),...
                fliplr(firingRateL3_NC(ii,:)-(1.96.*firingRateL3_sem_NC(ii,:)))];
            fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL3_NC(ii,:)+(1.96.*firingRateL3_sem_NC(ii,:)),...
                fliplr(firingRateL3_NC(ii,:))];
            fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p2 = plot(ts, firingRateL3_NC(ii,:), 'Color', 'b', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL3_NC(ii)), peakfiringRateL3_NC(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 14,...
                num2str(ts(ipeakfiringRateL3_NC(ii))'),...
                'HorizontalAlignment','center', 'Color', 'b') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            title(sprintf('Neuron %d', ii))
            legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
            ylim([0 maxVn])
            xlim([min(ts), max(ts)])
            xlabel('Time from saccade')
            set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
        end

        % -- mean across neurons
        maxV = max([mfiringRateL3_Go +(1.96.*firingRateL3_all_sem_Go);...
            mfiringRateL3_NC +(1.96.*firingRateL3_all_sem_NC)],[],'all')+5;

        subplot(rows,cols,size(firingRateL3_Go,1)+1)
        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL3_Go,...
            fliplr(mfiringRateL3_Go-(1.96.*firingRateL3_all_sem_Go))];
        fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL3_Go+(1.96.*firingRateL3_all_sem_Go),...
            fliplr(mfiringRateL3_Go)];
        fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p1 = plot(ts, mfiringRateL3_Go, 'Color', 'r', 'linewidth',1.5);
        hold on;
        [M_Go, i_Go] = max(mfiringRateL3_Go);
        plot(ts(i_Go).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_Go), M_Go+3, num2str(ts(i_Go)),...
            'Color','k','FontWeight','bold')

        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL3_NC,...
            fliplr(mfiringRateL3_NC-(1.96.*firingRateL3_all_sem_NC))];
        fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL3_NC+(1.96.*firingRateL3_all_sem_NC),...
            fliplr(mfiringRateL3_NC)];
        fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p2 = plot(ts, mfiringRateL3_NC, 'Color', 'b', 'linewidth',1.5);
        hold on;
        [M_NC, i_NC] = max(mfiringRateL3_NC);
        plot(ts(i_NC).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_NC), M_NC+3, num2str(ts(i_NC)),...
            'Color','k','FontWeight','bold')

        ylim([0 maxV])
        xlim([min(ts), max(ts)])
        legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
        title('Mean across neurons')
        xlabel('Time from saccade')
        set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')

        if noShowFigures == 1
            % Set CreateFcn callback
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
        end

        % Save Fig file
        if saveFigures==1; saveas(hFig, fullfile(savePathFig, ...
                ['frL3units_Session_' num2str(SessNumb(n)) '.jpg'])); end

        if noShowFigures == 1; close all; end
    end
    %% --- get firing rate of neurons in layer 6
    % get mean firing rate across trials of all neurons in layer 6 and
    % compute the mean across neurons
    firingRateL6_Go = cell2mat(rate_dataGo.firingRatet(:,depths>=L6(1)&depths<=L6(2))');
    firingRateL6_NC = cell2mat(rate_dataNC.firingRatet(:,depths>=L6(1)&depths<=L6(2))');

    mfiringRateL6_Go = mean(firingRateL6_Go, 1,'omitnan');
    mfiringRateL6_NC = mean(firingRateL6_NC, 1,'omitnan');

    % get the firing rate per trial of all layer 6 neurons
    firingRateL6_trials_Go = rate_dataGo.firingRatei(:,depths>=L6(1)&depths<=L6(2));
    firingRateL6_trials_NC = rate_dataNC.firingRatei(:,depths>=L6(1)&depths<=L6(2));
    firingRateL6_trials3D_Go = cell2mat(vertcat(firingRateL6_trials_Go'));
    firingRateL6_trials3D_NC = cell2mat(vertcat(firingRateL6_trials_NC'));

    % calculate the SEM across trails for each neuron
    % -- go trials
    firingRateL6_sem_Go = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL6_trials_Go, 'UniformOutput', false)';
    firingRateL6_sem_Go = cell2mat(firingRateL6_sem_Go);
%     firingRateL6_all_sem_Go = std(firingRateL6_trials3D_Go,'omitnan')...
%         ./sqrt(size(firingRateL6_trials3D_Go,1));
    % -- NC trials
    firingRateL6_sem_NC = cellfun(@(x) std(x,'omitnan')./sqrt(size(x,1)),...
        firingRateL6_trials_NC, 'UniformOutput', false)';
    firingRateL6_sem_NC = cell2mat(firingRateL6_sem_NC);
%     firingRateL6_all_sem_NC = std(firingRateL6_trials3D_NC,'omitnan')...
%         ./sqrt(size(firingRateL6_trials3D_NC,1));

    % calculate the SEM of across the mean firing rate of the neurons
    firingRateL6_all_sem_Go = std(firingRateL6_Go, 0, 1,...
        'omitnan')./sqrt(size(firingRateL6_Go,1));
    firingRateL6_all_sem_NC = std(firingRateL6_NC, 0, 1,...
        'omitnan')./sqrt(size(firingRateL6_NC,1));

    [peakfiringRateL6_Go,ipeakfiringRateL6_Go] = max(firingRateL6_Go, [], 2);
    [peakfiringRateL6_NC,ipeakfiringRateL6_NC] = max(firingRateL6_NC, [], 2);

    %% ----------- figure: mean rate L6 neurons across trials
    if ~isempty(firingRateL6_Go) || ~isempty(firingRateL6_NC)

        if noShowFigures == 1
            hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 8]);
        else
            hFig = figure('Units', 'inches','Position',[0.05 0.05 10 8]);
        end

        % filling the area between these curve
        %     colors = {'r', 'g', 'b', 'c', 'm', 'k',...
        %         [0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250],...
        %         [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330],...
        %         [0.6350 0.0780 0.1840]};

        cols = 3;
        rows = 2; %round((size(firingRateL6_Go,1) + 1)/cols);


        %     hold on;
        for ii=1:size(firingRateL6_Go,1)
            subplot(rows,cols,ii)
            maxVn = max([firingRateL6_Go(ii,:) +(1.96.*firingRateL6_sem_Go(ii,:));...
                firingRateL6_NC(ii,:) +(1.96.*firingRateL6_sem_NC(ii,:))],[],'all')+15;

            % -- Go trace
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL6_Go(ii,:),...
                fliplr(firingRateL6_Go(ii,:)-(1.96.*firingRateL6_sem_Go(ii,:)))];
            fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL6_Go(ii,:)+(1.96.*firingRateL6_sem_Go(ii,:)),...
                fliplr(firingRateL6_Go(ii,:))];
            fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p1 = plot(ts, firingRateL6_Go(ii,:), 'Color', 'r', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL6_Go(ii)), peakfiringRateL6_Go(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 8,...
                num2str(ts(ipeakfiringRateL6_Go(ii))'),...
                'HorizontalAlignment','center', 'Color', 'r') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            % -- NC trace
            hold on;
            x2 = [ts, fliplr(ts)];
            inBetween1 = [firingRateL6_NC(ii,:),...
                fliplr(firingRateL6_NC(ii,:)-(1.96.*firingRateL6_sem_NC(ii,:)))];
            fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            inBetween2 = [firingRateL6_NC(ii,:)+(1.96.*firingRateL6_sem_NC(ii,:)),...
                fliplr(firingRateL6_NC(ii,:))];
            fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            hold on;
            p2 = plot(ts, firingRateL6_NC(ii,:), 'Color', 'b', 'linewidth',1.5);
            hold on;
            plot(ts(ipeakfiringRateL6_NC(ii)), peakfiringRateL6_NC(ii), 'vb',...
                'MarkerFaceColor','b');
            text(0.5, maxVn - 14,...
                num2str(ts(ipeakfiringRateL6_NC(ii))'),...
                'HorizontalAlignment','center', 'Color', 'b') %,...
            %             'BackgroundColor',[0.5, 0.5, 0.5])

            title(sprintf('Neuron %d', ii))
            legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
            ylim([0 maxVn])
            xlim([min(ts), max(ts)])
            xlabel('Time from saccade')
            set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
        end

        % -- mean across neurons
        maxV = max([mfiringRateL6_Go +(1.96.*firingRateL6_all_sem_Go);...
            mfiringRateL6_NC +(1.96.*firingRateL6_all_sem_NC)],[],'all')+5;

        subplot(rows,cols,size(firingRateL6_Go,1)+1)
        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL6_Go,...
            fliplr(mfiringRateL6_Go-(1.96.*firingRateL6_all_sem_Go))];
        fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL6_Go+(1.96.*firingRateL6_all_sem_Go),...
            fliplr(mfiringRateL6_Go)];
        fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p1 = plot(ts, mfiringRateL6_Go, 'Color', 'r', 'linewidth',1.5);
        hold on;
        [M_Go, i_Go] = max(mfiringRateL6_Go);
        plot(ts(i_Go).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_Go), M_Go+3, num2str(ts(i_Go)),...
            'Color','k','FontWeight','bold')

        hold on;
        x2 = [ts, fliplr(ts)];
        inBetween1 = [mfiringRateL6_NC,...
            fliplr(mfiringRateL6_NC-(1.96.*firingRateL6_all_sem_NC))];
        fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        inBetween2 = [mfiringRateL6_NC+(1.96.*firingRateL6_all_sem_NC),...
            fliplr(mfiringRateL6_NC)];
        fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on;
        p2 = plot(ts, mfiringRateL6_NC, 'Color', 'b', 'linewidth',1.5);
        hold on;
        [M_NC, i_NC] = max(mfiringRateL6_NC);
        plot(ts(i_NC).*ones(1,length(0:maxV)), 0:maxV, '-k')
        text(ts(i_NC), M_NC+3, num2str(ts(i_NC)),...
            'Color','k','FontWeight','bold')

        ylim([0 maxV])
        xlim([min(ts), max(ts)])
        legend([p1, p2], {'Go trial', 'NC trial'},'Location','northwest')
        title('Mean across neurons')
        xlabel('Time from saccade')
        set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')

        if noShowFigures == 1
            % Set CreateFcn callback
            set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
        end

        % Save Fig file
        if saveFigures==1; saveas(hFig, fullfile(savePathFig, ...
                ['frL6units_Session_' num2str(SessNumb(n)) '.jpg'])); end

        if noShowFigures == 1; close all; end
    end
end

%% save 

mean_fr_L5 = mean_fr_L5(cellfun(@(x) size(x,1)>0, mean_fr_L5(:,1)), :);
mean_fr_L3 = mean_fr_L3(cellfun(@(x) size(x,1)>0, mean_fr_L3(:,1)), :);

mean_fr_L3_Go = cell2mat(mean_fr_L3(:,1));
mean_fr_L5_Go = cell2mat(mean_fr_L5(:,1));
mean_fr_L3_NC = cell2mat(mean_fr_L3(:,2));
mean_fr_L5_NC = cell2mat(mean_fr_L5(:,2));

L3_mean_fr_Go = mean(mean_fr_L3_Go, 1, 'omitnan');
L5_mean_fr_Go = mean(mean_fr_L5_Go, 1, 'omitnan'); 
L3_mean_fr_NC = mean(mean_fr_L3_NC, 1, 'omitnan');
L5_mean_fr_NC = mean(mean_fr_L5_NC, 1, 'omitnan'); 

L3_sem_fr_Go = std(mean_fr_L3_Go, 0, 1, ...
    'omitnan')./sqrt(size(mean_fr_L3_Go,1));
L5_sem_fr_Go = std(mean_fr_L5_Go, 0, 1, ...
    'omitnan')./sqrt(size(mean_fr_L5_Go,1)); 

L3_sem_fr_NC = std(mean_fr_L3_NC, 0, 1, ...
    'omitnan')./sqrt(size(mean_fr_L3_NC,1));
L5_sem_fr_NC = std(mean_fr_L5_NC, 0, 1, ...
    'omitnan')./sqrt(size(mean_fr_L5_NC,1)); 

save(fullfile(savePath, 'fr_error_cells_saccade.mat'), 'mean_fr_L5_Go', ...
    'mean_fr_L5_NC', 'mean_fr_L3_Go', 'mean_fr_L3_NC', 'L3_mean_fr_Go', ...
    'L3_mean_fr_NC', 'L5_mean_fr_Go', 'L5_mean_fr_NC', 'L3_sem_fr_Go', ...
    'L3_sem_fr_NC', 'L5_sem_fr_Go', 'L5_sem_fr_NC', 'ts_sess')

%% L3
ts_plot = ts_sess{1};

maxV = max([L3_mean_fr_Go +(1.96.*L3_sem_fr_Go);...
            L3_mean_fr_NC +(1.96.*L3_sem_fr_NC)],[],'all')+5;

figure;
hold on;
x2 = [ts_plot, fliplr(ts_plot)];
inBetween1 = [L3_mean_fr_Go,...
    fliplr(L3_mean_fr_Go-(1.96.*L3_sem_fr_Go))];
fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
inBetween2 = [L3_mean_fr_Go+(1.96.*L3_sem_fr_Go),...
    fliplr(L3_mean_fr_Go)];
fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(ts_plot, L3_mean_fr_Go, 'Color', 'r', 'linewidth',1.5);
hold on;
[M_Go, i_Go] = max(L3_mean_fr_Go);
plot(ts_plot(i_Go).*ones(1,length(0:maxV)), 0:maxV, '-k')
text(ts_plot(i_Go), M_Go+3, num2str(ts_plot(i_Go)),...
    'Color','k','FontWeight','bold')
hold on;
x2 = [ts_plot, fliplr(ts_plot)];
inBetween1 = [L3_mean_fr_NC,...
    fliplr(L3_mean_fr_NC-(1.96.*L3_sem_fr_NC))];
fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
inBetween2 = [L3_mean_fr_NC+(1.96.*L3_sem_fr_NC),...
    fliplr(L3_mean_fr_NC)];
fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
p1 = plot(ts_plot, L3_mean_fr_NC, 'Color', 'b', 'linewidth',1.5);
hold on;
[M_NC, i_NC] = max(L3_mean_fr_NC);
plot(ts_plot(i_NC).*ones(1,length(0:maxV)), 0:maxV, '-k')
text(ts_plot(i_NC), M_NC+3, num2str(ts_plot(i_NC)),...
    'Color','k','FontWeight','bold')
xlim([ts_plot(1) ts_plot(end)])
xlabel('Time from saccade onset')
ylabel('Spks/s  ')
set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')

%% L5

maxV = max([L5_mean_fr_Go +(1.96.*L5_sem_fr_Go);...
            L5_mean_fr_NC +(1.96.*L5_sem_fr_NC)],[],'all')+5;

figure;
hold on;
x2 = [ts_plot, fliplr(ts_plot)];
inBetween1 = [L5_mean_fr_Go,...
    fliplr(L5_mean_fr_Go-(1.96.*L5_sem_fr_Go))];
fill(x2, inBetween1, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
inBetween2 = [L5_mean_fr_Go+(1.96.*L5_sem_fr_Go),...
    fliplr(L5_mean_fr_Go)];
fill(x2, inBetween2, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(ts_plot, L5_mean_fr_Go, 'Color', 'r', 'linewidth',1.5);
hold on;
[M_Go, i_Go] = max(L5_mean_fr_Go);
plot(ts_plot(i_Go).*ones(1,length(0:maxV)), 0:maxV, '-k')
text(ts_plot(i_Go), M_Go+3, num2str(ts_plot(i_Go)),...
    'Color','k','FontWeight','bold')
hold on;
x2 = [ts_plot, fliplr(ts_plot)];
inBetween1 = [L5_mean_fr_NC,...
    fliplr(L5_mean_fr_NC-(1.96.*L5_sem_fr_NC))];
fill(x2, inBetween1, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
inBetween2 = [L5_mean_fr_NC+(1.96.*L5_sem_fr_NC),...
    fliplr(L5_mean_fr_NC)];
fill(x2, inBetween2, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(ts_plot, L5_mean_fr_NC, 'Color', 'b', 'linewidth',1.5);
hold on;
[M_NC, i_NC] = max(L5_mean_fr_NC);
plot(ts_plot(i_NC).*ones(1,length(0:maxV)), 0:maxV, '-k')
text(ts_plot(i_NC), M_NC+3, num2str(ts_plot(i_NC)),...
    'Color','k','FontWeight','bold')
xlim([ts_plot(1) ts_plot(end)])
xlabel('Time from saccade onset')
ylabel('Spks/s  ')
set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
