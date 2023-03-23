%% calculate the ISI distribution of error neurons in the pre-target window

clear
clc

noShowFigures = 1; % if =1; figures don't pop out
saveFigures = 1; % if =1; save figures (jpeg)

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

%% electrodes layer assignment

num_layers = 4;
L12 = [1 4]; % mine 1-5
L3 = [5 8]; % mine 6-8
L5 = [9 12];
L6 = [13 16];

%% ISI dist calculation parameters

max_isi = 0.1; % sec, max value for isi
bin_size = 0.002; % sec, bin size
latensy_str = 'prestim'; % latency

%% create saving directory

savePath = fullfile('data', 'monkeys_data','isi_pretarget_errorBS');

if ~exist(savePath, 'dir') % checks if the folder already exists
    mkdir(savePath);  % creates a folder named 'file'
end

savePathSubFigs = fullfile(savePath,'Figs');

if ~exist(savePathSubFigs, 'dir') % checks if the folder already exists
    mkdir(savePathSubFigs);  % creates a folder named 'file'
end

%% selected neurons based on their spike rate

selL5unitsBasedonSpking = cell(1, length(SessNumb));
selL5unitsBasedonSpking(ismember(SessNumb,[15 17 18 19])) = {[1 2],...
    1, [1 3], [1 2]}; % for sessions with error neurons in layer 5. 
% Sessions 15, 17, 18 and 19 in that order.

selL3unitsBasedonSpking = cell(1, length(SessNumb));
selL3unitsBasedonSpking(ismember(SessNumb,14:19)) = {[1 3 4], [1 2], 1,...
    [1 2], 1:3, 1:3}; % for sessions with
% error neurons in layer 3.

% selL6unitsBasedonSpking = cell(1, length(SessNumb));
% selL6unitsBasedonSpking(ismember(SessNumb,[14 15])) = {1:3, 1}; % for 
% sessions with error neurons in layer 6.

%% pre-allocating memory

minISI_allcells = cell(1,length(SessNumb)); % minimum isi;
% {session_num}[num_units x trial_type]
L3cells_names = cell(1,length(SessNumb)); % name of the cells in layer 3, 
% all sessions
L5cells_names = cell(1,length(SessNumb)); % name of the cells in layer 5, 
% all sessions

isih_L3 = []; % isi dist of L3 BS cells
isidens_L3 = []; % isi_n vs isi_n+1 dist map
isih_L5 = []; % isi dist of L5 BS cells
isidens_L5 = []; % isi_n vs isi_n+1 dist map

mean_preTarget_wind = zeros(length(SessNumb),2); % mean pre target window length

%% analyzing each session individually

for n=1:length(SessNumb)

    sprintf('session %d', SessNumb(n))

    %% --- load data
    load(fullfile(loadPath,['Session_' num2str(SessNumb(n)) '_Spks_Target.mat']))

    %% select BS error neurons
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
        spikesGo_Target.label); % BS and error cells within
    % spikesGo_Target.label
    depths = BSerrorcells_depths(i_BSerrorcells);    

    %% --- select BS error cells before calculating the isi dist
    cfg = [];
    cfg.spikechannel = BSerrorCells;
    spksGo_Target = ft_spike_select(cfg, spikesGo_Target);
    spksNC_Target = ft_spike_select(cfg, spikesNC_Target);

    %% --- get pre target window 

    mean_preTarget_wind(n,:) = [mean(abs(spksGo_Target.trialtime(:,1))) ...
        mean(abs(spksNC_Target.trialtime(:,1)))];

    %% --- calculate the isi for each neuron

    clearvars isih_Go isih_NC

    cfg = [];
    cfg.bins = 0:bin_size:max_isi;
    cfg.latency = latensy_str;
    isih_Go = ft_spike_isi(cfg, spksGo_Target);
    isih_NC = ft_spike_isi(cfg, spksNC_Target);

    myFind = @(time_stamp, x) time_stamp(find(x,1,'first'));
    minisi_Go = cellfun(myFind, repmat({isih_Go.time}, ...
        size(isih_Go.avg,1), 1),num2cell(isih_Go.avg,2) ...
        ,'UniformOutput',false);
    minisi_NC = cellfun(myFind, repmat({isih_NC.time}, ...
        size(isih_NC.avg,1), 1),num2cell(isih_NC.avg,2) ...
        ,'UniformOutput',false);

    % ceheck if there are cells with no spikes
    i_emptyGo = cell2mat(cellfun(@isempty, minisi_Go, 'UniformOutput',false));
    i_emptyNC = cell2mat(cellfun(@isempty, minisi_NC, 'UniformOutput',false));
    if nnz(i_emptyGo) > 0
        minisi_Go(i_emptyGo) = {nan};
    end
    if nnz(i_emptyNC) > 0
        minisi_NC(i_emptyNC) = {nan};
    end
    minISI_allcells{n} = [cell2mat(minisi_Go), cell2mat(minisi_NC)];

    %% -------- save isi for the session
    save(fullfile(savePathSubFolder,['isi_Sess' num2str(SessNumb(n)) '.mat']), ...
        'isih_Go','isih_NC', 'depths')

    %% --- devide cells per layer

    L5cells = BSerrorCells((depths>=L5(1)&depths<=L5(2)));
    L3cells = BSerrorCells((depths>=L3(1)&depths<=L3(2)));

    L3cells_names{n} = L3cells;
    L5cells_names{n} = L5cells;

    %% plots
    %% L5
    if ~isempty(L5cells(selL5unitsBasedonSpking{n}))
        clearvars hdl_Go_L5 hdl_NC_L5 dens_Go_L5 dens_NC_L5

        for ii = selL5unitsBasedonSpking{n}
            
            if noShowFigures == 1
                hFig = figure('Visible', 'off','Units', 'inches', ...
                    'Position',[0.05 0.05 10 5]);
            else
                hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
            end

            cfg = [];
            cfg.spikechannel = L5cells{ii};
            cfg.interpolate = 5;
            cfg.window = 'gausswin';
            cfg.winlen = 0.005;
            cfg.scatter = 'no';
            cfg.colormap = jet(300);

            subplot(1,2,1)
            [~, hdl_Go_L5.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_Go);

            subplot(1,2,2)
            [~, hdl_NC_L5.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_NC);

            if noShowFigures == 1
                % Set CreateFcn callback
                set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            end

            % Save Fig file
            if saveFigures==1; saveas(hFig, fullfile(savePathSubFigs, ...
                    ['isiL5_Sess_' num2str(SessNumb(n)) '_' ...
                    L5cells{ii} '.jpg']));
            end

        end

        dens_Go_L5 = structfun(@(x) x.density.CData, hdl_Go_L5,...
            "UniformOutput",false);
        dens_Go_L5 = struct2cell(dens_Go_L5);
        dens_Go_L5 = cat(3,dens_Go_L5{:});

        dens_NC_L5 = structfun(@(x) x.density.CData, hdl_NC_L5,...
            "UniformOutput",false);
        dens_NC_L5 = struct2cell(dens_NC_L5);
        dens_NC_L5 = cat(3,dens_NC_L5{:});

        if n==1
            isih_L5 = mean(cat(3, isih_Go.avg(ismember(isih_Go.label, ...
                L5cells(selL5unitsBasedonSpking{n})),:), ...
                isih_NC.avg(ismember(isih_Go.label, ...
                L5cells(selL5unitsBasedonSpking{n})),:)), 3, ...
                'omitnan')./size(spksGo_Target.trialtime, 1);
            
            isidens_L5 = mean(cat(4, dens_Go_L5, dens_NC_L5)./ ...
                size(spksGo_Target.trialtime, 1),4,'omitnan');
        else
            isih_L5 = cat(1, isih_L5, ...
                mean(cat(3, isih_Go.avg(ismember(isih_Go.label, ...
                L5cells(selL5unitsBasedonSpking{n})),:), ...
                isih_NC.avg(ismember(isih_Go.label, ...
                L5cells(selL5unitsBasedonSpking{n})),:)), 3, ...
                'omitnan')./size(spksGo_Target.trialtime, 1));

            isidens_L5 = cat(3, isidens_L5, ...
                mean(cat(4, dens_Go_L5, dens_NC_L5)./ ...
                size(spksGo_Target.trialtime, 1), 4,'omitnan'));
        end

        if noShowFigures == 1; close all; end
    end

    %% L3
    if ~isempty(L3cells(selL3unitsBasedonSpking{n}))
        clearvars hdl_Go_L3 hdl_NC_L3 dens_Go_L3 dens_NC_L3

        for ii=selL3unitsBasedonSpking{n}

            if noShowFigures == 1
                hFig = figure('Visible', 'off','Units', ...
                    'inches','Position',[0.05 0.05 10 5]);
            else
                hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
            end

            cfg = [];
            cfg.spikechannel = L3cells{ii};
            cfg.interpolate = 5;
            cfg.window = 'gausswin';
            cfg.winlen = 0.005;
            cfg.scatter = 'no';
            cfg.colormap = jet(300);

            subplot(1,2,1)
            [~, hdl_Go_L3.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_Go);

            subplot(1,2,2)
            [~, hdl_NC_L3.(cfg.spikechannel)] = ft_spike_plot_isireturn(cfg,isih_NC);

            if noShowFigures == 1
                % Set CreateFcn callback
                set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
            end

            % Save Fig file
            if saveFigures==1; saveas(hFig, fullfile(savePathSubFigs, ...
                    ['isiL3_Sess_' num2str(SessNumb(n)) '_' ...
                    L3cells{ii} '.jpg']));
            end

        end

        dens_Go_L3 = structfun(@(x) x.density.CData, hdl_Go_L3,...
            "UniformOutput",false);
        dens_Go_L3 = struct2cell(dens_Go_L3);
        dens_Go_L3 = cat(3,dens_Go_L3{:});

        dens_NC_L3 = structfun(@(x) x.density.CData, hdl_NC_L3,...
            "UniformOutput",false);
        dens_NC_L3 = struct2cell(dens_NC_L3);
        dens_NC_L3 = cat(3,dens_NC_L3{:});

        if n==1
            isih_L3 = mean(cat(3, isih_Go.avg(ismember(isih_Go.label, ...
                L3cells(selL3unitsBasedonSpking{n})),:), ...
                isih_NC.avg(ismember(isih_Go.label, ...
                L3cells(selL3unitsBasedonSpking{n})),:)), 3, ...
                'omitnan')./size(spksGo_Target.trialtime, 1);

            isidens_L3 = mean(cat(4, dens_Go_L3, dens_NC_L3)./ ...
                size(spksGo_Target.trialtime, 1),4,'omitnan');
        else
            isih_L3 = cat(1, isih_L3, ...
                mean(cat(3, isih_Go.avg(ismember(isih_Go.label, ...
                L3cells(selL3unitsBasedonSpking{n})),:), ...
                isih_NC.avg(ismember(isih_Go.label, ...
                L3cells(selL3unitsBasedonSpking{n})),:)), 3, ...
                'omitnan')./size(spksGo_Target.trialtime, 1));

            isidens_L3 = cat(3, isidens_L3, ...
                mean(cat(4, dens_Go_L3, dens_NC_L3)./ ...
                size(spksGo_Target.trialtime, 1),4,'omitnan'));
        end

        if noShowFigures == 1; close all; end
    end

    %% L6
%     if ~isempty(L6cells)
% 
%         for ii=1:length(L6cells)
%             if noShowFigures == 1
%                 hFig = figure('Visible', 'off','Units', 'inches','Position',[0.05 0.05 10 5]);
%             else
%                 hFig = figure('Units', 'inches','Position',[0.05 0.05 10 5]);
%             end
% 
%             cfg = [];
%             cfg.spikechannel = L6cells{ii};
%             cfg.interpolate = 5;
%             cfg.window = 'gausswin';
%             cfg.winlen = 0.005;
%             cfg.scatter = 'no';
%             cfg.colormap = jet(300);
% 
%             %             sgtitle(sprintf('Neuron %d', ii))
%             subplot(1,2,1)
%             ft_spike_plot_isireturn(cfg,isih_Go);
%             %             title('Go trials')
% %             set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
% 
%             subplot(1,2,2)
%             ft_spike_plot_isireturn(cfg,isih_NC);
%             %             title('NC trials')
% %             set(gca, 'box', 'off','linewidth',1,'fontsize',12,'fontweight','bold')
% 
%             if noShowFigures == 1
%                 % Set CreateFcn callback
%                 set(hFig, 'CreateFcn', 'set(gcbo,''Visible'',''on'')');
%             end
% 
%             % Save Fig file
%             if saveFigures==1; saveas(hFig, fullfile(savePathSubFigs, ...
%                     ['isiL6_Sess_' num2str(SessNumb(n)) '_' L6cells{ii} '.jpg'])); end
% 
%             if noShowFigures == 1; close all; end
% 
%         end
% 
%     end

end

%% mean dist across cells, and trials

isih_L3_mean = mean(isih_L3,'omitnan');
isih_L5_mean = mean(isih_L5,'omitnan');

isidens_L3_mean = mean(isidens_L3,3,'omitnan');
isidens_L5_mean = mean(isidens_L5,3,'omitnan');

%% --- save 

% -- min isi of all cells
save(fullfile(savePathSubFolder,'isi_min.mat'), 'minISI_allcells')

% -- mean ISI distribution 
save(fullfile(savePathSubFolder,'mean_isi_dist.mat'), 'isih_L3_mean', ...
    'isih_L5_mean', 'isih_L3', 'isih_L5', 'isidens_L3', 'isidens_L5')

% -- mean pre-target window length per session and trial type
save(fullfile(savePathSubFolder,'mean_preTarget_wind.mat'), 'mean_preTarget_wind')

close all
