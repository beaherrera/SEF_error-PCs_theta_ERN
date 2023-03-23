%% cal ERN and Pe magnitude as a func of the cortical column diameter

clear
clc

tspan = -500:1000; % time span in ms

ERN = median([184 224]); % ERN peak time in recorded EEG
Pe = round(median([302 327])); % Pe peak time in recorded EEG

%% perpendicular sessions number

Sess_EuP1 = 14:19; % Eu, site P1
Sess_XP2P3 = 20:29; % X, 20-25 site P2, 26-29 site P3

SessNumb = [Sess_EuP1, Sess_XP2P3];

%% load CSDs

load("data\monkeys_data\lfp_csd\avg_lfp_csd_sess.mat", 'avg_LFP')
avg_LFP_sess_Go = structfun(@(x) x.avg_LFP_Go, avg_LFP, ...
    'UniformOutput',false);
avg_LFP_sess_Go = struct2cell(avg_LFP_sess_Go);
avg_LFP_sess_Go = cat(3, avg_LFP_sess_Go{:});
grand_avg_LFP_Go_Eu = mean(avg_LFP_sess_Go(:,:,ismember(SessNumb, Sess_EuP1)), ...
    3, 'omitnan');

avg_LFP_sess_NC = structfun(@(x) x.avg_LFP_NC, avg_LFP, ...
    'UniformOutput',false);
avg_LFP_sess_NC = struct2cell(avg_LFP_sess_NC);
avg_LFP_sess_NC = cat(3, avg_LFP_sess_NC{:});
grand_avg_LFP_NC_Eu = mean(avg_LFP_sess_NC(:,:,ismember(SessNumb, Sess_EuP1)), ...
    3, 'omitnan');

load("data\head_model_NMTv2\surfaces\electrodes.mat")

%% load lead field matrices

selected_electrodes = [6 11 15 31];

% average reference
Ne = 4;
Hn = eye(Ne) - (ones(Ne,1)*ones(Ne,1)')/Ne;

load("data\head_model_NMTv2\lead_fields\leadFieldDipBEMVert12454.mat", 'Keoo')
Kedip_left = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load("data\head_model_NMTv2\lead_fields\leadFieldDipBEMVert12535.mat", 'Keoo')
Kedip_right = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load('data\head_model_NMTv2\lead_fields\leadFieldQuadBEMVert12454.mat', 'Keoo')
Kequad_left = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load('data\head_model_NMTv2\lead_fields\leadFieldQuadBEMVert12535.mat', 'Keoo')
Kequad_right = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load('data\head_model_NMTv2\lead_fields\leadFieldMonBEMVert12454Cdip.mat', 'Keoo')
Kemon_left = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load('data\head_model_NMTv2\lead_fields\leadFieldMonBEMVert12535Cdip.mat', 'Keoo')
Kemon_right = Hn*Keoo(selected_electrodes,:); clearvars Keoo

load("data\head_model_NMTv2\SEF_vert_and_orint.mat", 'orientation')

%% Eu's EEG

load("data\monkeys_data\ERN_Eu_4ch_SessAvg.mat", 'ERN_Go', 'ERN_NC')

ts = -299:(900-300);

ERN_Go_avg = squeeze(mean(ERN_Go, 1, 'omitnan'))';
ERN_Go_avg = Hn*ERN_Go_avg;

ERN_NC_avg = squeeze(mean(ERN_NC, 1, 'omitnan'))';
ERN_NC_avg = Hn*ERN_NC_avg;

%% normalization

norm_factor_EEG = max(abs([ERN_Go_avg(:,(ts>=-150 & ts <=500)) ...
    ERN_NC_avg(:,(ts>=-150 & ts <=500))]), [], 2, 'omitnan');

ERN_Go_avg = ERN_Go_avg./norm_factor_EEG;
ERN_NC_avg = ERN_NC_avg./norm_factor_EEG;

%% calculate CSD

% -- Parameters CSD calculation
Ne = 16; % number of electrodes in the shank
a = 1e-6; % [mm] position of the first electrode
elec_spacing = 0.15; % [mm] electrode spacing
ze = a:elec_spacing:((Ne-1)*elec_spacing + a); % position of the electrodes
% along z
el_pos = ze*1e-3;  % [m] electrode positions with respect to the pia surface
cond = 0.33; %[S/m] gray matter conductance 0.4
cond_top = 0; % conductance outside the gray matter (above pia matter)
this_tol = 1e-6;  % tolerance
gauss_sigma = 0.1e-3;   %[m] Gaussian filter std
filter_range = 5*gauss_sigma; % numeric filter must be finite in extent
diam_list = (1:6)*1e-3; % [m] cylinder diameter diam_list(ii)*1e-3; %

%% loop over diameters

ERN_exp_CSD = zeros(size(diam_list));
Pe_exp_CSD = zeros(size(diam_list));

for diam = diam_list

    % -- calculate cubic splines
    Fcs = F_cubic_spline(el_pos,diam,cond,cond_top);

    % -- get CSD
    [zs_Go,CSD_cs_Go] = make_cubic_splines(el_pos, ...
        grand_avg_LFP_Go_Eu, Fcs);
    if ~isempty(gauss_sigma) && gauss_sigma~=0 %filter iCSD
        [zs_Go,CSD_cs_Go] = gaussian_filtering(zs_Go, ...
            CSD_cs_Go, gauss_sigma, filter_range);
    end
    grand_avg_iCSD_Go = CSD_cs_Go; % [nA/mm3] current source density
    zs_Go = zs_Go*1e3; % [mm]

    [zs_NC,CSD_cs_NC] = make_cubic_splines(el_pos, ...
        grand_avg_LFP_NC_Eu, Fcs);
    if ~isempty(gauss_sigma) && gauss_sigma~=0 %filter iCSD
        [zs_NC,CSD_cs_NC] = gaussian_filtering(zs_NC, ...
            CSD_cs_NC, gauss_sigma, filter_range);
    end
    grand_avg_iCSD_NC = CSD_cs_NC; % [nA/mm3] current source density
    zs_NC = zs_NC*1e3; % [mm]

    %% get dipoles
    rc = diam.*1e3/2; % cortical column diameter in mm

    d_expCSD_Go = cal_dip_CSD(zs_NC', grand_avg_iCSD_Go, rc); % nA*m
    d_expCSD_NC = cal_dip_CSD(zs_NC', grand_avg_iCSD_NC, rc); % nA*m

    %% get monopole

    mon_expCSD_Go = sum(mean(abs(diff(zs_Go)))*grand_avg_iCSD_Go*(pi*(rc^2)), 1); % nA
    mon_expCSD_NC = sum(mean(abs(diff(zs_Go)))*grand_avg_iCSD_NC*(pi*(rc^2)), 1); % nA

    %% get quadrupoles

    quad_expCSD_Go = cal_quad_CSD(zs_NC', grand_avg_iCSD_Go, rc).*1e-6; % nA*m2
    quad_expCSD_NC = cal_quad_CSD(zs_NC', grand_avg_iCSD_NC, rc).*1e-6; % nA*m2

    %% calculate EEG
    %% dipolar component

    EEG_dip_exp_Go = (Kedip_right*(orientation(2,:)') + ...
        Kedip_left*(orientation(1,:)'))*(d_expCSD_Go.*1e-3); % uV
    
    EEG_dip_exp_NC = (Kedip_right*(orientation(2,:)') + ...
        Kedip_left*(orientation(1,:)'))*(d_expCSD_NC.*1e-3); % uV

    %% quadripolar component

    Q = blkdiag(zeros(3,1), zeros(3,1), [0; 0; 1]);

    Kequad = (Kequad_right*Q + Kequad_left*Q);
    Kequad = Kequad(:,end);

    EEG_quad_exp_Go = Kequad*(quad_expCSD_Go.*1e-3); % uV
    EEG_quad_sim_Go = Kequad*(quad_simCSD_Go.*1e-3); % uV

    EEG_quad_exp_NC = Kequad*(quad_expCSD_NC.*1e-3); % uV
    EEG_quad_sim_NC = Kequad*(quad_simCSD_NC.*1e-3); % uV

    %% monopolar component

    EEG_mon_exp_Go = (Kemon_left + Kemon_right)*(mon_expCSD_Go.*1e-3); % uV
    
    EEG_mon_exp_NC = (Kemon_left + Kemon_right)*(mon_expCSD_NC.*1e-3); % uV
    
    %% sum

    EEG_exp_Go = EEG_mon_exp_Go + EEG_quad_exp_Go + EEG_dip_exp_Go;
    EEG_exp_NC = EEG_mon_exp_NC + EEG_quad_exp_NC + EEG_dip_exp_NC;

    EEG_exp_Go =  EEG_exp_Go - mean(EEG_exp_Go(1,tspan>=-150 & tspan<=-100));
    EEG_exp_NC =  EEG_exp_NC - mean(EEG_exp_NC(1,tspan>=-150 & tspan<=-100));

    %% find ERN and Pe peak

    ERN_exp_CSD(ismember(diam_list, diam)) = EEG_exp_NC(1, tspan==ERN) -...
        EEG_exp_Go(1, tspan==ERN);
    Pe_exp_CSD(ismember(diam_list, diam)) = EEG_exp_NC(1, tspan==Pe) -...
        EEG_exp_Go(1, tspan==Pe);

end

%% linear regression
[f_ERN, gof_ERN] = fit(diam_list'.*1e3, ERN_exp_CSD','poly1');
[f_Pe, gof_Pe] = fit(diam_list'.*1e3, Pe_exp_CSD','poly1');

%% 

font = 10;
width = 1.5;

ylim_min = -3;
ylim_max = 0;

figure('Units', 'inches','Position',[0 0 6 2]);
tiledlayout(1, 2,'TileSpacing','Compact','Padding','Compact');

nexttile
h1 = plot(f_ERN, '-', diam_list.*1e3, ERN_exp_CSD,'r.');
h1(1).MarkerSize = 10;  
h1(2).Color = [0 204 0]./255;
xlabel('Column diameter (mm)')
ylabel({'Eu CSD EEG', 'mag at ERN (\muV)'})
ylim([ylim_min ylim_max])
legend({'Data',sprintf('Fitted Curve \n Slope = %.2f \n  R^2=%.2f', ...
    f_ERN.p1, gof_ERN.rsquare)}, ...
    'Box','off','Location','southwest')
set(gca, 'box', 'off','linewidth',width,'fontsize',font,'fontweight','bold')

nexttile
h2 = plot(f_Pe, '-', diam_list.*1e3, Pe_exp_CSD,'r.');
h2(1).MarkerSize = 10;  
h2(2).Color = [0 204 0]./255;
xlabel('Column diameter (mm)')
ylabel({'Eu CSD EEG', 'mag at Pe (\muV)'})
legend({'Data',sprintf('Fitted Curve \n Slope = %.2f \n  R^2=%.2f', ...
    f_Pe.p1, gof_ERN.rsquare)}, ...
    'Box','off','Location','southwest')
ylim([ylim_min ylim_max])
set(gca, 'box', 'off','linewidth',width,'fontsize',font,'fontweight','bold')
