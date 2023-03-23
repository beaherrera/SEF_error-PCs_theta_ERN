function rate_data = cal_spkrate_ftspkstr(spikeTrain, timesp, w, flt)
% calculate peri-stimulus time histogram and firing rate
% Input:
%       spikeTrain: fieldtrip structure of spike data
%       timesp: histogram edges to creater spike_count hist
%       w: SD of the Gaussian windows used to smooth the spike rate
%       flt: boolean, if True or 1, Gaussian filter through smoothdata
%       function is used instead of conv function.
%
% Output:
%       rate_data: structure
%           - psthi: peri-stimulus time histogram of each neuron per trial.
%           cell-array 1 x nNeurons, each cell element is an array
%           nTrials x nBins
%           - pstht: peri-stimulus time histogram per neuron across trials.
%           cell-array 1 x nNeurons
%           - spk_counti: hist spk_count of each neuron per trial.
%           cell-array 1 x nNeurons, each cell element is an array
%           nTrials x nBins
%           - firingRatei: spike rate of each neuron per trial.
%           cell-array 1 x nNeurons, each cell element is an array
%           nTrials x nBins 
%           - firingRatet: mean spike rate per neuron across trials. 
%           cell-array 1 x nNeurons
%
% @Copyright 2021- , Beatriz Herrera
% Created Sep 27, 2021
%
%%

NumTrials = size(spikeTrain.trialtime,1); % number of trials
spkstrialNum = spikeTrain.trial; % trial number to which the spike times in
                            % spkstrial spkstimes belong. 
spkstimes = spikeTrain.time; % spike times per neuron in all trials.

% pre-allocating 
psthi = cell(NumTrials,length(spikeTrain.label));
spk_counti = cell(NumTrials,length(spikeTrain.label));
firingRatei = cell(NumTrials,length(spikeTrain.label));

for ii=1:NumTrials
    sepTrialsFun = @(x, y) x(y==ii); % separate spike times into 
                                    %  the corresponding trials 
    [spk_counti(ii,:), psthi(ii,:)] = cellfun(@(x, y) PSTH(sepTrialsFun(x, y), timesp),...
        spkstimes, spkstrialNum, 'UniformOutput', false);
    firingRatei(ii,:) = cellfun(@(x) firingrate(x, w, flt),...
        psthi(ii,:), 'UniformOutput', false); 
end

psthi = arrayfun(@(col) cell2mat(psthi(:,col)),(1:size(psthi,2))',...
    'UniformOutput', false)'; % psth for each trial per each neuron
spk_counti = arrayfun(@(col) cell2mat(spk_counti(:,col)),(1:size(spk_counti,2))',...
    'UniformOutput', false)'; % spike count for each trial per each neuron
firingRatei = arrayfun(@(col) cell2mat(firingRatei(:,col)),(1:size(firingRatei,2))',...
    'UniformOutput', false)'; % spike rate for each trial per each neuron

pstht = cellfun(@(x) mean(x,'omitnan'), psthi,'UniformOutput', false); % psth across trials per neuron
% spk_countt = cellfun(@nanmean, spk_counti,...
%     'UniformOutput', false);
firingRatet = cellfun(@(x) mean(x,'omitnan'),...
    firingRatei,'UniformOutput', false); % mean spike rate across trials per neuron

% create structure to store data
rate_data.psthi = psthi;
rate_data.pstht = pstht;
rate_data.spk_counti = spk_counti;
rate_data.firingRatei = firingRatei;
rate_data.firingRatet = firingRatet;

end
