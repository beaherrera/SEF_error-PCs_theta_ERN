function ft_data = convSimLFPs2ftRawDataStr(data, Fs_lfp, Fs_MUA)
%% converts data to FieldTrip Raw data structure
% Inputs:
%       data: Nchans x Ntrials cell-array in mVolts
%       Fs_lfp: lfp sampling rate in Hz
% Outputs:
%       ft_data: fieldtrip LFP data structure
%           Fields of the structure:
%               - trials: lfp separated by trials | cell-array nTrials
%               - times: time span of the trials | cell array nTrials
%               - labels: name of the electrodes in the linear probe | cell array nElectrodes
%               - fsample: lfp sampling rate
%               - sampleinfo: trials information (begin, end) in recording samples 
%               - dimord = '{rpt}_chan_time': dimension of the lfp data
%               - cfg.trl = [trl; [beginSamp endSamp offset]]; trial definition for the raw data
%               - hdr: header info structure to link to spiking data | essential for spk-lfp
%               analyses
%                   - hdr.nChans: number of electrodes 
%                   - hdr.Fs: lfp sampling rate
%                   - hdr.label: electrodes' label
%                   - hdr.nTrials: number of trials
%                   - hdr.nSamplesPre
%                   - hdr.FirstTimeStamp
%                   - hdr.TimeStampPerSample: number of Mua samples per lfp time point
% Copyright 2022 @ Beatriz Herrera
%%

nTrials = size(data, 2); % number of trials 
nEle = size(data, 1); % number of electrodes

% concatenating across electrodes and converting from mV -> V
trials = arrayfun(@(col)[(data{:,col})]'.*1e-3,...
    (1:size(data,2))','UniformOutput', false); % volts

timeGenFun = @(x) 0:(1/Fs_lfp):((size(x,2)-1)/Fs_lfp); % assumes x -> nchs
% x ntrials | timespan generator function
times = cellfun(timeGenFun, trials, 'UniformOutput', false); % in seconds | generates the trials timespan 
lenghtTrials = cellfun(@(x) size(x,2), trials); % computes the length of each trial in time points
sampleInfo = [(1:lenghtTrials(1,1):(nTrials*lenghtTrials(1,1)))' ...
    (lenghtTrials(1,1):lenghtTrials(1,1):(nTrials*lenghtTrials(1,1)))'];

% generating the electrodes labels
labels = {'e01'};
for ii=2:nEle
    if ii<10
        labels = cat(1,labels,{sprintf('e0%.0f',ii)});
    else
        labels = cat(1,labels,{sprintf('e%.0f',ii)});
    end
end
ft_data.trial = trials'; % lfp separated by trials | cell-array nTrials
ft_data.time = times'; % time span of the trials | cell array nTrials
ft_data.label = labels; % name of the electrodes in the linear probe | cell array nElectrodes 
ft_data.fsample = Fs_lfp; % lfp sampling rate
ft_data.sampleinfo = sampleInfo; % trials information (begin, end) in recording samples 
ft_data.dimord = '{rpt}_chan_time'; % dimension of the lfp data

% - header info structure to link to spiking data | essential for spk-lfp
% analyses
ft_data.hdr.nChans = nEle; % number of electrodes 
ft_data.hdr.Fs = Fs_lfp; % lfp sampling rate
ft_data.hdr.label = labels; % electrodes' label
ft_data.hdr.nTrials = nTrials; % number of trials
ft_data.hdr.nSamplesPre = 0;
ft_data.hdr.FirstTimeStamp = 0;
ft_data.hdr.TimeStampPerSample = Fs_MUA/ft_data.fsample;  % number of Mua samples per lfp time point

% - specifying the trial definition for the raw data
beginSamp = sampleInfo(:,1);
offset = zeros(nTrials,1);
endSamp = sampleInfo(:,2);
trl = [];
ft_data.cfg.trl = [trl; [beginSamp endSamp offset]];

end