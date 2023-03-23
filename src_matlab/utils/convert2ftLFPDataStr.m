function ft_data = convert2ftLFPDataStr(data, Fs_lfp, Fs_Mua, TrialStartSamples)
%% converts data to FieldTrip Raw data structure
% Inputs:
%       data: Nchans x Ntrials cell-array in mVolts
%       Fs_lfp: lfp sampling rate in Hz
%       Fs_Mua: mua sampling rate in Hz
%       TrialStartSamples: trial start in recording samples, vector
% Outputs:
%       ft_data: fieldtrip LFP data structure
%           Fields of the structure:
%               - trials: lfp separated by trials | cell-array nTrials
%               - times: time span of the trials | cell array nTrials
%               - labels: name of the electrodes in the linear probe | 
%                   cell array nElectrodes
%               - fsample: lfp sampling rate
%               - sampleinfo: trials information (begin, end) in 
%                   recording samples 
%               - dimord = '{rpt}_chan_time': dimension of the lfp data
%               - cfg.trl = [trl; [beginSamp endSamp offset]]; trial 
%                   definition for the raw data
%               - hdr: header info structure to link to spiking data | 
%                   essential for spk-lfp
%               analyses
%                   - hdr.nChans: number of electrodes 
%                   - hdr.Fs: lfp sampling rate
%                   - hdr.label: electrodes' label
%                   - hdr.nTrials: number of trials
%                   - hdr.nSamplesPre
%                   - hdr.FirstTimeStamp
%                   - hdr.TimeStampPerSample: number of Mua samples per 
%                       lfp time point
% Copyright 2021 @ Beatriz Herrera
%%

nTrials = 1; % number of selected trials 
nEle = size(data,1); % number of electrodes

datas = data(:,TrialStartSamples:end).*1e-3; % converting from mV -> V

% -- low-pass filtering the data at 100 Hz using a two-pass fourth order 
% Butterworth filter.
N = 4; % filter order
fsample = Fs_lfp; % lfp sampling rate
type = 'but'; % Butterworth filter
dir='twopass'; % two-pass filter 
instabilityfix = 'reduce';
filterFreq = 100; % cutoff frequency

trials = ft_preproc_lowpassfilter(datas, fsample, ...
    filterFreq, N, type, dir, instabilityfix); % filtered data
trials = {trials};

timeGenFun = @(x) 0:(1/Fs_lfp):((size(x,2)-1)/Fs_lfp); % assumes 
% x -> nchs x ntrials | timespan generator function
times = cellfun(timeGenFun, trials, 'UniformOutput', false); % in 
% seconds | generates the trials timespan 
lenghtTrials = cellfun(@(x) size(x,2), trials); % computes the length of
% each trial in time points
sampleInfo = [TrialStartSamples TrialStartSamples+lenghtTrials-1]; 

% generating the electrodes labels
labels = {'e1'};
for ii=2:nEle
    labels = cat(1,labels,{sprintf('e%.0f',ii)});
end
ft_data.trial = trials'; % lfp separated by trials | cell-array nTrials
ft_data.time = times'; % time span of the trials | cell array nTrials
ft_data.label = labels; % name of the electrodes in the linear 
% probe | cell array nElectrodes 
ft_data.fsample = Fs_lfp; % lfp sampling rate
ft_data.sampleinfo = sampleInfo; % trials information (begin, end)
% in recording samples 
ft_data.dimord = '{rpt}_chan_time'; % dimension of the lfp data

% - header infor structure to link to spiking data | essential for spk-lfp
% analyses
ft_data.hdr.nChans = nEle; % number of electrodes 
ft_data.hdr.Fs = Fs_lfp; % lfp sampling rate
ft_data.hdr.label = labels; % electrodes' label
ft_data.hdr.nTrials = nTrials; % number of trials
ft_data.hdr.nSamplesPre = 0;
ft_data.hdr.FirstTimeStamp = 0;
ft_data.hdr.TimeStampPerSample = Fs_Mua/ft_data.fsample;  % number
% of Mua samples per lfp time point

% - specifying the trial definition for the raw data
beginSamp = zeros(nTrials,1);
offset = zeros(nTrials,1);
endSamp = lenghtTrials;
trl = [];
ft_data.cfg.trl = [trl; [beginSamp endSamp offset]];

end