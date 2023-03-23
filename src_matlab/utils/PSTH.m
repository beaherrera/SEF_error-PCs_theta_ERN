function [spk_count, psth] = PSTH(spiketimes, timesp)
% peri-stimulus time histogram
% Input:
%       spiketimes: spike train times
%       timesp: time intervals or histogram edges for computing the spike
%       rate
% Output:
%       spk_count: number of spikes per bin
%       psth: peri-stimulus time histogram

[spk_count, ~]= histcounts(spiketimes, timesp);
psth = spk_count./(diff(timesp));

end