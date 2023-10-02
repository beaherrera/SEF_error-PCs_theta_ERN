
# Simulated Data

## Data description

Sipiking activity and local field potentials (LFPs) obtained from biophysical simulations.

Processed LFPs are available through the Open Science Framework (OSF): https://osf.io/xs9nq/; location: OSF Storage folder data/sim_data/processed_lfps. Please, download and add the 'processed_lfps' folder to this directory before running the analysis scripts. The MATLAB project is already configured to add the data to the MATLAB path.


### .mat files
- events_timing_*.mat: time onset of the trial events (saccade and target onset) considered for each simulated error (NC - noncanceled) and correct (Go) trial. Saccade times are relative to the target onset.
- *_isi_dist_sim.mat: inter-spike interval distribution of the simulated neurons.
- *_spk_ft_struct_sim.mat: fieldtrip data structure containing the spiking activity of the simulated neurons during the whole simulation.
- *_spk_rate_sim.mat: spike rate of the simulated neurons.
- *_spk_tone_ft_struct_sim.mat: fieldtrip data structure containing the spiking activity of the simulated neurons until tone onset.
- tone_times.mat: time onset of the feedback tones considered for each simulated error (NC - noncanceled) and correct (Go) trial. Tones times are relative to the target onset. 