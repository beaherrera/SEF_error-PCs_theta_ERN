
# SEF Monkeys Data

## Data description

### folders

Available through the Open Science Framework (OSF): https://osf.io/xs9nq/; location: OSF Storage folder data/monkeys_data. All files are fieldtrip data structures. Please, download and add the data to this folder before running the analysis scripts. The MATLAB project is already configured to add the data to the MATLAB path.
- laminar_data: laminar field potentials and spiking activity recorded per session relative to target and saccade onset.
    - Session_*_Saccade.mat: LFPs relative to saccade onset for Go (correct) and NC (error) trials.
    - Session_*_Saccade_filt.mat: filtered LFPs relative to saccade onset for Go (correct) and NC (error) trials.
    - Session_*_Spks_Saccade.mat: spiking activity of all recorded single units relative to saccade onset for Go (correct) and NC (error) trials.
    - Session_*_Spks_Target.mat: spiking activity of all recorded single units relative to target onset for Go (correct) and NC (error) trials.
    - Session_*_Target_filt.mat: filtered LFPs relative to target onset for Go (correct) and NC (error) trials.
- laminar_power_saccade & laminar_power_target: laminar time-varying power maps for each frequency band (theta, alpha, beta, gamma) relative saccade and target onset, respectively.
    - sessions_*_power.mat: cell arrays dimensions: (electrodes x time x trials).

### .mat files
- avg_lfp_csd_sess.mat: averaged local field potential and current source density across trials for each session and trial type. Sessions 14-19: monkey Eu; Sessions 20-29: monkey X. tspan = -500:1000 ms relative to saccade onset.
    - avg_LFP.Session_*: averaged local field potentials. Dimensions: (electrode x time).
    - iCSD.Session_*.iCSD_*: current source density maps. Dimensions: (cortical depth x time).
    - iCSD.Session_*.zs_*: cortical depth in mm where the CSDs were estimated.
- eleAlignment.mat: electrodes aligment across sessions. Columns 1-6 == Sessions 14-19: monkey Eu; columns 7-16 == sessions 20-29: monkey X.
- ERN_Eu_4ch_SessAvg.mat: Monkey Eu's EEG (FpFz, Cz, F3, and F4 ) recorded -300 to 500 ms relative to saccade onset. Dimensions: (channel x time). FpFz - channel 1.
- grand average event-related potential of the error-related negativity (ERN) for monkey Eu. Dimensions: (electrode x time).
- fr_error_cells_saccade.mat and fr_error_cells_target.mat: spiking activity of L3 and L5 pyramidal error neurons relative to saccade and target onset, respectively. 
    - L*_mean_fr_*: mean spike rate across recorded neurons per layer and trial type.
    - L*_sem_fr_*: spike rate standard error of the mean across recorded neurons per layer and trial type.
    - mean_fr_*_*: mean spike rate across trials for each recorded neuron per layer and trial type.
    - ts_sess: time stamps relative to the trial events (saccade and target onset).
- grand_avg_lfp_csd.mat: grand average laminar local field potentials and current source density maps across all sessions per trial type (NC - error noncanceled trials and Go - correct go trials). 
    - grand_avg_LFP_*: grand average local field potentials. Dimensions: (electrode x time).
    - grand_iCSD_*: grand average current source density maps. Dimensions: (cortical depth x time).
    - tspan: time stamps relative to saccade onset.
    - zs_*: cortical depth in mm where the CSDs were estimated.
- idx_L*neurons_sac.mat: indices of selected L3 and L5 pyramidal error neurons for constraining the biophysical models.
- neuronsInfoAmirSteven.mat: information of the recorded neurons (e.g., name, depth, unit type (BS - broad spiking - or NS - narrow spiking unit), etc.).
- sess_avg_power.mat: laminar time-varying power for each frequency band (theta, alpha, beta and gamma) of the local field potentials for each monkey (Eu and X) and trial type (NC - error noncanceled trials and Go - correct go trials).