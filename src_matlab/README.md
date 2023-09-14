
Main analysis codes for generating the figures in Herrera et al. 2023.

## Description
- [ana_monkeys_data](ana_monkeys_data): scripts for analyzing the spiking activity and local field potentials recorded from the macaque monkeys while performing the stop-signal task. 
- [ana_sim_data](ana_sim_data): scripts for analyzing the simulated data (e.i., pyramidal error neurons spiking activity and their evoked local field potentials).
- [data/head_model_NMTv2](data\head_model_NMTv2): head model.
    - [lead_fields](data\head_model_NMTv2\lead_fields): lead field matrices used in the SEF EEG forward model.
    - [surfaces](data\head_model_NMTv2\surfaces): surfaces of the head model.
    - [SEF_vert_and_orient.mat](data\head_model_NMTv2\SEF_vert_and_orient.mat): location and orientation of SEF dipoles.
- [data/monkeys_data](data\monkeys_data): processed experimental data. Check [README.md](data\monkeys_data\README.md) in this folder for more details.
- [data/sim_data](data\sim_data): simulated data. Check [README.md](data\sim_data\README.md) in this folder for more details.
- [make_figs](make_figs): codes for making the figures in Herrera et al. 2023.
- [utils](utils): toolboxes used in the analyses.
- [Matlab.prj](Matlab.prj): Matlab project.
- [plot_headmodel_surfaces.m](plot_headmodel_surfaces.m): run to plot the monkey's head model used in the EEG forward model.
- [project_startup.m](project_startup.m): Matlab project startup script, initializes fieldtrip and brainstorm without the GUI.

## Dependencies
All dependencies EXCEPT Brainstorm are included in the Matlab project and will be added to the Matlab path after running Matlab.prj. 
- MATLAB (release 2021b, [The MathWorks](https://www.mathworks.com/?s_tid=gn_logo)). Toolboxes: Signal Processing Toolbox, Statistics and Machine Learning Toolbox, and Wavelet Toolbox.
- [CSDplotter](src_matlab\utils\CSDplotter-0.1.1) [[Pettersen et al. 2006](10.1016/j.jneumeth.2005.12.005)]
- [FieldTrip](https://www.fieldtriptoolbox.org/) (release [fieldtrip-20221126](src_matlab\utils\fieldtrip-20221126)). The fieldtrip release is already included in the Matlab project and is added to the Matlab path when you run Matlab.prj. The function [ft_spike_plot_isireturn.m](utils\fieldtrip-20221126\contrib\spike\ft_spike_plot_isireturn.m) was modified to output figure handle of the ISI heatmap (Figs 2E & 3E). Using your preferred fieldtrip version without including this modification will result in an error when running [cal_isi_dist_simL3errorneurons.m](src_matlab\ana_sim_data\cal_isi_dist_simL3errorneurons.m) and [cal_isi_dist_simL5errorneurons.m](src_matlab\ana_sim_data\cal_isi_dist_simL5errorneurons.m).
- [Brainstorm](https://neuroimage.usc.edu/brainstorm/Introduction). This toolbox is NOT INCLUDED in the Matlab project. You must download and install it before running Matlab.prj. 