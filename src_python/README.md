## Biophysical Simulations in Herrera et al. 2023

Codes for replicating simulations in Herrera et al. 2023.

Note: mod files must be compiled before running simulation scripts.

### Description
- [cell_models](cell_models): cells model files.
- [Data](Data): data files.
- [generate_tone_times.py](generate_tone_times.py): run to generate the tone times for each simulated trial.
- [hayetalmodel_active_declarations.py](hayetalmodel_active_declarations.py): biophysics of the Hay et al. 2011 model (ModelDB, accession #139653, “cell #1”).
- make_fig5_data_*.py: run to generate the data for creating Fig. 5.
- [make_synapses.py](make_synapses.py): functions package to insert synapses into the L3 and L5 PC models.
- [modhayetalmodel_active_declarations.py](modhayetalmodel_active_declarations.py): biophysics of the Hay et al. 2011 model, incorporating the modifications of voltage-gated calcium channel densities as in Shai et al. (2015) and Ih channel density distribution as in Labarrera et al. (2018) (Leleo and Segev 2021).
- [params_synapses.py](params_synapses.py): synaptic models parameters, and number of synapses per region per cell type.
- [PopulationL3PCsEyal.py](PopulationL3PCsEyal.py) and [PopulationL5PCsHay.py](PopulationL5PCsHay.py): Python classes to create populations of N (unconnected) L3 and L5 PCs discribed by Eyal et al. (2018) and Hay et al. (2011) model, respectively.
- run_*_sim_Go.py: run to simulate the L3 or L5 PCs population for correct (Go) trials. Check the script for the input parameters.
- run_*_sim_NC.py: run to simulate the L3 or L5 PCs population for error (NC) trials. Check the script for the input parameters.
- [utils](utils): util functions.

### Dependencies
