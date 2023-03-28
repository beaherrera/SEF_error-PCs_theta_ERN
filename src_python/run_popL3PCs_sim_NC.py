# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:50:03 2022

@author: Beatriz Herrera

This script simulates the activity of L3 pyramidal error neurons on Error trials.
Error trials == NC trials.

    * Generate data for Fig 2: `num_trials = 116` and `POPULATION_SIZE = 10`
    * Generate data for Figs 4, 6-8: `num_trials = 20` and `POPULATION_SIZE = 625`

Only these parameters were changed. All other parameters are the same
in all the simulations.

Default parameters are set for generating data for Fig 2.
    
"""
from __future__ import division
import gc
import sys
import os
import LFPy
import numpy as np
import scipy.stats as stats
from scipy import io
from mpi4py import MPI
from neuron import h
from os.path import join
from PopulationL3PCsEyal import PopulationL3PC
from utils import create_output_folder_L3Pop, gen_file_name_L3Pop


# #### MPI ####
# instantize the communication world
COMM = MPI.COMM_WORLD
# get the size of the communication world
SIZE = COMM.Get_size()
# get this particular processes' `rank` ID
RANK = COMM.Get_rank()


""" Functions """


def main():

    # Generate data for Fig 2: `num_trials = 116` and `POPULATION_SIZE = 10`
    # Generate data for Figs 4, 6-8: `num_trials = 20` and `POPULATION_SIZE = 625`

    dt = 2 ** -4  # ms | time step for simulation
    SEED = 12  # seed for random number generator
    num_trials = 116  # number of trials to be simulated
    trials_indices = np.arange(1, num_trials + 1)  # trial indices
    POPULATION_SIZE = 10  # number of L3 PCs in the population
    warmup_period = 500  # ms | warm up period, not icnluded in analyses
    sim_length = 10000  # ms | length of the simulation, excluding
    # warm up period

    trial_type = "NC_trial"
    child_folder = trial_type + "_Sept9_mpi"  # name of the folder
    # where results will be saved

    if "win32" in sys.platform:  # if running on windows local machine

        trials_ID = trials_indices
        job_0 = True

        main_path_folder = os.path.join(r"E:\results_L3PCsPopMky", child_folder)
        # path to the main folder where results will be saved

    else:  # if running on a cluster
        # distribute trials across slurm job array tasks

        # slurm job array index value
        task_ID = int(os.environ["SLURM_ARRAY_TASK_ID"])
        # total number of slurm job array tasks
        num_tasks = int(os.environ["SLURM_ARRAY_TASK_COUNT"])

        trials_ID = trials_indices[trials_indices % num_tasks == task_ID]  # trials
        # indices to be simulated in this task (SLURM_ARRAY_TASK_ID)

        if task_ID == 0:
            job_0 = True
        else:
            job_0 = False

        main_path_folder = os.path.join("results_L3PCsPopMky", child_folder)  # path
        # to the main folder where results will be saved

    # ----  Options for cell model and simulation:
    cellParameters = {
        "tstop": warmup_period + sim_length,  # ms | total simulation time
        "dt": dt,  # ms | simulation time step
        "cell_model": "cell0603_08_model_602",  # cell model to be used
        "morphology": "2013_03_06_cell03_789_H41_03.ASC",  # morphology file
        "save_neurons_data": True,  # save neuron data
        "plot_synapse_locs": False,  # plot synapse locations
        "show_plot_synapse_locs": False,  # show plot of synapse locations
        "show_lfp_plot": False,  # show LFP plot
    }

    # ----  Options for synapse configuration:
    # Case 0: no synapses in that region
    # Case 1: NMDA + AMPA + GABAA + GABAB
    # Case 2: NMDA + AMPA + GABAA
    # Case 3: NMDA + AMPA + GABAB
    # Case 4: NMDA + AMPA
    # Case 5: NMDA + GABAA + GABAB
    # Case 6: NMDA + GABAA
    # Case 7: NMDA + GABAB
    # Case 8: AMPA + GABAA + GABAB
    # Case 9: AMPA + GABAA
    # Case 10: AMPA + GABAB
    # Case 11: NMDA
    # Case 12: AMPA
    stimulusType = {
        "trial_type": trial_type,  # trial type being simulated. affects the
        # value of the time-locked stimulus.
        "synapse_type": "clustered",  # Opts: distributed or clustered
        "dend_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
        "apic_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
    }

    # load experimental distribution of target times
    if RANK == 0:  # only load once
        event_times = io.loadmat(join("Data", "event_times_rel2target.mat"))
        xk_target = event_times["edges_target_time"][:, 1:] - 0.5
        xk_saccade = event_times["edges_saccade_time"][:, 1:] - 0.5
        if stimulusType["trial_type"] == "Go_trial":
            pk_target = event_times["p_target_times_Go"]
            pk_saccadeT = event_times["p_saccade_times_Go"]
        else:
            pk_target = event_times["p_target_times_NC"]
            pk_saccadeT = event_times["p_saccade_times_NC"]

        target_dist = stats.rv_discrete(
            name="target_dist", values=(xk_target, pk_target)
        )
        saccade_dist = stats.rv_discrete(
            name="target_dist", values=(xk_saccade, pk_saccadeT)
        )

        target_times = np.tile(
            target_dist.rvs(size=(1, num_trials), random_state=SEED),
            (POPULATION_SIZE, 1),
        )
        saccade_times = np.tile(
            saccade_dist.rvs(size=(1, num_trials), random_state=SEED),
            (POPULATION_SIZE, 1),
        )

    else:  # other ranks
        target_times = None
        saccade_times = None

    # broadcast the target and saccade times to all ranks
    target_times = COMM.bcast(target_times, root=0)
    saccade_times = COMM.bcast(saccade_times, root=0)

    # ----  Options for time-locked input:
    inputParameters = {
        "warmup_period": warmup_period,  # length of the warmup period
        "target_times": target_times,  # target times
        "saccade_times": saccade_times,  # saccade times
        #
        "Gaussian_Input": True,  # turn-on time-locked inputs, otherwise:
        # only background inputs are inserted
        "Skewed_dist": True,  # use a normal gaussian or skewed Gaussian dist
        "a": -1,  # shape parameter for skewed Gaussian dist
        "mu": (
            (warmup_period + 216.6)
            if stimulusType["trial_type"] == "Go_trial"
            else (warmup_period + 298.6)
        ),  # ms, center of the Gaussian
        # defined as warup + target_time + saccade_time + delay
        "sigma": 141.6 if stimulusType["trial_type"] == "Go_trial" else 178.6,
        # ms, width of the Gaussian. 141.6: mean peak time of
        # Go trials firing rate relative to saccade. 178.6: NC trials
        "dend_syn": True,  # add time-locked input to dendrites
        "oblq_syn": False,  # add time-locked input to proximal apical dendrites
        "apic_syn": True,  # add time-locked input to distal apical dendrites
        "n_dend": 2 if stimulusType["trial_type"] == "Go_trial" else 4,
        # number of pre-synaptic spikes for basal inputs
        "n_apic": 3,  # number of pre-synaptic spikes for apical inputs
    }

    # - background input spike rates
    # ---- excitatory synapses ----
    r_AMPA_dend = 3.5
    r_NMDA_dend = r_AMPA_dend

    r_AMPA_apic = 2
    r_NMDA_apic = r_AMPA_apic

    # # ---- inhibitory synapses ----
    # r_GABAA_dend = 0
    # r_GABAB_dend = 0
    # r_GABAA_apic = 0
    # r_GABAB_apic = 0

    backgrd_inputs_rates = {
        "dend": dict(
            r_NMDA_dend=r_NMDA_dend,
            r_AMPA_dend=r_AMPA_dend,
            # r_GABAA_dend=r_GABAA_dend,
            # r_GABAB_dend=r_GABAB_dend,
        ),
        "apic": dict(
            r_NMDA_apic=r_NMDA_apic,
            r_AMPA_apic=r_AMPA_apic,
            # r_GABAA_apic=r_GABAA_apic,
            # r_GABAB_apic=r_GABAB_apic,
        ),
    }

    # ---- Define electrode geometry corresponding to a laminar probe:
    a = 0  # location of the first electrode relative to the grey matter
    # / CSF boundary
    electrode_spacing = 150  # inter-electrodes space
    Ne = 16  # number of electrodes inside the cortex
    # z position of the electrodes in microns
    z = -np.mgrid[a : (Ne * electrode_spacing + a) : electrode_spacing]
    electrodeParameters = {
        "x": np.zeros(z.size),
        "y": np.zeros(z.size),
        "z": z,
        "sigma": 0.33,  # S/m
        "method": "pointsource",  # method used to compute the LFP
    }

    rates_dend = {
        "r_AMPA": r_AMPA_dend,
        "r_NMDA": r_NMDA_dend,
        # "r_GABAA": r_GABAA_dend,
        # "r_GABAB": r_GABAB_dend,
    }
    rates_apic = {
        "r_AMPA": r_AMPA_apic,
        "r_NMDA": r_NMDA_apic,
        # "r_GABAA": r_GABAA_apic,
        # "r_GABAB": r_GABAB_apic,
    }

    if RANK == 0:  # only rank 0 saves the events timing
        file_name = gen_file_name_L3Pop(
            stimulusType, **{"rates_dend": rates_dend, "rates_apic": rates_apic}
        )
        data_folder = create_output_folder_L3Pop(
            main_path_folder, POPULATION_SIZE, stimulusType, inputParameters
        )

        save_events = {
            "target_times": target_times,  # target times
            "saccade_times": saccade_times,  # saccade times
        }

        if job_0:
            io.savemat(
                join(data_folder, ("events_timing_" + file_name + ".mat")), save_events
            )
    else:  # other ranks do not save the events timing
        data_folder = None
        file_name = None

    # ---- broadcast data_folder and file_name to all ranks
    data_folder = COMM.bcast(data_folder, root=0)
    file_name = COMM.bcast(file_name, root=0)

    # ---- run
    for runNumb in trials_ID:  # loop over trials

        # delete old sections from NEURON namespace
        h("forall delete_section()")
        LFPy.cell.neuron.h("forall delete_section()")

        # ---- INITIALIZE POPULATION
        population = PopulationL3PC(
            POPULATION_SIZE,
            cellParameters,
            electrodeParameters,
            stimulusType,
            backgrd_inputs_rates,
            inputParameters,
            runNumb,
            SEED,
            data_folder,
        )
        population.run(file_name)
        population.save_simData(file_name)
        # population.plotstuff(file_name)
        del population
        gc.collect()


if __name__ == "__main__":
    # load some required neuron-interface files
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")

    # population1 = main()

    main()

    # plt.plot(population.somav[0,5000:10000].T)
    # plt.plot(population1.somav[0,5000:10000].T)
