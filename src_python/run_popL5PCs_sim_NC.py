# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:50:03 2022

@author: Beatriz Herrera

This script simulates the activity of L5 pyramidal error neurons on Correct trials.
Correct trials == Go trials.

    * Generate data for Fig 3: `num_trials = 116` and `POPULATION_SIZE = 5`
    * Generate data for Figs 4, 6-8: `num_trials = 20` and `POPULATION_SIZE = 1000`

Only these parameters were changed. All other parameters are the same
in all the simulations.

Default parameters are set for generating data for Fig 3.

"""
from __future__ import division
import gc
import sys
import os
import LFPy
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from scipy import io
from mpi4py import MPI
from neuron import h
from os.path import join
from utils import create_output_folder_L5Pop, gen_file_name_L5Pop
from PopulationL5PCsHay import PopulationL5PC
from LFPy.inputgenerators import get_activation_times_from_distribution


# #### MPI ####
# instantize the communication world
COMM = MPI.COMM_WORLD
# get the size of the communication world
SIZE = COMM.Get_size()
# get this particular processes' `rank` ID
RANK = COMM.Get_rank()


""" Functions """


def main():

    # Generate data for Fig 3: `num_trials = 116` and `POPULATION_SIZE = 5`
    # Generate data for Figs 4, 6-8: `num_trials = 20` and `POPULATION_SIZE = 1000`

    dt = 2 ** -4  # ms | simulation time step
    SEED = 12  # seed for random number generator
    num_trials = 116  # number of trials to be simulated
    trials_indices = np.arange(1, num_trials + 1)
    POPULATION_SIZE = 5  # number of L5 PCs in the population
    warmup_period = 500  # ms | warm up period, not included in analyses
    sim_length = 10000  # ms | length of the simulation, excluding
    # warm up period

    trial_type = "Go_trial"
    child_folder = trial_type + "_Oct10_2022"  # name of the folder
    # where results will be saved

    if "win32" in sys.platform:  # if running on windows local machine

        trials_ID = trials_indices
        job_0 = True

        main_path_folder = os.path.join(
            r"D:\results_L5PCsPopMky", child_folder
        )  # path to the folder where results will be saved

    else:  # if running on the cluster
        # distribute trials across slurm job array tasks

        # slurm job array index value
        task_ID = int(os.environ["SLURM_ARRAY_TASK_ID"])
        # total number of tasks in the slurm job array
        num_tasks = int(os.environ["SLURM_ARRAY_TASK_COUNT"])

        trials_ID = trials_indices[trials_indices % num_tasks == task_ID]  # trials
        # indices to be simulated in this task (SLURM_ARRAY_TASK_ID)

        if task_ID == 0:
            job_0 = True
        else:
            job_0 = False

        main_path_folder = os.path.join("results_L5PCsPopMky", child_folder)  # path
        # to the folder where results will be saved

    # ----  Options for cell model and simulation:
    cellParameters = {
        "tstop": warmup_period + sim_length,  # ms | simulation time
        "dt": dt,  # ms | simulation time step
        "mod_Hay_model": True,  # use the modified Hay model
        "save_neurons_data": True,  # save data from each neuron
        "plot_synapse_locs": True,  # plot synapse locations
        "show_plot_synapse_locs": False,  # show plot of synapse locations
        "show_lfp_plot": True,  # show plot of LFP
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
        "oblq_synp": 0,  # case number within synp configuration options;
        # if 0, no synapses are inserted
        "apic_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
        # synaptic locations activated during baseline period
        "background_inputs": "dend",
        # baseline | "dend" == only basal inputs during baseline period
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

        target_dist = st.rv_discrete(name="target_dist", values=(xk_target, pk_target))
        saccade_dist = st.rv_discrete(
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

    # broadcast target and saccade times to all ranks
    target_times = COMM.bcast(target_times, root=0)
    saccade_times = COMM.bcast(saccade_times, root=0)

    # ----  Options for input configuration:
    inputParameters = {
        "target_times": target_times,  # target times
        "saccade_times": saccade_times,  # saccade times
        "warmup_period": warmup_period,  # length of the warmup period
        #
        # Inputs Definition
        "Gaussian_Input": True,  # time-locked activation of synaptic inputs
        # following Gaussian probability distribution
        "increase_baseline_inputs": True,  # if True, the mean of the Poisson
        # process activating background synapses is increased after target
        # onset
        #
        # Time-locked inputs definition
        # -- Time-locked synaptic inputs location
        "dend_syn": True,  # add time-locked input to basal dendrites
        "oblq_syn": False,  # add time-locked input to proximal apical dendrites
        "apic_syn": True,  # add time-locked input to distal apical dendrites
        "sep_hemitrees": True,  # separate synaptic inputs based on the
        # hemitree they belong to
        "hemitree1": True,  # activate synapses on hemitree 1 (compartments
        # [617,691])
        "hemitree2": True,  # activate synapses on hemitree 2 (compartments
        # [695,737])
        "hemitree3": True,  # activate synapses on hemitree 3 (compartments
        # [738,894])
        #
        # -- Gaussian Shape Parameters
        "a_dend": [0, 5],  # basal synapses | original value -5
        # apical
        "a_apic": {"a_apic": 2, "a_h1": [2], "a_h2": [2], "a_h3": [2]},
        #  synapses
        # "a_oblq": 0,  # oblique synapses
        # -- center of Gaussian defined as warup + target_time + saccade_time
        #  + delay, ms
        "mu_dend": [warmup_period - 70, warmup_period + 120,],  # ms, basal synapses
        "mu_apic": {
            "mu_apic": warmup_period + 100,
            "mu_h1": [warmup_period + 100],
            "mu_h2": [warmup_period + 100],
            "mu_h3": [warmup_period + 100],
        },  # ms, apical synapses
        # "mu_oblq":  None,  # warmup_period + dalay,  # ms, oblique synapses
        # -- Standard deviation of the Gaussian
        "sigma_dend": [140, 250],  # ms, width of the Gaussian, Basal Inputs
        # "sigma_oblq": None,  # ms, width of the Gaussian, Oblique Inputs
        "sigma_apic": {
            "sigma_apic": 200,
            "sigma_h1": [200],
            "sigma_h2": [200],
            "sigma_h3": [200],
        },  # ms, width of the Gaussian, Apical Inputs
        # -- Number of pre-synaptic spikes
        "n_dend": [2, 2],  # basal synapses
        # "n_oblq": None,  # oblique synapses
        # apical
        "n_apic": {"n_apic": 1, "n_h1": [1], "n_h2": [1], "n_h3": [1]},
        #
        # Parameters to increase background inputs after target
        "spTimesFun": get_activation_times_from_distribution,
        "args_post_event_bkg_input": dict(
            tstop=cellParameters["tstop"],
            distribution=st.expon,
            rvs_args=dict(loc=0.0, scale=(1 / 4) * 1e3),
        ),
    }

    # - stimulation spike rate
    # ---- excitatory synapses ----
    r_AMPA_dend = 2
    r_NMDA_dend = r_AMPA_dend

    r_AMPA_oblq = 0.5
    r_NMDA_oblq = r_AMPA_oblq

    r_AMPA_apic = 0.5
    r_NMDA_apic = r_AMPA_apic

    # # ---- inhibitory synapses ----
    # r_GABAA_dend = 0
    # r_GABAB_dend = 0
    # r_GABAA_oblq = 0
    # r_GABAB_oblq = 0
    # r_GABAA_apic = 0
    # r_GABAB_apic = 0

    backgrd_inputs_rates = {
        "dend": dict(
            r_NMDA_dend=r_NMDA_dend,
            r_AMPA_dend=r_AMPA_dend,
            # r_GABAA_dend=r_GABAA_dend,
            # r_GABAB_dend=r_GABAB_dend,
        ),
        "oblq": dict(
            r_NMDA_oblq=r_NMDA_oblq,
            r_AMPA_oblq=r_AMPA_oblq,
            # r_GABAA_oblq=r_GABAA_oblq,
            # r_GABAB_oblq=r_GABAB_oblq,
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
    rates_oblq = {
        "r_AMPA": r_AMPA_oblq,
        "r_NMDA": r_NMDA_oblq,
        # "r_GABAA": r_GABAA_oblq,
        # "r_GABAB": r_GABAB_oblq,
    }
    rates_apic = {
        "r_AMPA": r_AMPA_apic,
        "r_NMDA": r_NMDA_apic,
        # "r_GABAA": r_GABAA_apic,
        # "r_GABAB": r_GABAB_apic,
    }

    if RANK == 0:
        file_name = gen_file_name_L5Pop(
            stimulusType,
            **{
                "rates_dend": rates_dend,
                "rates_oblq": rates_oblq,
                "rates_apic": rates_apic,
            }
        )
        data_folder = create_output_folder_L5Pop(
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
    else:  # RANK != 0
        data_folder = None
        file_name = None

    # ---- broadcast data_folder and file_name
    data_folder = COMM.bcast(data_folder, root=0)
    file_name = COMM.bcast(file_name, root=0)

    # ---- run
    for runNumb in trials_ID:

        # delete old sections from NEURON namespace
        LFPy.cell.neuron.h("forall delete_section()")

        # ---- INITIALIZE POPULATION
        population = PopulationL5PC(
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
        population.plotstuff(file_name)
        del population
        gc.collect()


if __name__ == "__main__":
    # load some required neuron-interface files
    h.load_file("stdrun.hoc")
    h.load_file("import3d.hoc")

    main()
