# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:50:03 2022

@author: Beatriz Herrera

This script simulates a population of L5 PCs under random 
synaptic inputs (Poisson pre-synaptic spike trains) 
with or without a synchronized input to all synapses at 
1s after the beginning.

Synapse Configuration: distributed synapses. 
    In this configuration when '"Gaussian_Input": True', all synapses are 
    synchronously activated 'mu' secobds after the beginning of the simulation.
    Pre-synaptic spike times (number of spikes given by 'n') are 
    generated from a generalized Gaussian with mean
    'mu' and standard deviation 'sigma'.

Default values for the parameters are set for the no synchronized input case.

"""
from __future__ import division

import gc
import os
import sys
import LFPy
import numpy as np
from mpi4py import MPI
from neuron import h
from scipy import io
from os.path import join
from PopulationL5PCsHay import PopulationL5PC
from utils import create_output_folder_L5Pop, gen_file_name_L5Pop


# #### MPI ####
# instantize the communication world
COMM = MPI.COMM_WORLD
# get the size of the communication world
SIZE = COMM.Get_size()
# get this particular processes' `rank` ID
RANK = COMM.Get_rank()


""" Functions """


# @profile
def main():

    dt = 2 ** -4  # ms | simulation time step
    SEED = 12  # seed for random number generator
    num_trials = 1  # number of trials to be simulated
    trials_indices = np.arange(1, num_trials + 1)  # trial indices
    POPULATION_SIZE = 100  # number of L5 PCs in the population
    warmup_period = 500  # ms | warm up period, not included in analyses
    sim_length = 10000  # ms | length of the simulation, excluding
    # warm up period

    trial_type = "nosynchinput"
    child_folder = trial_type + "_Feb21_2023"

    if "win32" in sys.platform:  # if running on windows local machine

        trials_ID = trials_indices
        job_0 = True

        main_path_folder = os.path.join(
            r"D:\test_theta_paper\results_L5PCsPopMky", child_folder
        )  # path to folder where results will be saved

    else:  # if running on cluster
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

        main_path_folder = os.path.join("results_L5PCsPopMky", child_folder)  # path to
        # folder where results will be saved

    # ----  Options for cell model and simulation:
    cellParameters = {
        "tstop": warmup_period + sim_length,  # ms | simulation time
        "dt": dt,  # ms | simulation time step
        "mod_Hay_model": True,  # use the modified Hay model
        "save_neurons_data": True,  # whether to save each neuron data
        "plot_synapse_locs": True,  # whether to plot synapse locations
        "show_plot_synapse_locs": False,  # whether to show synapse location plot
        "show_lfp_plot": True,  # whether to show LFP plot
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
        "synapse_type": "distributed",  # Opts: distributed or clustered
        "dend_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
        "oblq_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
        "apic_synp": 4,  # case number within synp configuration options;
        # if 0, no synapses are inserted
    }

    # ----  Options for time-locked inputs:
    inputParameters = {
        "target_times": None,  # target times
        "saccade_times": None,  # saccade times
        "warmup_period": warmup_period,  # length of the warmup period
        #
        # Time-locked inputs definition
        "Gaussian_Input": False,  # time-locked activation of synaptic inputs
        # following Gaussian probability distribution
        #
        # Gaussian Shape Parameters
        "a": 2,
        # center of Gaussian defined as warup + target_time + saccade_time
        #  + delay, ms
        "mu": warmup_period + 1000,  # ms, center of the Gaussian
        "sigma": 200,  # ms, width of the Gaussian
        "n": 1,  # number of pre-synaptic spikes
    }

    # - stimulation spike rate
    # ---- excitatory synapses ----
    r_AMPA_dend = 5
    r_NMDA_dend = r_AMPA_dend

    r_AMPA_oblq = 4
    r_NMDA_oblq = r_AMPA_oblq

    r_AMPA_apic = 1
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

    if RANK == 0:  # if rank 0, create output folder and file name
        file_name = gen_file_name_L5Pop(
            stimulusType,
            **{
                "rates_dend": rates_dend,
                "rates_oblq": rates_oblq,
                "rates_apic": rates_apic,
            }
        )
        data_folder = create_output_folder_L5Pop(
            main_path_folder, POPULATION_SIZE, stimulusType
        )

    else:  # if not rank 0, set data_folder and file_name to None
        data_folder = None
        file_name = None

    data_folder = COMM.bcast(data_folder, root=0)  # broadcast data_folder to all ranks
    file_name = COMM.bcast(file_name, root=0)  # broadcast file_name to all ranks

    # ---- run
    for runNumb in trials_ID:

        # delete old sections from NEURON namespace
        LFPy.cell.neuron.h("forall delete_section()")

        # runNumb = 1

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
