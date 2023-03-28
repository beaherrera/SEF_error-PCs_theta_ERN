# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 17:44:39 2022

@author: Beatriz Herrera
"""

from __future__ import division
import sys
import neo
import LFPy
import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from scipy import io
from mpi4py import MPI
from neuron import h
from os.path import join
from lfpykit import CellGeometry
from elephant import spike_train_generation
from params_synapses import get_synp_params, get_synp_dist_L5PCs
from matplotlib.collections import PolyCollection
from matplotlib.collections import LineCollection

# #### MPI ####
# instantize the communication world
COMM = MPI.COMM_WORLD
# get the size of the communication world
SIZE = COMM.Get_size()
# get this particular processes' `rank` ID
RANK = COMM.Get_rank()


class PopulationL5PC:
    """
    L5-PC population.

    Create a population (unconnected cells) of 'POPULATION_SIZE' L5 PCs
    described by the Hay et al. 2011 model, consisting of LFPy.Cell objects.

    Based on the prototype cell population in example_mpi.py from LFPy package
    examples/.

    """

    def __init__(
        self,
        POPULATION_SIZE,
        cellParameters,
        electrodeParameters,
        stimulusType,
        backgrd_inputs_rates,
        inputParameters,
        runNumb,
        SEED,
        data_folder,
    ):
        """
        Class initialization.

        Parameters
        ----------
        POPULATION_SIZE : int
            number of cells.
        cellParameters : dict
            simulation parameters of the LFPy.Cell object.
        populationParameters : dict
        electrodeParameters : dict
            dictonary containing the coordinates of the sim extracellular
            electrodes, the conductance of the cortex and the method to
            calculate the extracellular potentials using LFPy. Structure for
            LFPy use.
        stimulusType : dict
            dictionary containing the stimulus definition and parameters.
        inputParameters : dict
            dict containing time-locked inputs definition.
        runNumb : int
            run number. Used to simulate different trials.
        SEED: float
            SEED for random numbers generator.
        data_folder : str
            name of the folder to store results.

        Returns
        -------
        None.

        """
        self.POPULATION_SIZE = POPULATION_SIZE
        self.cellParameters = cellParameters
        self.electrodeParameters = electrodeParameters
        self.stimulusType = stimulusType
        self.backgrd_inputs_rates_dend = backgrd_inputs_rates["dend"]
        self.backgrd_inputs_rates_oblq = backgrd_inputs_rates["oblq"]
        self.backgrd_inputs_rates_apic = backgrd_inputs_rates["apic"]
        self.inputParameters = inputParameters
        self.runNumb = runNumb
        self.data_folder = data_folder
        self.SEED = SEED

        # get cell positions, store in self-object
        self.cellPositions = self.setRandCellPositions()
        # self.cellRotations_z = self.drawRandCellRotations()

        self.CELLINDICES = np.arange(self.POPULATION_SIZE)
        self.RANK_CELLINDICES = self.CELLINDICES[self.CELLINDICES % SIZE == RANK]

        # container for single-cell output generated on this RANK
        self.results = dict((i, {}) for i in self.RANK_CELLINDICES)

    def run(self, file_name):
        """Execute the proper simulation and collect simulation results.

        Parameters
        ----------
        file_name : _type_
            path to saving directory.
        """
        # produce simulation results
        self.distribute_cellsims(file_name)

        RANK_CELLINDICES = []
        for i in range(SIZE):
            RANK_CELLINDICES += [self.CELLINDICES[self.CELLINDICES % SIZE == i]]

        if self.POPULATION_SIZE <= 10:
            # gather data on this RANK
            if RANK_CELLINDICES[RANK].size > 0:
                cell_geo_temp = dict(
                    (i, {}) for i in range(RANK_CELLINDICES[RANK].size)
                )
                for i, cellindex in enumerate(RANK_CELLINDICES[RANK]):
                    if i == 0:
                        lfp_temp = np.zeros(
                            [RANK_CELLINDICES[RANK].size]
                            + list(self.results[cellindex]["LFP"].shape),
                            dtype=np.float32,
                        )
                        somav_temp = np.zeros(
                            [RANK_CELLINDICES[RANK].size]
                            + list(self.results[cellindex]["somav"].shape),
                            dtype=np.float32,
                        )
                        v_mbp_temp = np.zeros(
                            [RANK_CELLINDICES[RANK].size]
                            + list(self.results[cellindex]["v_mbp"].shape),
                            dtype=np.float32,
                        )

                    lfp_temp[i,] = self.results[cellindex]["LFP"]  #
                    # Ncells x Nelectrodes x Ntimes
                    somav_temp[i,] = self.results[cellindex]["somav"]  #
                    # Ncells xNtimes
                    v_mbp_temp[i,] = self.results[cellindex]["v_mbp"]  #
                    # Ncells xNtimes
                    cell_geo_temp[i] = self.results[cellindex]["cell_geo"]
                    # Ncells -> one object per cell

            if RANK == 0:
                # container of all output
                self.lfp = np.zeros(
                    [self.POPULATION_SIZE] + list(self.results[cellindex]["LFP"].shape),
                    dtype=np.float32,
                )
                self.somav = np.zeros(
                    [self.POPULATION_SIZE]
                    + list(self.results[cellindex]["somav"].shape),
                    dtype=np.float32,
                )
                self.v_mbp = np.zeros(
                    [self.POPULATION_SIZE]
                    + list(self.results[cellindex]["v_mbp"].shape),
                    dtype=np.float32,
                )
                self.cell_geo = dict((i, {}) for i in range(self.POPULATION_SIZE))

                # fill in values from this RANK
                if RANK_CELLINDICES[0].size > 0:
                    for j, k in enumerate(RANK_CELLINDICES[0]):
                        self.lfp[k,] = lfp_temp[
                            j,
                        ]
                        self.somav[k,] = somav_temp[
                            j,
                        ]
                        self.v_mbp[k,] = v_mbp_temp[
                            j,
                        ]
                        self.cell_geo[k] = cell_geo_temp[j]

                # iterate over all other RANKs
                for i in range(1, len(RANK_CELLINDICES)):
                    if RANK_CELLINDICES[i].size > 0:
                        # receive on RANK 0 from all other RANK
                        lfp_temp = np.zeros(
                            [RANK_CELLINDICES[i].size]
                            + list(self.results[cellindex]["LFP"].shape),
                            dtype=np.float32,
                        )
                        somav_temp = np.zeros(
                            [RANK_CELLINDICES[i].size]
                            + list(self.results[cellindex]["somav"].shape),
                            dtype=np.float32,
                        )
                        v_mbp_temp = np.zeros(
                            [RANK_CELLINDICES[i].size]
                            + list(self.results[cellindex]["v_mbp"].shape),
                            dtype=np.float32,
                        )

                        COMM.Recv([lfp_temp, MPI.FLOAT], source=i, tag=13)
                        COMM.Recv([somav_temp, MPI.FLOAT], source=i, tag=14)
                        COMM.Recv([v_mbp_temp, MPI.FLOAT], source=i, tag=15)
                        cell_geo_temp = COMM.recv(source=i, tag=16)

                        # fill in values
                        for j, k in enumerate(RANK_CELLINDICES[i]):
                            self.lfp[k,] = lfp_temp[
                                j,
                            ]
                            self.somav[k,] = somav_temp[
                                j,
                            ]
                            self.v_mbp[k,] = v_mbp_temp[
                                j,
                            ]
                            self.cell_geo[k] = cell_geo_temp[j]
                self.lfp = np.array(self.lfp).sum(axis=0)  # summing for the cells
                # -> Nelectrodes x Ntimes
            else:
                self.lfp = None
                self.somav = None
                self.v_mbp = None
                self.cell_geo = None

                if RANK_CELLINDICES[RANK].size > 0:
                    # send to RANK 0
                    COMM.Send([lfp_temp, MPI.FLOAT], dest=0, tag=13)
                    COMM.Send([somav_temp, MPI.FLOAT], dest=0, tag=14)
                    COMM.Send([v_mbp_temp, MPI.FLOAT], dest=0, tag=15)
                    COMM.send(cell_geo_temp, dest=0, tag=16)

        else:
            # gather data on this RANK
            if RANK_CELLINDICES[RANK].size > 0:
                for i, cellindex in enumerate(RANK_CELLINDICES[RANK]):
                    if i == 0:
                        lfp_temp = np.zeros(
                            [RANK_CELLINDICES[RANK].size]
                            + list(self.results[cellindex]["LFP"].shape),
                            dtype=np.float32,
                        )
                    lfp_temp[i,] = self.results[cellindex]["LFP"]  #
                    # Ncells x Nelectrodes x Ntimes

            if RANK == 0:
                # container of all output
                self.lfp = np.zeros(
                    [self.POPULATION_SIZE] + list(self.results[cellindex]["LFP"].shape),
                    dtype=np.float32,
                )

                # fill in values from this RANK
                if RANK_CELLINDICES[0].size > 0:
                    for j, k in enumerate(RANK_CELLINDICES[0]):
                        self.lfp[k,] = lfp_temp[
                            j,
                        ]

                # iterate over all other RANKs
                for i in range(1, len(RANK_CELLINDICES)):
                    if RANK_CELLINDICES[i].size > 0:
                        # receive on RANK 0 from all other RANK
                        lfp_temp = np.zeros(
                            [RANK_CELLINDICES[i].size]
                            + list(self.results[cellindex]["LFP"].shape),
                            dtype=np.float32,
                        )
                        COMM.Recv([lfp_temp, MPI.FLOAT], source=i, tag=13)

                        # fill in values
                        for j, k in enumerate(RANK_CELLINDICES[i]):
                            self.lfp[k,] = lfp_temp[
                                j,
                            ]

                self.lfp = np.array(self.lfp).sum(axis=0)  # summing for the cells
                # -> Nelectrodes x Ntimes
            else:
                self.lfp = None

                if RANK_CELLINDICES[RANK].size > 0:
                    # send to RANK 0
                    COMM.Send([lfp_temp, MPI.FLOAT], dest=0, tag=13)

        COMM.Barrier()

    def distribute_cellsims(self, file_name):
        """Will run cell simulations."""
        # start unique cell simulation on every RANK,
        # and store the electrode and cell objects in dicts indexed by
        # cellindex

        for cellindex in self.RANK_CELLINDICES:
            self.cellsim(cellindex, file_name)

        COMM.Barrier()

    def cellsim(self, cellindex, file_name):
        """
        Run cell and LFP simulation procedure.

        Parameters
        ----------
        cellindex : int
            Index of the cell in population being simulated.

        Returns
        -------
        None.

        """
        # Initialize cell instance
        cell = self.make_cellHayModelL5PC()

        # set the position of midpoint in soma
        cell.set_pos(
            x=self.cellPositions[cellindex, 0],
            y=self.cellPositions[cellindex, 1],
            z=self.cellPositions[cellindex, 2],
        )

        # attach synapse with parameters and set spikes time
        # cell, synapse, isyn = make_backgroundNoise(cell=cell)
        if self.stimulusType["synapse_type"] == "clustered":
            cell = self.make_clustered_synapses(cell, cellindex)
        elif self.stimulusType["synapse_type"] == "distributed":
            cell = self.make_distributed_synapses(cell, cellindex)

        # create extracellular electrode object
        electrode = LFPy.RecExtElectrode(cell, **self.electrodeParameters)

        # Parameters for the cell.simulate() call,
        # recording membrane- and syn.-currents
        simulationParameters = {
            "probes": [electrode],
            # 'rec_imem' : True,  # Record Membrane currents during simulation
            "rec_vmem": True,  # record membrane voltage
        }

        # perform NEURON simulation, results saved as attributes in cell
        cell.simulate(**simulationParameters)

        # store synapse type
        syntype = []
        for i in range(len(cell.synapses)):
            syntype.append(str(cell.synapses[i].syntype))

        # store cells geo for the plots
        zips = []
        for x, z in cell.get_idx_polygons():
            zips.append(list(zip(x, z)))

        cell_geo = CellGeometry(x=cell.x, y=cell.y, z=cell.z, d=cell.d)
        cell_geo.zips = zips

        # extract spike times: soma and dendrites
        somav = cell.somav[np.where(cell.tvec >= self.inputParameters["warmup_period"])]
        dendv = cell.vmem[
            616, np.where(cell.tvec >= self.inputParameters["warmup_period"]),
        ]

        signal_somav = neo.core.AnalogSignal(
            somav, units="mV", sampling_rate=(1 / (cell.dt * 1e-3)) * pq.Hz
        )
        signal_dendv = neo.core.AnalogSignal(
            dendv.T, units="mV", sampling_rate=(1 / (cell.dt * 1e-3)) * pq.Hz
        )

        soma_spkTimes = spike_train_generation.peak_detection(
            signal_somav, threshold=0.0 * pq.mV, sign="above", as_array=False
        )

        dend_spkTimes = spike_train_generation.peak_detection(
            signal_dendv, threshold=-25.0 * pq.mV, sign="above", as_array=False
        )

        # save individual neuron data
        if self.cellParameters["save_neurons_data"]:
            saveData = {
                "cell_geo": cell_geo,  # geometry L5 PCs
                # [mV] somatic membrane potatential
                "Vs": np.array(cell.somav),
                # [mV] distal dendrites memberane potential
                "v_mbp": np.array(cell.vmem[616, :]),
                # [ms] time of presynaptic spikes Ncells x Nsynp x Ntskp
                "pre_synSpkTimes": np.array(cell._sptimeslist, dtype="object"),
                # type of synapse associated with pre_synSpkTimes, same order
                "synType": np.array(syntype, dtype="object"),
                "syn_indx": cell.synidx,  # location of synapses
                "soma_spkTimes": soma_spkTimes,  # [s] times of somatic APs
                "dend_spkTimes": dend_spkTimes,  # [s]
                # [mV] lfp produced by the neuron
                "LFP_neuron": electrode.data,
            }

            # mat files
            io.savemat(
                join(
                    self.data_folder,
                    "NeuronsData_r"
                    + str(self.runNumb)
                    + "_n#"
                    + str(cellindex)
                    + "_"
                    + file_name
                    + ".mat",
                ),
                saveData,
            )

        # store primary results from simulation
        if self.POPULATION_SIZE <= 10:
            self.results[cellindex]["LFP"] = electrode.data
            self.results[cellindex]["cell_geo"] = zips
            self.results[cellindex]["somav"] = cell.somav
            self.results[cellindex]["v_mbp"] = cell.vmem[616, :]
        else:
            self.results[cellindex]["LFP"] = electrode.data

    def setRandCellPositions(self):
        """
        set random cell positions within a cylinder constraints.

        Returns
        -------
        ndarray
            cells location within the cylinder.

        """
        if RANK == 0:
            data = np.load("cellsPops_SEF.npz")
            cellPositions = data["cellPositions_L5"]
        else:
            cellPositions = None

        return COMM.bcast(cellPositions, root=0)

    def make_cellHayModelL5PC(self):
        """
        Create LFPy cell object using Hay et al. 2011 model.

        Parameters
        ----------
        None.

        Returns
        -------
        cell : TYPE
            DESCRIPTION.

        """
        holding_potential = -80  # [mV] resting membrane potential

        # define cell parameters used as input to cell-class
        model_pth = join("cell_models", "HayModel")

        if self.cellParameters["mod_Hay_model"]:
            custom_codes_file = "custom_codes_2.hoc"
            from modhayetalmodel_active_declarations import active_declarations
        else:
            custom_codes_file = "custom_codes.hoc"
            from hayetalmodel_active_declarations import active_declarations

        cell_Parameters = {
            "morphology": join(model_pth, "morphologies", "cell1.hoc"),
            "v_init": holding_potential,  # initial crossmembrane potential
            "passive": False,  # switch on passive mechs
            "nsegs_method": "lambda_f",  # method for setting number of
            # segments
            "lambda_f": 100,  # segments are isopotential at this frequency
            "dt": self.cellParameters["dt"],  # dt of LFP and NEURON simulation
            "tstart": -250,  # start time, recorders start at t=0
            "tstop": self.cellParameters["tstop"],  # stop time,
            # end of the simulation
            "custom_code": [join(model_pth, "morphologies", custom_codes_file)],
            "custom_fun": [active_declarations],  # will execute this function
            "custom_fun_args": [{}],
        }

        cell = LFPy.Cell(**cell_Parameters)
        cell.set_rotation(z=np.pi)

        return cell

    def make_distributed_synapses(self, cell, cell_index):
        """
        Create randomly distributed synapses.

        Parameters
        ----------
        cell : TYPE
            DESCRIPTION.
        cell_index : TYPE
            DESCRIPTION.

        Returns
        -------
        cell : TYPE
            DESCRIPTION.

        """
        from make_synapses import insert_dist_synapses

        # get synapse parameters
        synp_models_params = get_synp_params()
        syn_models_params_AMPA = synp_models_params["syn_models_params_AMPA"]
        syn_models_params_NMDA = synp_models_params["syn_models_params_NMDA"]
        syn_models_params_GABAA = synp_models_params["syn_models_params_GABAA"]
        syn_models_params_GABAB = synp_models_params["syn_models_params_GABAB"]

        # get synapse distribution
        syn_params = get_synp_dist_L5PCs(
            cell.tstop,
            **self.backgrd_inputs_rates_dend,
            **self.backgrd_inputs_rates_oblq,
            **self.backgrd_inputs_rates_apic
        )
        insert_synapses_NMDA_args = syn_params["insert_synapses_NMDA_args"]
        insert_synapses_AMPA_args = syn_params["insert_synapses_AMPA_args"]
        insert_synapses_GABAA_args = syn_params["insert_synapses_GABAA_args"]
        insert_synapses_GABAB_args = syn_params["insert_synapses_GABAB_args"]

        if self.cellParameters["plot_synapse_locs"]:
            idx_synp_NMDA = np.array(int(0))
            idx_synp_AMPA = np.array(int(0))
            idx_synp_GABAA = np.array(int(0))
            idx_synp_GABAB = np.array(int(0))

        if self.stimulusType["dend_synp"] > 0:
            # ---- dend inputs

            # condition for NMDA synapses
            NMDA_dend = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if NMDA_dend:
                syn_params_NMDA_dend = {
                    "section": "dend",
                    "n": insert_synapses_NMDA_args["n_dend"],
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_dend"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_NMDA,
                    syn_params_NMDA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    idx_synp_NMDA = idx_synp_NMDA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_dend = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if AMPA_dend:
                syn_params_AMPA_dend = {
                    "section": "dend",
                    "n": insert_synapses_AMPA_args["n_dend"],
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_dend"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_AMPA,
                    syn_params_AMPA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    idx_synp_AMPA = idx_synp_AMPA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_dend = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if GABAA_dend:
                syn_params_GABAA_dend = {
                    "section": "dend",
                    "n": insert_synapses_GABAA_args["n_dend"],
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_dend"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAA,
                    syn_params_GABAA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    idx_synp_GABAA = idx_synp_GABAA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_dend = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if GABAB_dend:
                syn_params_GABAB_dend = {
                    "section": "dend",
                    "n": insert_synapses_GABAB_args["n_dend"],
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_dend"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAB,
                    syn_params_GABAB_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    idx_synp_GABAB = idx_synp_GABAB[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAB = []

        if self.stimulusType["oblq_synp"] > 0:
            # ---- oblq inputs

            # condition for NMDA synapses
            NMDA_oblq = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if NMDA_oblq:
                syn_params_NMDA_oblq = {
                    "section": "oblq",
                    "n": insert_synapses_NMDA_args["n_dend"],
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_oblq"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_NMDA,
                    syn_params_NMDA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    idx_synp_NMDA = idx_synp_NMDA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_oblq = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if AMPA_oblq:
                syn_params_AMPA_oblq = {
                    "section": "oblq",
                    "n": insert_synapses_AMPA_args["n_dend"],
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_oblq"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_AMPA,
                    syn_params_AMPA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    idx_synp_AMPA = idx_synp_AMPA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_oblq = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if GABAA_oblq:
                syn_params_GABAA_oblq = {
                    "section": "oblq",
                    "n": insert_synapses_GABAA_args["n_dend"],
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_oblq"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAA,
                    syn_params_GABAA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    idx_synp_GABAA = idx_synp_GABAA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_oblq = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if GABAB_oblq:
                syn_params_GABAB_oblq = {
                    "section": "oblq",
                    "n": insert_synapses_GABAB_args["n_dend"],
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_oblq"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAB,
                    syn_params_GABAB_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    idx_synp_GABAB = idx_synp_GABAB[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAB = []

        if self.stimulusType["apic_synp"] > 0:
            # ---- apic inputs

            # condition for NMDA synapses
            NMDA_apic = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if NMDA_apic:
                syn_params_NMDA_apic = {
                    "section": "apic",
                    "n": insert_synapses_NMDA_args["n_apic"],
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_apic"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_NMDA,
                    syn_params_NMDA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    if idx_synp_NMDA[0] == 0:
                        idx_synp_NMDA = idx_synp_NMDA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_NMDA.size == 1:
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_apic = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if AMPA_apic:
                syn_params_AMPA_apic = {
                    "section": "apic",
                    "n": insert_synapses_AMPA_args["n_apic"],
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_apic"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_AMPA,
                    syn_params_AMPA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    if idx_synp_AMPA[0] == 0:
                        idx_synp_AMPA = idx_synp_AMPA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_AMPA.size == 1:
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_apic = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if GABAA_apic:
                syn_params_GABAA_apic = {
                    "section": "apic",
                    "n": insert_synapses_GABAA_args["n_apic"],
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_apic"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAA,
                    syn_params_GABAA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    if idx_synp_GABAA[0] == 0:
                        idx_synp_GABAA = idx_synp_GABAA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_GABAA.size == 1:
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_apic = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if GABAB_apic:
                syn_params_GABAB_apic = {
                    "section": "apic",
                    "n": insert_synapses_GABAB_args["n_apic"],
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_apic"],
                }
                cell, idx_plot = insert_dist_synapses(
                    cell,
                    cell_index,
                    syn_models_params_GABAB,
                    syn_params_GABAB_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    if idx_synp_GABAB[0] == 0:
                        idx_synp_GABAB = idx_synp_GABAB[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_GABAB.size == 1:
                idx_synp_GABAB = []

        # plotting synaptic locs
        if self.cellParameters["plot_synapse_locs"]:
            # figsize=[0.73,1.38],
            fig = plt.figure(figsize=[1.88, 4.2], dpi=600)
            plt.axes([0.175, 0.0, 0.65, 1], aspect="equal")
            plt.axis("off")
            zips = []
            for x, z in cell.get_idx_polygons():
                zips.append(list(zip(x, z)))
            linecol = LineCollection(
                zips, edgecolor="none", facecolor="k", rasterized=False,
            )
            ax = plt.gca()
            ax.add_collection(linecol)
            # adding dots at synaptic locs
            ax.plot(
                np.median(cell.x[idx_synp_NMDA, :], axis=1),
                np.median(cell.z[idx_synp_NMDA, :], axis=1),
                "o",
                markeredgecolor="none",
                markerfacecolor="b",
                markersize=1.5,
                label="NMDA",
            )
            ax.plot(
                np.median(cell.x[idx_synp_AMPA, :], axis=1),
                np.median(cell.z[idx_synp_AMPA, :], axis=1),
                "o",
                markeredgecolor="none",
                markerfacecolor="r",
                markersize=1.5,
                label="AMPA",
            )
            ax.scatter(
                np.median(cell.x[idx_synp_GABAA, :], axis=1),
                np.median(cell.z[idx_synp_GABAA, :], axis=1),
                edgecolor="g",
                color="none",
                s=5,
                label="GABAA",
            )
            ax.scatter(
                np.median(cell.x[idx_synp_GABAB, :], axis=1),
                np.median(cell.z[idx_synp_GABAB, :], axis=1),
                edgecolor="k",
                color="none",
                s=5,
                label="GABAB",
            )
            # ax.legend(bbox_to_anchor=(0.3, -0.07, 0.1, 0.5),
            #           loc="lower center",
            #           ncol=2)

            if ("win32" in sys.platform) and self.cellParameters[
                "show_plot_synapse_locs"
            ]:
                plt.show()
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_n#"
                        + str(cell_index)
                        + "_SynpLocs.png",
                    )
                )
            else:
                plt.ioff()
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_n#"
                        + str(cell_index)
                        + "_SynpLocs.png",
                    )
                )
                fig.clear()
                plt.close(fig)

        return cell

    def make_clustered_synapses(self, cell, cell_index):
        """
        Create randomly distributed clustered synapses.

        Parameters
        ----------
        cell : TYPE
            DESCRIPTION.
        cell_index : TYPE
            DESCRIPTION.

        Returns
        -------
        cell : TYPE
            DESCRIPTION.

        """
        from make_synapses import fill_clustered_synapses_vectors
        from make_synapses import insert_clustered_synapses_l5
        from params_synapses import get_cluster_parameters

        cluster_parameters = get_cluster_parameters()

        # get sysnapse parameters
        synp_models_params = get_synp_params()
        syn_models_params_AMPA = synp_models_params["syn_models_params_AMPA"]
        syn_models_params_NMDA = synp_models_params["syn_models_params_NMDA"]
        syn_models_params_GABAA = synp_models_params["syn_models_params_GABAA"]
        syn_models_params_GABAB = synp_models_params["syn_models_params_GABAB"]

        # get synapse distribution
        syn_params = get_synp_dist_L5PCs(
            cell.tstop,
            **self.backgrd_inputs_rates_dend,
            **self.backgrd_inputs_rates_oblq,
            **self.backgrd_inputs_rates_apic
        )
        insert_synapses_NMDA_args = syn_params["insert_synapses_NMDA_args"]
        insert_synapses_AMPA_args = syn_params["insert_synapses_AMPA_args"]
        insert_synapses_GABAA_args = syn_params["insert_synapses_GABAA_args"]
        insert_synapses_GABAB_args = syn_params["insert_synapses_GABAB_args"]

        rd = h.Random(self.SEED * (cell_index + 1))  # generate different
        # synapse locations in each run

        if self.cellParameters["plot_synapse_locs"]:
            idx_synp_NMDA = np.array(int(0))
            idx_synp_AMPA = np.array(int(0))
            idx_synp_GABAA = np.array(int(0))
            idx_synp_GABAB = np.array(int(0))

        if self.stimulusType["dend_synp"] > 0:
            # ---- dend inputs
            # condition for NMDA synapses
            NMDA_dend = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if NMDA_dend:
                sec_name_NMDA, seg_x_NMDA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_NMDA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "dend",
                )

                syn_params_NMDA_dend = {
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_dend"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "dend",
                    sec_name_NMDA,
                    seg_x_NMDA,
                    syn_models_params_NMDA,
                    syn_params_NMDA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    idx_synp_NMDA = idx_synp_NMDA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_dend = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if AMPA_dend:
                sec_name_AMPA, seg_x_AMPA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_AMPA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "dend",
                )

                syn_params_AMPA_dend = {
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_dend"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "dend",
                    sec_name_AMPA,
                    seg_x_AMPA,
                    syn_models_params_AMPA,
                    syn_params_AMPA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    idx_synp_AMPA = idx_synp_AMPA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_dend = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if GABAA_dend:
                sec_name_GABAA, seg_x_GABAA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "dend",
                )

                syn_params_GABAA_dend = {
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_dend"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "dend",
                    sec_name_GABAA,
                    seg_x_GABAA,
                    syn_models_params_GABAA,
                    syn_params_GABAA_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    idx_synp_GABAA = idx_synp_GABAA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_dend = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["dend_synp"],
                )
            )
            if GABAB_dend:
                sec_name_GABAB, seg_x_GABAB = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAB_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "dend",
                )

                syn_params_GABAB_dend = {
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_dend"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "dend",
                    sec_name_GABAB,
                    seg_x_GABAB,
                    syn_models_params_GABAB,
                    syn_params_GABAB_dend,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    idx_synp_GABAB = idx_synp_GABAB[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAB = []

        if self.stimulusType["oblq_synp"] > 0:
            # ---- oblq inputs
            # condition for NMDA synapses
            NMDA_oblq = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if NMDA_oblq:
                sec_name_NMDA, seg_x_NMDA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_NMDA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "oblq",
                )

                syn_params_NMDA_oblq = {
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_oblq"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "oblq",
                    sec_name_NMDA,
                    seg_x_NMDA,
                    syn_models_params_NMDA,
                    syn_params_NMDA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    idx_synp_NMDA = idx_synp_NMDA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_oblq = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if AMPA_oblq:
                sec_name_AMPA, seg_x_AMPA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_AMPA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "oblq",
                )

                syn_params_AMPA_oblq = {
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_oblq"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "oblq",
                    sec_name_AMPA,
                    seg_x_AMPA,
                    syn_models_params_AMPA,
                    syn_params_AMPA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    idx_synp_AMPA = idx_synp_AMPA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_oblq = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if GABAA_oblq:
                sec_name_GABAA, seg_x_GABAA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAA_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "oblq",
                )

                syn_params_GABAA_oblq = {
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_oblq"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "oblq",
                    sec_name_GABAA,
                    seg_x_GABAA,
                    syn_models_params_GABAA,
                    syn_params_GABAA_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    idx_synp_GABAA = idx_synp_GABAA[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_oblq = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["oblq_synp"],
                )
            )
            if GABAB_oblq:
                sec_name_GABAB, seg_x_GABAB = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAB_args["n_dend"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "oblq",
                )

                syn_params_GABAB_oblq = {
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_oblq"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "oblq",
                    sec_name_GABAB,
                    seg_x_GABAB,
                    syn_models_params_GABAB,
                    syn_params_GABAB_oblq,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    idx_synp_GABAB = idx_synp_GABAB[1:]
            elif (
                self.cellParameters["plot_synapse_locs"]
                and self.stimulusType["apic_synp"] == 0
            ):
                idx_synp_GABAB = []

        if self.stimulusType["apic_synp"] > 0:
            # ---- apic inputs

            # condition for NMDA synapses
            NMDA_apic = any(
                np.isin(
                    np.concatenate((np.arange(1, 8), 11), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if NMDA_apic:
                sec_name_NMDA, seg_x_NMDA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_NMDA_args["n_apic"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "apic",
                )

                syn_params_NMDA_apic = {
                    "spTimesFun": insert_synapses_NMDA_args["spTimesFun"],
                    "args": insert_synapses_NMDA_args["args_apic"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "apic",
                    sec_name_NMDA,
                    seg_x_NMDA,
                    syn_models_params_NMDA,
                    syn_params_NMDA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_NMDA = np.concatenate((idx_synp_NMDA, idx_plot), axis=None)
                    if idx_synp_NMDA[0] == 0:
                        idx_synp_NMDA = idx_synp_NMDA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_NMDA.size == 1:
                idx_synp_NMDA = []

            # condition for AMPA synapses
            AMPA_apic = any(
                np.isin(
                    np.concatenate((np.arange(1, 5), 8, 9, 10, 12), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if AMPA_apic:
                sec_name_AMPA, seg_x_AMPA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_AMPA_args["n_apic"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "apic",
                )

                syn_params_AMPA_apic = {
                    "spTimesFun": insert_synapses_AMPA_args["spTimesFun"],
                    "args": insert_synapses_AMPA_args["args_apic"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "apic",
                    sec_name_AMPA,
                    seg_x_AMPA,
                    syn_models_params_AMPA,
                    syn_params_AMPA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_AMPA = np.concatenate((idx_synp_AMPA, idx_plot), axis=None)
                    if idx_synp_AMPA[0] == 0:
                        idx_synp_AMPA = idx_synp_AMPA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_AMPA.size == 1:
                idx_synp_AMPA = []

            # condition for GABAA synapses
            GABAA_apic = any(
                np.isin(
                    np.concatenate((1, 2, 5, 6, 8, 9), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if GABAA_apic:
                sec_name_GABAA, seg_x_GABAA = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAA_args["n_apic"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "apic",
                )

                syn_params_GABAA_apic = {
                    "spTimesFun": insert_synapses_GABAA_args["spTimesFun"],
                    "args": insert_synapses_GABAA_args["args_apic"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "apic",
                    sec_name_GABAA,
                    seg_x_GABAA,
                    syn_models_params_GABAA,
                    syn_params_GABAA_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAA = np.concatenate(
                        (idx_synp_GABAA, idx_plot), axis=None
                    )
                    if idx_synp_GABAA[0] == 0:
                        idx_synp_GABAA = idx_synp_GABAA[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_GABAA.size == 1:
                idx_synp_GABAA = []

            # condition for GABAB synapses
            GABAB_apic = any(
                np.isin(
                    np.concatenate((1, 3, 5, 7, 8, 10), axis=None),
                    self.stimulusType["apic_synp"],
                )
            )
            if GABAB_apic:
                sec_name_GABAB, seg_x_GABAB = fill_clustered_synapses_vectors(
                    cell,
                    round(
                        insert_synapses_GABAB_args["n_apic"]
                        / cluster_parameters["CLUSTER_SIZE"]
                    ),
                    rd,
                    "apic",
                )

                syn_params_GABAB_apic = {
                    "spTimesFun": insert_synapses_GABAB_args["spTimesFun"],
                    "args": insert_synapses_GABAB_args["args_apic"],
                }
                cell, idx_plot = insert_clustered_synapses_l5(
                    cell,
                    cell_index,
                    "apic",
                    sec_name_GABAB,
                    seg_x_GABAB,
                    syn_models_params_GABAB,
                    syn_params_GABAB_apic,
                    self.inputParameters,
                    self.runNumb,
                    self.SEED,
                    self.stimulusType["background_inputs"],
                )

                if self.cellParameters["plot_synapse_locs"]:
                    idx_synp_GABAB = np.concatenate(
                        (idx_synp_GABAB, idx_plot), axis=None
                    )
                    if idx_synp_GABAB[0] == 0:
                        idx_synp_GABAB = idx_synp_GABAB[1:]
            elif self.cellParameters["plot_synapse_locs"] and idx_synp_GABAB.size == 1:
                idx_synp_GABAB = []

        # plotting synaptic locs
        if self.cellParameters["plot_synapse_locs"]:
            # figsize=[0.73,1.38],
            fig = plt.figure(figsize=[1.88, 4.2], dpi=600)
            plt.axes([0.175, 0.0, 0.65, 1], aspect="equal")
            plt.axis("off")
            zips = []
            for x, z in cell.get_idx_polygons():
                zips.append(list(zip(x, z)))
            linecol = LineCollection(
                zips, edgecolor="none", facecolor="k", rasterized=False,
            )
            ax = plt.gca()
            ax.add_collection(linecol)
            # adding dots at synaptic locs
            ax.plot(
                np.median(cell.x[idx_synp_NMDA, :], axis=1),
                np.median(cell.z[idx_synp_NMDA, :], axis=1),
                "o",
                markeredgecolor="none",
                markerfacecolor="b",
                markersize=1.5,
                label="NMDA",
            )
            ax.plot(
                np.median(cell.x[idx_synp_AMPA, :], axis=1),
                np.median(cell.z[idx_synp_AMPA, :], axis=1),
                "o",
                markeredgecolor="none",
                markerfacecolor="r",
                markersize=1.5,
                label="AMPA",
            )
            ax.scatter(
                np.median(cell.x[idx_synp_GABAA, :], axis=1),
                np.median(cell.z[idx_synp_GABAA, :], axis=1),
                edgecolor="g",
                color="none",
                s=5,
                label="GABAA",
            )
            ax.scatter(
                np.median(cell.x[idx_synp_GABAB, :], axis=1),
                np.median(cell.z[idx_synp_GABAB, :], axis=1),
                edgecolor="k",
                color="none",
                s=5,
                label="GABAB",
            )

            if ("win32" in sys.platform) and self.cellParameters[
                "show_plot_synapse_locs"
            ]:
                plt.show()
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_n#"
                        + str(cell_index)
                        + "_ClustSynpLocs.png",
                    )
                )
            else:
                plt.ioff()
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_n#"
                        + str(cell_index)
                        + "_ClustSynpLocs.png",
                    )
                )
                fig.clear()
                plt.close(fig)

        return cell

    def save_simData(self, file_name):
        """Save simulations data."""
        if RANK == 0:
            saveData = {
                "ze": self.electrodeParameters["z"],  # [um] electrodes position
                "LFP": self.lfp,  # [mV] LFP values
            }

            # mat files
            io.savemat(
                join(
                    self.data_folder,
                    "SimData_r"
                    + str(self.runNumb)
                    + "_PS"
                    + str(self.POPULATION_SIZE)
                    + "_"
                    + file_name
                    + ".mat",
                ),
                saveData,
            )

    def plotstuff(self, file_name):
        """Plot LFPs and membrane voltage traces."""
        if RANK == 0:
            color = iter(plt.cm.rainbow(np.linspace(0, 1, self.POPULATION_SIZE)))

            fig = plt.figure(figsize=(12, 8), dpi=600)  #

            if self.POPULATION_SIZE <= 10:
                ax = fig.add_axes(
                    [0.05, 0.0, 0.45, 1.0],
                    aspect="equal",
                    frameon=False,
                    xticks=[],
                    xticklabels=[],
                    yticks=[],
                    yticklabels=[],
                )
                for cellindex in range(self.POPULATION_SIZE):
                    c = next(color)
                    polycol = PolyCollection(
                        self.cell_geo[cellindex],
                        edgecolors="none",
                        facecolors=c,
                        zorder=self.cellPositions[cellindex, 1],
                    )

                    ax.add_collection(polycol)

                ax.plot(
                    self.electrodeParameters["x"],
                    self.electrodeParameters["z"],
                    marker="o",
                    color="k",
                    clip_on=False,
                    zorder=0,
                )

                ax2 = fig.add_axes([0.5, 0.7, 0.3, 0.25])
                ax1 = fig.add_axes([0.5, 0.385, 0.3, 0.25])
                tvec = (
                    np.arange(np.size(self.somav, axis=1)) * self.cellParameters["dt"]
                )
                for i, key in enumerate(np.arange(np.size(self.somav, axis=0))):
                    ax1.plot(tvec, self.somav[i,], label="cell %i" % key)
                    ax2.plot(tvec, self.v_mbp[i,], label="cell %i" % key)
                ax1.legend()
                # plt.ylim(-100,50)
                ax1.set_ylabel("$V_{soma}$ (mV)")
                ax1.set_title("somatic potentials")
                ax2.set_title("dendritic potentials")
                ax1.set_ylim(-100, 50)
                ax2.set_ylim(-100, 50)
                ax1.set_xlim(0, self.cellParameters["tstop"])
                ax2.set_xlim(0, self.cellParameters["tstop"])

            ax = fig.add_axes([0.5, 0.075, 0.3, 0.25])
            cax = fig.add_axes([0.825, 0.075, 0.02, 0.25])
            # cax = fig.add_axes([0.91, 0.55, 0.02, 0.40])
            tvec = np.arange(np.size(self.lfp, axis=1)) * self.cellParameters["dt"]
            im = ax.pcolormesh(
                tvec,
                self.electrodeParameters["z"],
                self.lfp,
                cmap="PRGn",
                vmin=-self.lfp.std() * 3,
                vmax=self.lfp.std() * 3,
                shading="auto",
            )
            ax.axis(ax.axis("tight"))
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label("LFP (mV)")
            ax.set_title("superimposed LFP")
            ax.set_xlabel("time (ms)")
            ax.set_ylabel(r"$z$ ($\mu$m)")
            ax.set_xlim(0, self.cellParameters["tstop"])

            if ("win32" in sys.platform) and self.cellParameters["show_lfp_plot"]:
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_PS"
                        + str(self.POPULATION_SIZE)
                        + "+"
                        + file_name
                        + ".png",
                    )
                )
                plt.show()
            else:
                plt.ioff()
                fig.savefig(
                    join(
                        self.data_folder,
                        "Fig_run#"
                        + str(self.runNumb)
                        + "_PS"
                        + str(self.POPULATION_SIZE)
                        + "+"
                        + file_name
                        + ".png",
                    )
                )
                fig.clear()
                plt.close(fig)
