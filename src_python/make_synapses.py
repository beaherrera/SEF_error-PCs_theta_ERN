# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 12:21:11 2022

@author: Beatriz Herrera

Synapses Configuration.

"""
import LFPy
import utils
import neuron
import itertools
import numpy as np
import scipy.stats as st
from neuron import h


def _basal_random_synapse(HCell, synaptic_loc):
    """
    Return a random location in the basal tree of this cell.

    Codes taken from Eyal et al. 2018.

    Parameters
    ----------
    synaptic_loc : array
        synaptic locations.

    Returns
    -------
    sec : list
        List with the name of the sections where synapses will be inserted.
    x : list
        List of seg.x values associated with sec.

    """
    len0 = 0
    len1 = 0
    for sec in HCell.basal:
        len1 += sec.L
        if len1 >= synaptic_loc:
            x = (synaptic_loc - len0) / sec.L
            return sec, x
        h.pop_section()
        len0 = len1


def _apical_random_synapse(HCell, synaptic_loc, distal_apic_x0):
    """
    Return a random location in the apical tree of this cell.

    Codes taken from Eyal et al. 2018 and adapted for our purposes.

    Parameters
    ----------
    synaptic_loc : array
        synaptic locations.
    distal_apic_x0: float
        closest distance, in um, of the distal apical tree to the soma.

    Returns
    -------
    sec : list
        List with the name of the sections where synapses will be inserted.
    x : list
        List of seg.x values associated with sec.

    """
    len0 = 0
    len1 = 0
    for sec in HCell.apical:
        if h.distance(0, sec=sec) > distal_apic_x0:
            len1 += sec.L
            if len1 >= synaptic_loc:
                x = (synaptic_loc - len0) / sec.L
                return sec, x
            h.pop_section()
            len0 = len1


def _oblique_random_synapse(HCell, synaptic_loc, oblique_end):
    """
    Return a random location in the apical tree of this cell.

    Codes taken from Eyal et al. 2018 and adapted for our purposes.

    Parameters
    ----------
    synaptic_loc : array
        synaptic locations.
    oblique_end: float
        distance in um relative to the soma representing the aproximate end
        of the oblique dendrites.

    Returns
    -------
    sec : list
        List with the name of the sections where synapses will be inserted.
    x : list
        List of seg.x values associated with sec.

    """
    len0 = 0
    len1 = 0
    for sec in HCell.apical:
        if h.distance(sec(1)) < oblique_end:
            len1 += sec.L
            if len1 >= synaptic_loc:
                x = (synaptic_loc - len0) / sec.L
                return sec, x
            h.pop_section()
            len0 = len1


def _random_synapse(HCell, rd, section, total_L, dist_x=None):
    """
    Return a random location in the neuron model.

    Codes taken from Eyal et al. 2018 and adapted for our purposes.

    Parameters
    ----------
    rd : NEURON random obj
        NEURON random number generator.
    section : str
        Name of the section of the neuron where the synaptic inputs arrive.
    total_L : int or np.array
        total dendritic length.
    dist_x : float, um
        distance to the soma representing the end or beginning of the oblique
        and distal apical dendritic regions, respectively.

    Returns
    -------
    sec : list
        List with the name of the sections where synapses will be inserted.
    x : list
        List of seg.x values associated with sec.

    use to choose a synaptic location out of the uniform
    distribution of dendritic locations that give the same probability
    to any point on the dendritic tree.
    note that just choosing segments randomly would ignore
    the segments physical length and would bias more synapses
    on shorter segments

    """
    if section == "dend":
        synaptic_loc = rd.uniform(0, total_L)
        return _basal_random_synapse(HCell, synaptic_loc)
    elif section == "apic":
        synaptic_loc = rd.uniform(0, total_L)
        return _apical_random_synapse(HCell, synaptic_loc, dist_x)
    else:
        synaptic_loc = rd.uniform(0, total_L)
        return _oblique_random_synapse(HCell, synaptic_loc, dist_x)


def fill_clustered_synapses_vectors(
    cell, number_of_clusters, rd, section, L5PC_model=True
):
    """
    Choose random number_of_clusters locations.

    The center of each cluster is chosen randomly, and synapses
    are distributed within the cluster length, CLUSTER_L,
    which is 20 um in this work.

    Codes taken from Eyal et al. 2018 and adapted for our purposes.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    number_of_clusters : int
        total number of synaptic cluster to be inserted.
    rd : NEURON random obj
        NEURON random number generator.
    section : str
        Name of the section of the neuron where the synaptic inputs arrive.
    L5PC_model : str, optional
        Pyramidal cell model being used. Default is True, L5 PC model.

    Returns
    -------
    syn_segments : list
        DESCRIPTION.

    """
    from params_synapses import get_cluster_parameters

    cluster_parameters = get_cluster_parameters()

    if L5PC_model:
        """Hay et al. L5 PC model."""
        oblique_end = 420  # um, end of the oblique dendrites | dist relative
        # to the soma
        distal_apic_x0 = 620  # um, beginning of the distal apical tree
        # relative to the soma
        number_of_basal_clusters = round(2 * number_of_clusters / 3)
        HCell = neuron.h
    else:
        """Eyal et al. L3 PC model."""
        oblique_end = 250  # um, end of the oblique dendrites | dist relative
        # to the soma
        distal_apic_x0 = 258.7  # um, beginning of the distal apical tree
        # relative to the soma
        number_of_basal_clusters = round((0.34 / (0.34 + 0.18)) * number_of_clusters)
        HCell = cell.template

    if section == "dend":
        total_L = sum([sec.L for sec in HCell.basal])
        number_of_clusters = number_of_basal_clusters
        distance = None

    elif section == "oblq":
        total_L = 0
        for sec in HCell.apical:
            if h.distance(sec(1)) < oblique_end:
                # print(sec)
                # print(h.distance(0, sec=sec))
                total_L = total_L + sec.L
        number_of_clusters = number_of_clusters - number_of_basal_clusters
        distance = oblique_end

    else:
        total_L = 0
        for sec in HCell.apical:
            if h.distance(0, sec=sec) > distal_apic_x0:
                # print(sec)
                # print(h.distance(0, sec=sec))
                total_L = total_L + sec.L
        distance = distal_apic_x0

    sec_name = []
    seg_x = []
    for i in range(number_of_clusters):
        sec, X_center = _random_synapse(HCell, rd, section, total_L, distance)
        for j in range(cluster_parameters["CLUSTER_SIZE"]):
            if sec.L < cluster_parameters["CLUSTER_L"]:
                x = rd.uniform(0, 1)
            elif X_center < cluster_parameters["CLUSTER_L"] / sec.L:
                x = rd.uniform(0, cluster_parameters["CLUSTER_L"] / sec.L)
            elif X_center > (1 - cluster_parameters["CLUSTER_L"] / sec.L):
                x = rd.uniform(1 - cluster_parameters["CLUSTER_L"] / sec.L, 1)
            else:  # the standard case
                x = rd.uniform(
                    X_center - cluster_parameters["CLUSTER_L"] / 2.0 / sec.L,
                    X_center + cluster_parameters["CLUSTER_L"] / 2.0 / sec.L,
                )

            # syn_segments.append(sec(x))
            sec_name.append(sec)
            seg_x.append(x)

    return sec_name, seg_x  # syn_segments


def insert_clustered_synapses_l5(
    cell,
    cell_index,
    section,
    sec_name,
    seg_x,
    syn_model_params,
    syn_params,
    input_params,
    runNumb,
    SEED,
    bkg_inputs,
):
    """
    Insert synapses onto the sections 'sec_name(seg_x)' of L5 PC model.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    cell_index : int
        Index of the cell being simulated.
    section : str
        area of the neuron where synapses will be insereted
    sec_name : list
        Name of the sections where synapses will be inserted.
    seg_x : list
        seg.x values within the sec_name where synapses will be inserted.
    syn_model_params : dict
        Dictionary of parameters of synaptic mechanisms. E.g. 'e', 'tau', etc.
    syn_params : dict
        Dictionary of parameters for the pre-synaptic spike trains generation.
    input_params : dict
        Dictionary specifying time-locked inputs properties.
    runNumb : int
        number of the trial being simulated.
    SEED : float
        seed for the random number generator.
    bkg_inputs : str
        name of the sections for background input only.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    idx_int : int
        Cell index where the synaptic input arrives.

    """
    from params_synapses import get_cluster_parameters

    cluster_parameters = get_cluster_parameters()

    # total number of synapses
    number_of_synapses = len(sec_name)

    np.random.seed(SEED * cell_index + runNumb)  # synapse activation SEED
    n_synchro = round(number_of_synapses / cluster_parameters["CLUSTER_SIZE"])

    """ Background Synaptic Activation == Baseline Inputs. """
    if any(
        section in sec_list_bkg_inputs for sec_list_bkg_inputs in bkg_inputs
    ):  # Basal Inputs Only
        if input_params["increase_baseline_inputs"]:
            """ Increase background inputs after saccade onset. """

            syn_params["args"].update(
                {
                    "tstop": (
                        input_params["target_times"][cell_index, runNumb - 1]
                        + input_params["warmup_period"]
                        + input_params["saccade_times"][cell_index, runNumb - 1]
                    ),
                }
            )

            input_params["args_post_event_bkg_input"].update(
                {
                    "tstart": (
                        input_params["target_times"][cell_index, runNumb - 1]
                        + input_params["warmup_period"]
                        + input_params["saccade_times"][cell_index, runNumb - 1]
                    )
                }
            )
            spktrain_clusteri = input_params["spTimesFun"](
                n=n_synchro, **input_params["args_post_event_bkg_input"]
            )
            spiketrain = list(
                itertools.chain.from_iterable(
                    itertools.repeat(x, cluster_parameters["CLUSTER_SIZE"])
                    for x in spktrain_clusteri
                )
            )

        spkTimes_clusteri = syn_params["spTimesFun"](n=n_synchro, **syn_params["args"])
        spkTimes_all = list(
            itertools.chain.from_iterable(
                itertools.repeat(x, cluster_parameters["CLUSTER_SIZE"])
                for x in spkTimes_clusteri
            )
        )

    """ Generate time-locked synaptic inputs. """
    if input_params["Gaussian_Input"]:

        syn_loc_time_input = np.random.randint(
            0, number_of_synapses, round(number_of_synapses / 2)
        )

        """ Basal time-locked inputs """
        if (section == "dend") and input_params["dend_syn"]:
            stimLocked_spkTimes = []

            for i in range(len(input_params["a_dend"])):
                stimLocked_spkTimesi = st.skewnorm.rvs(
                    a=input_params["a_dend"][i],
                    loc=(
                        input_params["mu_dend"][i]
                        + input_params["target_times"][cell_index, runNumb - 1]
                        + input_params["saccade_times"][cell_index, runNumb - 1]
                    ),
                    scale=input_params["sigma_dend"][i],
                    size=(number_of_synapses, input_params["n_dend"][i]),
                )

                if i == 0:
                    stimLocked_spkTimes = stimLocked_spkTimesi
                else:
                    stimLocked_spkTimes = np.append(
                        stimLocked_spkTimes, stimLocked_spkTimesi, axis=1
                    )

        """ Oblique time-locked inputs """
        if (section == "oblq") and input_params["oblq_syn"]:
            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a_oblq"],
                loc=(
                    input_params["mu_oblq"]
                    + input_params["target_times"][cell_index, runNumb - 1]
                    + input_params["saccade_times"][cell_index, runNumb - 1]
                ),
                scale=input_params["sigma_oblq"],
                size=(number_of_synapses, input_params["n_oblq"]),
            )

        """ Distal Apical time-locked inputs """
        if (section == "apic") and input_params["apic_syn"]:
            # activating the whole distal tree
            if not input_params["sep_hemitrees"]:
                stimLocked_spkTimes = st.skewnorm.rvs(
                    a=input_params["a_apic"]["a_apic"],
                    loc=(
                        input_params["mu_apic"]["mu_apic"]
                        + input_params["target_times"][cell_index, runNumb - 1]
                        + input_params["saccade_times"][cell_index, runNumb - 1]
                    ),
                    scale=input_params["sigma_apic"]["sigma_apic"],
                    size=(number_of_synapses, input_params["n_apic"]["n_apic"]),
                )
            # distinct activation of apical hemi-trees
            else:
                if input_params["hemitree1"]:
                    stimLocked_spkTimes_hemitree1 = []

                    for i in range(len(input_params["a_apic"]["a_h1"])):
                        stimLocked_spkTimes_hemitree1i = st.skewnorm.rvs(
                            a=input_params["a_apic"]["a_h1"][i],
                            loc=(
                                input_params["mu_apic"]["mu_h1"][i]
                                + input_params["target_times"][cell_index, runNumb - 1]
                                + input_params["saccade_times"][cell_index, runNumb - 1]
                            ),
                            scale=input_params["sigma_apic"]["sigma_h1"][i],
                            size=(
                                number_of_synapses,
                                input_params["n_apic"]["n_h1"][i],
                            ),
                        )
                        if i == 0:
                            stimLocked_spkTimes_hemitree1 = (
                                stimLocked_spkTimes_hemitree1i
                            )
                        else:
                            stimLocked_spkTimes_hemitree1 = np.append(
                                stimLocked_spkTimes_hemitree1,
                                stimLocked_spkTimes_hemitree1i,
                                axis=1,
                            )

                if input_params["hemitree2"]:
                    stimLocked_spkTimes_hemitree2 = []

                    for i in range(len(input_params["a_apic"]["a_h2"])):
                        stimLocked_spkTimes_hemitree2i = st.skewnorm.rvs(
                            a=input_params["a_apic"]["a_h2"][i],
                            loc=(
                                input_params["mu_apic"]["mu_h2"][i]
                                + input_params["target_times"][cell_index, runNumb - 1]
                                + input_params["saccade_times"][cell_index, runNumb - 1]
                            ),
                            scale=input_params["sigma_apic"]["sigma_h2"][i],
                            size=(
                                number_of_synapses,
                                input_params["n_apic"]["n_h2"][i],
                            ),
                        )
                        if i == 0:
                            stimLocked_spkTimes_hemitree2 = (
                                stimLocked_spkTimes_hemitree2i
                            )
                        else:
                            stimLocked_spkTimes_hemitree2 = np.append(
                                stimLocked_spkTimes_hemitree2,
                                stimLocked_spkTimes_hemitree2i,
                                axis=1,
                            )

                if input_params["hemitree3"]:
                    stimLocked_spkTimes_hemitree3 = []

                    for i in range(len(input_params["a_apic"]["a_h3"])):
                        stimLocked_spkTimes_hemitree3i = st.skewnorm.rvs(
                            a=input_params["a_apic"]["a_h3"][i],
                            loc=(
                                input_params["mu_apic"]["mu_h3"][i]
                                + input_params["target_times"][cell_index, runNumb - 1]
                                + input_params["saccade_times"][cell_index, runNumb - 1]
                            ),
                            scale=input_params["sigma_apic"]["sigma_h3"][i],
                            size=(
                                number_of_synapses,
                                input_params["n_apic"]["n_h3"][i],
                            ),
                        )
                        if i == 0:
                            stimLocked_spkTimes_hemitree3 = (
                                stimLocked_spkTimes_hemitree3i
                            )
                        else:
                            stimLocked_spkTimes_hemitree3 = np.append(
                                stimLocked_spkTimes_hemitree3,
                                stimLocked_spkTimes_hemitree3i,
                                axis=1,
                            )

    idx_int = []

    for i in range(number_of_synapses):
        idx = utils.get_index(cell, seg_x[i], sec_name[i])
        syn_model_params.update({"idx": int(idx)})
        idx_int.append(int(idx))

        spiketimes = []
        # Baseline Inputs
        if any(section in sec_list_bkg_inputs for sec_list_bkg_inputs in bkg_inputs):
            spiketimes = spkTimes_all[i]

            if input_params["increase_baseline_inputs"]:
                spiketimes = np.append(spiketimes, spiketrain[i])

        # Time-locked Inputs
        if input_params["Gaussian_Input"]:
            if (
                ((section == "dend") and input_params["dend_syn"])
                or (
                    (section == "apic")
                    and input_params["apic_syn"]
                    and not input_params["sep_hemitrees"]
                )
                or ((section == "oblq") and input_params["oblq_syn"])
            ) and any(syn_loc_time_input == i):
                spiketimes = np.append(spiketimes, stimLocked_spkTimes[i])

            if (
                (section == "apic")
                and input_params["apic_syn"]
                and input_params["sep_hemitrees"]
            ) and any(syn_loc_time_input == i):
                if (idx >= 617 and idx <= 691) and input_params["hemitree1"]:
                    spiketimes = np.append(spiketimes, stimLocked_spkTimes_hemitree1[i])
                if (idx >= 695 and idx <= 737) and input_params["hemitree2"]:
                    spiketimes = np.append(spiketimes, stimLocked_spkTimes_hemitree2[i])
                if (idx >= 738 and idx <= 894) and input_params["hemitree3"]:
                    spiketimes = np.append(spiketimes, stimLocked_spkTimes_hemitree3[i])

            spiketimes.sort()

        if len(spiketimes) == 0:
            continue
        # Create synapse(s) and setting times using the Synapse class
        # in LFPy
        s = LFPy.Synapse(cell, **syn_model_params)
        s.set_spike_times(spiketimes)

    return cell, idx_int


def insert_clustered_synapses_l3(
    cell,
    cell_index,
    section,
    sec_name,
    seg_x,
    syn_model_params,
    syn_params,
    input_params,
    runNumb,
    SEED,
    # bkg_inputs,
):
    """
    Insert synapses onto the sections 'sec_name(seg_x)' of L3 PC model.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    cell_index : int
        Index of the cell being simulated.
    section : str
        area of the neuron where synapses will be insereted
    sec_name : list
        Name of the sections where synapses will be inserted.
    seg_x : list
        seg.x values within the sec_name where synapses will be inserted.
    syn_model_params : dict
        Dictionary of parameters of synaptic mechanisms. E.g. 'e', 'tau', etc.
    syn_params : dict
        Dictionary of parameters for the pre-synaptic spike trains generation.
    input_params : dict
        Dictionary specifying time-locked inputs properties.
    runNumb : int
        number of the trial being simulated.
    SEED : float
        seed for the random number generator.
    bkg_inputs : str
        name of the sections for background input only.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    idx_int : int
        Cell index where the synaptic input arrives.

    """
    from params_synapses import get_cluster_parameters

    cluster_parameters = get_cluster_parameters()

    # total number of synapses
    number_of_synapses = len(sec_name)

    np.random.seed(SEED * cell_index + runNumb)  # synapse activation SEED
    n_synchro = round(number_of_synapses / cluster_parameters["CLUSTER_SIZE"])

    """ Background Synaptic Activation == Baseline Inputs. """
    # if bkg_inputs == section or bkg_inputs == "allsec":
    spkTimes_clusteri = syn_params["spTimesFun"](n=n_synchro, **syn_params["args"])
    spkTimes_all = list(
        itertools.chain.from_iterable(
            itertools.repeat(x, cluster_parameters["CLUSTER_SIZE"])
            for x in spkTimes_clusteri
        )
    )

    """ Generate time-locked synaptic inputs. """
    if input_params["Gaussian_Input"]:

        """ Basal dendrites time-locked inputs. """
        if (section == "dend") and input_params["dend_syn"]:
            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a"] if input_params["Skewed_dist"] else 0,
                loc=(
                    input_params["mu"]
                    + input_params["target_times"][cell_index, runNumb - 1]
                    + input_params["saccade_times"][cell_index, runNumb - 1]
                ),
                scale=input_params["sigma"],
                size=(number_of_synapses, input_params["n_dend"]),
            )

            syn_loc_time_input = np.random.randint(
                0, number_of_synapses, round(number_of_synapses / 2)
            )

        """" Distal apical time-locked inputs. """
        if (section == "apic") and input_params["apic_syn"]:
            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a"] if input_params["Skewed_dist"] else 0,
                loc=(
                    input_params["mu"]
                    + input_params["target_times"][cell_index, runNumb - 1]
                    + input_params["saccade_times"][cell_index, runNumb - 1]
                ),
                scale=input_params["sigma"],
                size=(number_of_synapses, input_params["n_apic"]),
            )

    np.random.seed(SEED)
    idx_int = []

    for i in range(number_of_synapses):
        idx = utils.get_index(cell, seg_x[i], sec_name[i])
        syn_model_params.update({"idx": int(idx)})
        idx_int.append(int(idx))

        spiketimes = spkTimes_all[i]

        if input_params["Gaussian_Input"]:
            if (section == "apic") and input_params["apic_syn"]:
                spiketimes = np.append(spiketimes, stimLocked_spkTimes[i])

            if ((section == "dend") and input_params["dend_syn"]) and any(
                syn_loc_time_input == i
            ):
                spiketimes = np.append(spiketimes, stimLocked_spkTimes[i])

            spiketimes.sort()
        # Create synapse(s) and setting times using the Synapse class
        # in LFPy
        s = LFPy.Synapse(cell, **syn_model_params)
        s.set_spike_times(spiketimes)

    return cell, idx_int


def insert_dist_synapses(
    cell,
    cell_index,
    syn_model_params,
    syn_params,
    input_params,
    runNumb,
    SEED,
    L5PC_model=True,
):
    """
    Find n compartments to insert synapses onto.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    cell_index : int
        Index of the cell being simulated.
    syn_model_params : dict
        Dictionary of parameters of synaptic mechanisms. E.g. 'e', 'tau', etc.
    syn_params : dict
        Dictionary containing the specifics of synaptic inputs including
        number of synapses, n; distribution to draw pre-synaptic spike times,
        spTimesFun; arguments for the spike times distribution; etc.
    input_params: dict
        Dictionary specifying parameters for extra time-locked stimulation.
    runNumb: int
        number of trial being simulated.
    L5PC_model : str, optional
        Pyramidal cell model being used. Default is True, L5 PC model.

    Returns
    -------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    idx_int : int
        Cell index where the synaptic input arrives.

    """
    if L5PC_model:
        """Hay et al. L5 PC model."""
        oblique_end = 620  # um, end of the oblique dendrites | dist relative
        # to the soma
        distal_apic_x0 = 620  # um, beginning of the distal apical tree
        # relative to the soma
        number_of_basal_syn = round((2 / 3) * syn_params["n"])
    else:
        """Eyal et al. L3 PC model."""
        oblique_end = 250  # um, end of the oblique dendrites | dist relative
        # to the soma
        distal_apic_x0 = 258.7  # um, beginning of the distal apical tree
        # relative to the soma
        number_of_basal_syn = round((0.34 / (0.34 + 0.18)) * syn_params["n"])

    if syn_params["section"] == "apic" or (
        syn_params["section"] == "dend" and not L5PC_model
    ):

        np.random.seed(SEED * cell_index + runNumb)
        spkTimes_all = syn_params["spTimesFun"](n=syn_params["n"], **syn_params["args"])

        if input_params["Gaussian_Input"]:

            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a"],
                loc=input_params["mu"],
                scale=input_params["sigma"],
                size=(syn_params["n"], input_params["n"]),
            )

        np.random.seed(SEED)
        if syn_params["section"] == "apic":

            if input_params.get("rhythmicity_sim") == None:
                idx = cell.get_rand_idx_area_norm(
                    section="allsec",
                    nidx=syn_params["n"],
                    z_min=(distal_apic_x0 + np.median(cell.z[0, :])),
                )
            else:
                idx = cell.get_rand_idx_area_norm(section="apic", nidx=syn_params["n"],)
        else:
            if input_params.get("rhythmicity_sim") == None:
                idx = cell.get_rand_idx_area_norm(
                    section=["dend", "apic"],
                    nidx=syn_params["n"],
                    z_max=(oblique_end + np.median(cell.z[0, :])),
                )
            else:
                idx = cell.get_rand_idx_area_norm(section="dend", nidx=syn_params["n"],)

    elif syn_params["section"] == "oblq" and L5PC_model:

        np.random.seed(SEED * cell_index + runNumb)
        spkTimes_all = syn_params["spTimesFun"](
            n=(syn_params["n"] - number_of_basal_syn), **syn_params["args"]
        )

        if input_params["Gaussian_Input"]:

            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a"],
                loc=input_params["mu"],
                scale=input_params["sigma"],
                size=((syn_params["n"] - number_of_basal_syn), input_params["n"]),
            )

        np.random.seed(SEED)
        idx = cell.get_rand_idx_area_norm(
            section="apic",
            nidx=(syn_params["n"] - number_of_basal_syn),
            z_max=(oblique_end + np.median(cell.z[0, :])),
        )

    elif syn_params["section"] == "dend" and L5PC_model:

        np.random.seed(SEED * cell_index + runNumb)
        spkTimes_all = syn_params["spTimesFun"](
            n=number_of_basal_syn, **syn_params["args"]
        )

        if input_params["Gaussian_Input"]:
            stimLocked_spkTimes = st.skewnorm.rvs(
                a=input_params["a"],
                loc=input_params["mu"],
                scale=input_params["sigma"],
                size=(number_of_basal_syn, input_params["n"]),
            )

        np.random.seed(SEED)
        idx = cell.get_rand_idx_area_norm(section="dend", nidx=number_of_basal_syn)

    idx_int = []

    for i in range(len(idx)):
        syn_model_params.update({"idx": int(idx[i])})
        idx_int.append(int(idx[i]))

        spiketimes = spkTimes_all[i]
        if input_params["Gaussian_Input"]:
            spiketimes = np.append(spiketimes, stimLocked_spkTimes[i])
            spiketimes.sort()
        # Create synapse(s) and setting times using the Synapse class
        # in LFPy
        s = LFPy.Synapse(cell, **syn_model_params)
        s.set_spike_times(spiketimes)

    return cell, idx_int


def get_idx_name_obj(cell, idx=np.array([0], dtype=int)):
    """
    Return NEURON convention name of segments with index idx.

    The returned argument are an array with corresponding segment idx, and
    an object with the corresponding section name and position along the
    section, like; [neuron.h.soma[0](0.5)]

    Taken from LFPy cell.get_idx_name and adapted by Beatriz Herrera to
    return the section name as an object instead of as a string.
    May 26, 2021

    Parameters
    ----------
    idx: ndarray, dtype int
        segment indices, must be between 0 and cell.totnsegs

    Returns
    -------
    ndarray, dtype=object
        tuples with section names of segments
    """
    # ensure idx is array-like, or convert
    if isinstance(idx, int) or np.int64:
        idx = np.array([idx])
    elif len(idx) == 0:
        return
    else:
        idx = np.array(idx).astype(int)

    # ensure all idx are valid
    if np.any(idx >= cell.totnsegs):
        wrongidx = idx[np.where(idx >= cell.totnsegs)]
        raise Exception("idx %s >= number of compartments" % str(wrongidx))

    # create list of seg names:
    allsegnames = []
    allsegidx = []
    segidx = 0
    for sec in cell.allseclist:
        for seg in sec:
            allsegidx.append(segidx)
            allsegnames.append(sec(seg.x))
            segidx += 1

    return (
        np.array(allsegidx, dtype=list)[idx][0],
        np.array(allsegnames, dtype=list)[idx][0],
    )


def get_index(cell, seg_x, HCell_sec):
    """
    Return index of 'HCell_sec' Hoc.Obj section in LFPy.Cell Obj.

    Parameters
    ----------
    cell : TYPE
        DESCRIPTION.
    seg_x : TYPE
        DESCRIPTION.
    HCell_sec : TYPE
        DESCRIPTION.

    Returns
    -------
    i : TYPE
        DESCRIPTION.

    """
    # HCell = getattr(neuron.h,str(cell.template)[:-1])
    i = 0
    names = []
    for sec in cell.allseclist:
        # print(sec.name())
        for seg in sec:
            names.append(sec(seg.x))
            if sec.name() == HCell_sec.name() and (
                h.distance(sec(seg.x)) == h.distance(HCell_sec(seg_x))
            ):
                # print('Sec %s = %s' % (sec.name(), HCell_sec.name()))
                # print('Sec %d = %d' % (h.distance(sec(seg.x)),
                # h.distance(HCell_sec(seg_x))))
                break
            i += 1
        if sec.name() == HCell_sec.name():
            break
    return i
