# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 15:56:03 2022

@author: Beatriz Herrera

Synapse models parameters and number of synapses per cell.

"""
from __future__ import division

import scipy.stats as st
from LFPy.inputgenerators import get_activation_times_from_distribution


def get_synp_params():
    """
    Define synaptic models parameters.

    Returns
    -------
    dict

    """
    # Parmeters - synapse models
    syn_models_params_AMPA = {
        "e": 0,  # reversal potential
        "syntype": "AMPA",  # conductance based exponential synapse
        "tau_r_AMPA": 0.3,  # Time constant, rise
        "tau_d_AMPA": 1.8,  # Time constant, decay
        "weight": 0.00073027,  # Synaptic weight
        "record_current": False,  # record synaptic currents
    }

    syn_models_params_NMDA = {
        "e": 0,  # reversal potential
        "syntype": "NMDA",  # conductance based exponential synapse
        "tau_r_NMDA": 8.019,  # Time constant, rise
        "tau_d_NMDA": 34.9884,  # Time constant, decay
        "n_NMDA": 0.28011,
        "gama_NMDA": 0.0765685,
        "weight": 0.00131038,  # Synaptic weight
        "record_current": False,  # record synaptic currents
    }

    syn_models_params_GABAA = {
        "e": -80,  # reversal potential
        "syntype": "GABAA",  # conductance based exponential synapse
        "tau_r_GABAA": 0.2,  # Time constant, rise
        "tau_d_GABAA": 1.7,  # Time constant, decay
        "weight": 0.0001,  # Synaptic weight
        "record_current": False,  # record synaptic currents
        # 'record_potential': True
    }

    syn_models_params_GABAB = {
        "Erev": -95,  # reversal potential
        "syntype": "GABAB3",  # conductance based exponential synapse #
        "gmax": 0.0001,  # synaptic weight == maximal conductance
        "weight": 1,
        "record_current": False,  # record synaptic currents
        # 'record_potential': True
    }

    return dict(
        syn_models_params_AMPA=syn_models_params_AMPA,
        syn_models_params_NMDA=syn_models_params_NMDA,
        syn_models_params_GABAA=syn_models_params_GABAA,
        syn_models_params_GABAB=syn_models_params_GABAB,
    )


def get_synp_dist_L3PCs(tstop, **kwargs):
    """
    Define presynaptic inputs - parameters.

    where to insert, how many, and which input statistics.

    Number of synapses and ratios estimated from: Rapan et al. 2021. Neuroimage
    (226): 117574.

    Parameters
    ----------
    tstop : float
        end of the simulation.

    Returns
    -------
    dict

    """
    # - start time of stimulation
    tstart_AMPA = 0
    tstart_NMDA = tstart_AMPA
    tstart_GABAA = 0
    tstart_GABAB = 0

    # number of NMDA synapses in the ablique + basal dendrites, used as
    # reference number
    n_dend_NMDA_L5PCs = 890  # includes basal and oblique dendritic synapses
    n_AMPA_L3PCs = 0.296 * n_dend_NMDA_L5PCs
    n_NMDA_L3PCs = 3.27 * n_dend_NMDA_L5PCs

    insert_synapses_AMPA_args = {
        # 'section' : 'apic',
        "n_dend": int(round((0.18 + 0.39) * n_AMPA_L3PCs)),
        "n_apic": int(round(0.42 * n_AMPA_L3PCs)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_AMPA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_AMPA_dend"]) * 1e3
                    if "r_AMPA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_AMPA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_AMPA_apic"]) * 1e3
                    if "r_AMPA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    insert_synapses_NMDA_args = {
        "n_dend": int(round((0.18 + 0.39) * n_NMDA_L3PCs)),
        "n_apic": int(round(0.42 * n_NMDA_L3PCs)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_NMDA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_NMDA_dend"]) * 1e3
                    if "r_NMDA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_NMDA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_NMDA_apic"]) * 1e3
                    if "r_NMDA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    n_GABA = int(
        round(
            (
                insert_synapses_NMDA_args["n_dend"]
                + insert_synapses_NMDA_args["n_apic"]
                + insert_synapses_AMPA_args["n_apic"]
                + insert_synapses_AMPA_args["n_dend"]
            )
            / 4
        )
    )
    n_dend_GABA = int(round(n_GABA / 4.39))
    n_apic_GABA = n_GABA - n_dend_GABA

    insert_synapses_GABAA_args = {
        "n_dend": int(round(n_dend_GABA / 1.77)),
        "n_apic": int(round(n_apic_GABA / 1.67)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_GABAA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAA_dend"]) * 1e3
                    if "r_GABAA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_GABAA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAA_apic"]) * 1e3
                    if "r_GABAA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    insert_synapses_GABAB_args = {
        "n_dend": int(round(n_dend_GABA - insert_synapses_GABAA_args["n_dend"])),
        "n_apic": int(round(n_apic_GABA - insert_synapses_GABAA_args["n_apic"])),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_GABAB,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAB_dend"]) * 1e3
                    if "r_GABAB_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_GABAB,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAB_apic"]) * 1e3
                    if "r_GABAB_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    return dict(
        insert_synapses_AMPA_args=insert_synapses_AMPA_args,
        insert_synapses_NMDA_args=insert_synapses_NMDA_args,
        insert_synapses_GABAA_args=insert_synapses_GABAA_args,
        insert_synapses_GABAB_args=insert_synapses_GABAB_args,
    )


def get_synp_dist_L5PCs(tstop, **kwargs):
    """
    Define presynaptic inputs - parameters.

    where to insert, how many, and which input statistics.

    Parameters
    ----------
    tstop : float
        end of the simulation.

    Returns
    -------
    dict

    """
    # - start time of stimulation
    tstart_AMPA = 0
    tstart_NMDA = tstart_AMPA
    tstart_GABAA = 0
    tstart_GABAB = 0

    # number of NMDA synapses in the ablique + basal dendrites, used as
    # reference number
    n_dend_NMDA = 890  # includes basal and oblique dendritic synapses

    insert_synapses_AMPA_args = {
        # 'section' : 'apic',
        "n_dend": int(round(0.1045 * n_dend_NMDA)),
        "n_apic": int(round(0.397 * n_dend_NMDA)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_AMPA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_AMPA_dend"]) * 1e3
                    if "r_AMPA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_oblq": dict(
            tstart=tstart_AMPA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_AMPA_oblq"]) * 1e3
                    if "r_AMPA_oblq" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_AMPA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_AMPA_apic"]) * 1e3
                    if "r_AMPA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    insert_synapses_NMDA_args = {
        "n_dend": n_dend_NMDA,
        "n_apic": int(round(1.35 * n_dend_NMDA)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_NMDA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_NMDA_dend"]) * 1e3
                    if "r_NMDA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_oblq": dict(
            tstart=tstart_NMDA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_NMDA_oblq"]) * 1e3
                    if "r_NMDA_oblq" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_NMDA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_NMDA_apic"]) * 1e3
                    if "r_NMDA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    n_GABA = int(
        round(
            (
                n_dend_NMDA
                + insert_synapses_NMDA_args["n_apic"]
                + insert_synapses_AMPA_args["n_apic"]
                + insert_synapses_AMPA_args["n_dend"]
            )
            / 4
        )
    )
    n_dend_GABA = int(round(n_GABA / 4.39))
    n_apic_GABA = n_GABA - n_dend_GABA

    insert_synapses_GABAA_args = {
        "n_dend": int(round(n_dend_GABA / 1.77)),
        "n_apic": int(round(n_apic_GABA / 1.67)),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_GABAA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAA_dend"]) * 1e3
                    if "r_GABAA_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_oblq": dict(
            tstart=tstart_GABAA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAA_oblq"]) * 1e3
                    if "r_GABAA_oblq" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_GABAA,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAA_apic"]) * 1e3
                    if "r_GABAA_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    insert_synapses_GABAB_args = {
        "n_dend": int(round(n_dend_GABA - insert_synapses_GABAA_args["n_dend"])),
        "n_apic": int(round(n_apic_GABA - insert_synapses_GABAA_args["n_apic"])),
        "spTimesFun": get_activation_times_from_distribution,
        "args_dend": dict(
            tstart=tstart_GABAB,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAB_dend"]) * 1e3
                    if "r_GABAB_dend" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_oblq": dict(
            tstart=tstart_GABAB,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAB_oblq"]) * 1e3
                    if "r_GABAB_oblq" in kwargs.keys()
                    else 0
                ),
            ),
        ),
        "args_apic": dict(
            tstart=tstart_GABAB,
            tstop=tstop,
            distribution=st.expon,
            rvs_args=dict(
                loc=0.0,
                scale=(
                    (1 / kwargs["r_GABAB_apic"]) * 1e3
                    if "r_GABAB_apic" in kwargs.keys()
                    else 0
                ),
            ),
        ),
    }

    return dict(
        insert_synapses_AMPA_args=insert_synapses_AMPA_args,
        insert_synapses_NMDA_args=insert_synapses_NMDA_args,
        insert_synapses_GABAA_args=insert_synapses_GABAA_args,
        insert_synapses_GABAB_args=insert_synapses_GABAB_args,
    )


def get_cluster_parameters():
    """
    Return cluster parameters.

    Returns
    -------
    cluster_parameters : dict
        Dictionary containing the clusters' size and length,
        based on the literature, and other parameters.

    """
    cluster_parameters = {
        "CLUSTER_L": 20,  # length of the clusters in um
        "CLUSTER_SIZE": 20,  # number of synapses in a cluster
    }

    return cluster_parameters
