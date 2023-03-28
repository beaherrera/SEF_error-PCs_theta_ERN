# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:12:11 2022

@author:  Beatriz Herrera

Compartmental model of rat thick-tufted L5bPC developed by Hay
et al. (2011), including modifications of voltage gated calcium channel
densities as in Shai et al. (2015) and Ih channel density distribution as
in Labarrera et al. (2018).

active declarations of the model taken and adpated from Goldenberg and Segev
(2021): Burst Control to work on LFPy.

"""
from __future__ import division
from os.path import join
import neuron
from neuron import h
import sys

# import LFPy

""" Functions """


def biophys_active():
    """
    L5PCbiophys3 - Hay et al. 2011 model with modifications by Goldenberg and
    Segev 2021.

    Parameters
    ----------
    **kwargs : dict

    Returns
    -------
    None.

    """
    # loading mechanisms
    model_pth = join("cell_models", "HayModel")
    pth = join(model_pth, "mod")
    if not hasattr(h, "Ca_LVAst"):
        if "win32" in sys.platform:
            h.nrn_load_dll(pth + "/nrnmech.dll")
        else:
            neuron.load_mechanisms(pth)

    for sec in neuron.h.allsec():
        sec.insert("pas")
        sec.cm = 1.0
        sec.Ra = 100.0
        sec.e_pas = -90.0

    for sec in neuron.h.axon:
        sec.insert("Im")
        sec.insert("Ca_LVAst")
        sec.insert("Ca_HVA")
        sec.insert("SKv3_1")
        sec.insert("SK_E2")
        sec.insert("K_Tst")
        sec.insert("K_Pst")
        sec.insert("Nap_Et2")
        sec.insert("NaTa_t")
        sec.insert("CaDynamics_E2")
        sec.insert("Ih")
        sec.ek = -85
        sec.ena = 50
        sec.gIhbar_Ih = 0.0001 / 2
        sec.g_pas = 3e-5
        sec.gImbar_Im = 0.013322
        sec.decay_CaDynamics_E2 = 277.300774
        sec.gamma_CaDynamics_E2 = 0.000525
        sec.gCa_LVAstbar_Ca_LVAst = 0.000813
        sec.gCa_HVAbar_Ca_HVA = 0.000222
        sec.gSKv3_1bar_SKv3_1 = 0.473799
        sec.gSK_E2bar_SK_E2 = 0.000047
        sec.gK_Tstbar_K_Tst = 0.077274
        sec.gK_Pstbar_K_Pst = 0.188851
        sec.gNap_Et2bar_Nap_Et2 = 0.005834
        sec.gNaTa_tbar_NaTa_t = 3.89618

    for sec in neuron.h.soma:
        sec.insert("Im")
        sec.insert("Ca_LVAst")
        sec.insert("Ca_HVA")
        sec.insert("SKv3_1")
        sec.insert("SK_E2")
        sec.insert("NaTs2_t")
        sec.insert("CaDynamics_E2")
        sec.insert("Ih")
        sec.ek = -85
        sec.ena = 50
        sec.gIhbar_Ih = 0.0001 * 0.75
        sec.g_pas = 3e-5
        sec.gImbar_Im = 0.000008
        sec.decay_CaDynamics_E2 = 294.679571
        sec.gamma_CaDynamics_E2 = 0.000509
        sec.gCa_LVAstbar_Ca_LVAst = 0.000557
        sec.gCa_HVAbar_Ca_HVA = 0.000644
        sec.gSKv3_1bar_SKv3_1 = 0.338029
        sec.gSK_E2bar_SK_E2 = 0.09965
        sec.gNaTs2_tbar_NaTs2_t = 0.998912

    for sec in neuron.h.apic:
        sec.cm = 2
        sec.insert("Ih")
        sec.insert("SK_E2")
        sec.insert("Ca_LVAst")
        sec.insert("Ca_HVA")
        sec.insert("SKv3_1")
        sec.insert("NaTs2_t")
        sec.insert("Im")
        sec.insert("CaDynamics_E2")
        sec.ek = -85
        sec.ena = 50
        sec.decay_CaDynamics_E2 = 35.725651
        sec.gamma_CaDynamics_E2 = 0.000637
        sec.gSK_E2bar_SK_E2 = 0.000002
        sec.gCa_HVAbar_Ca_HVA = 0.000701
        sec.gSKv3_1bar_SKv3_1 = 0.001808
        sec.gNaTs2_tbar_NaTs2_t = 0.021489
        sec.gImbar_Im = 0.00099
        sec.gIhbar_Ih = 0.00015  # 0.00001*1.5
        sec.g_pas = 6e-5

    h.distribute_channels(
        "apic", "gIhbar_Ih", 4, -0.8696, 3.6161, 0.0, 2.0870, 0.00010000000
    )
    h.distribute_channels(
        "apic",
        "gCa_LVAstbar_Ca_LVAst",
        3,
        1.000000,
        0.010000,
        685.000000,
        885.000000,
        0.1419540000 * 1.6,
    )

    for sec in neuron.h.dend:
        sec.cm = 2
        sec.insert("Ih")
        sec.gIhbar_Ih = 0.0001 / 2
        sec.g_pas = 6e-5

    print("L5-PC inserted.")


def active_declarations(cell):
    """
    Set active conductances for Hay et al. 2011 model.

    Parameters
    ----------
    **kwargs : dict


    Returns
    -------
    None.

    """
    h.delete_axon()
    h.geom_nseg()
    h.define_shape()
    exec("biophys_active()")
