# -*- coding: utf-8 -*-
"""
Created on Tue Sep 28 10:36:18 2021

@author: Beatriz Herrera
"""

from __future__ import division
import os
import posixpath
import numpy as np
from neuron import h


def gen_file_name_L5Pop(stimulusType, **kwargs):
    """
    Generate output files' name for L3 PC Population.

    Parameters
    ----------
    stimulusType : dict
        stimulus definition.
    **kwargs : dict
        dictonary of additional stimulus parameters, input rates.

    Returns
    -------
    detName : str
        file name.

    """
    detName = ""
    if stimulusType["dend_synp"] > 0:
        firingRates_dend = kwargs["rates_dend"]

        detName = detName + ("Dend_r" + str(firingRates_dend["r_NMDA"]))

        if stimulusType["oblq_synp"] > 0 or stimulusType["apic_synp"] > 0:
            detName = detName + "_"

    if stimulusType["oblq_synp"] > 0:
        firingRates_oblq = kwargs["rates_oblq"]

        detName = detName + ("Oblq_" + str(firingRates_oblq["r_NMDA"]))

        if stimulusType["apic_synp"] > 0:
            detName = detName + "_"

    if stimulusType["apic_synp"] > 0:

        firingRates_apic = kwargs["rates_apic"]

        detName = detName + ("Apic_r" + str(firingRates_apic["r_NMDA"]))

    return detName


def gen_file_name_L3Pop(stimulusType, **kwargs):
    """
    Generate output files' name for L3 PC Population.

    Parameters
    ----------
    stimulusType : dict
        stimulus definition.
    **kwargs : dict
        dictonary of additional stimulus parameters, input rates.

    Returns
    -------
    detName : str
        file name.

    """
    detName = ""
    if stimulusType["dend_synp"] > 0:
        firingRates_dend = kwargs["rates_dend"]

        detName = detName + ("Dend_r" + str(firingRates_dend["r_NMDA"]))

        if stimulusType["apic_synp"] > 0:
            detName = detName + "_"
    if stimulusType["apic_synp"] > 0:

        firingRates_apic = kwargs["rates_apic"]

        detName = detName + ("Apic_r" + str(firingRates_apic["r_NMDA"]))

    return detName


def create_output_folder_L5Pop(main_path, POPULATION_SIZE, stimulusType):
    """
    Create output folder for L5 PC Population.

    Parameters
    ----------
    main_path : str
        path to the main folder.
    POPULATION_SIZE : int
        number of neurons in the population.
    stimulusType : dict
        stimulus definition.

    Returns
    -------
    data_folder : str
        path to the output folder.
    """

    data_folder = os.path.join(
        main_path,
        (
            "neurons#"
            + str(POPULATION_SIZE)
            + "_"
            + stimulusType["synapse_type"]
            + "_synp"
        ),
        (
            "StimDend#"
            + str(stimulusType["dend_synp"])
            + "_StimOblq#"
            + str(stimulusType["oblq_synp"])
            + "_StimApic#"
            + str(stimulusType["apic_synp"])
        ),
    )

    if not os.path.isdir(data_folder):
        try:
            os.makedirs(data_folder)
        except OSError:
            print("Creation of the directory %s failed" % data_folder)
        else:
            print("Successfully created the directory %s " % data_folder)

    return data_folder


def create_output_folder_L3Pop(main_path, POPULATION_SIZE, stimulusType):
    """
    Create output folder for L3 PC Population.

    Parameters
    ----------
    main_path : str
        path to the main folder.
    POPULATION_SIZE : int
        number of neurons in the population.
    stimulusType : dict
        stimulus definition.

    Returns
    -------
    data_folder : str
        path to the output folder.
    """

    data_folder = os.path.join(
        main_path,
        (
            "neurons#"
            + str(POPULATION_SIZE)
            + "_"
            + stimulusType["synapse_type"]
            + "_synp"
        ),
        (
            "StimDend#"
            + str(stimulusType["dend_synp"])
            + "_StimApic#"
            + str(stimulusType["apic_synp"])
        ),
    )

    if not os.path.isdir(data_folder):
        try:
            os.makedirs(data_folder)
        except OSError:
            print("Creation of the directory %s failed" % data_folder)
        else:
            print("Successfully created the directory %s " % data_folder)

    return data_folder


def get_index(cell, seg_x, HCell_sec):
    """
    Return the cell ind associated with the section=HCell_sec and seg.x=seg_x.

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    seg_x : int
        value of seg_x.
    HCell_sec : hoc.HocObject
        reference to the NEURON object representing the section of a neuron.

    Returns
    -------
    i : int
        index of the compartment associated with the section=HCell_sec
    and seg.x=seg_x.

    """
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
                #                        h.distance(HCell_sec(seg_x))))
                break
            i += 1
        if sec.name() == HCell_sec.name():
            break
    return i


def get_templatename(f):
    """
    Assess from hoc file the templatename being specified within.

    Arguments
    ---------
    f : file, mode 'r'

    Returns
    -------
    templatename : str

    """
    for line in f.readlines():
        if "begintemplate" in line.split():
            templatename = line.split()[-1]
            print("template {} found!".format(templatename))
            continue

    return templatename


def posixpth(pth):
    """Replace Windows path separators with posix style separators."""
    return pth.replace(os.sep, posixpath.sep)


def get_idx_name_modfun(cell, idx=np.array([0], dtype=int)):
    r"""
    Return NEURON convention name of segments with index idx.

    The returned arguments are an array with corresponding segment idx, and
    an object with the corresponding section name and position along the
    section. E.g., [0] and [neuron.h.soma[0](0.5)].

    Function (cell.get_idx_name()) taken from LFPy.cell class
    (https://github.com/LFPy/LFPy/blob/master/LFPy/cell.py) of the LFPy Python
    module (Hagen, E., Næss, S., Ness, T. V., & Einevoll, G. T. (2018).
            Multimodal modeling of neural network activity: Computing lfp,
            ecog, eeg, and meg signals with lfpy 2.0. Frontiers in
            neuroinformatics, 12, 92.).
    Adapted by Beatriz Herrera to return the section names as objects
    instead of as strings.
    May 26, 2021

    Parameters
    ----------
    cell : obj, LFPy.Cell object
        object built on top of NEURON representing biological neuron.
    idx : ndarray, dtype int, optional
        segment indices, must be between 0 and cell.totnsegs.
        The default is np.array([0], dtype=int).

    Raises
    ------
    Exception
        idx >= cell.totnsegs.

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

    return {
        np.array(allsegidx, dtype=list)[idx][0],
        np.array(allsegnames, dtype=list)[idx][0],
    }

