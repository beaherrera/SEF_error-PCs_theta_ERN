# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:21:25 2022

@author: Beatriz Herrera

This script generates the tone onset times for the simulations using 
the experimental data distribution.

"""
import numpy as np
import scipy.stats as st
from scipy import io
from os.path import join

SEED = 12
num_trials = 116
POPULATION_SIZE = 1000

# load experimental distribution of target times
event_times = io.loadmat(join("Data", "event_times_rel2target.mat"))
xk_tone = event_times["edges_toneS_time"][:, 1:] - 0.5

pk_tone_Go = event_times["p_toneS_times_Go"]
pk_tone_NC = event_times["p_toneS_times_NC"]

tone_dist_Go = st.rv_discrete(name="tone_dist_Go", values=(xk_tone, pk_tone_Go))
tone_dist_NC = st.rv_discrete(name="tone_dist_NC", values=(xk_tone, pk_tone_NC))

tone_times_Go = np.tile(
    tone_dist_Go.rvs(size=(1, num_trials), random_state=SEED), (POPULATION_SIZE, 1)
)
tone_times_NC = np.tile(
    tone_dist_NC.rvs(size=(1, num_trials), random_state=SEED), (POPULATION_SIZE, 1)
)

save_events = {
    "tone_times_Go": tone_times_Go,
    "tone_times_NC": tone_times_NC,
}

# io.savemat(join("../", "src_matlab", "data", "sim_data",
#                  "tone_times.mat"), save_events)
