"""
PyAPX - Sampling Module
Provides random sampling and Bayesian optimization sampling functionality
"""

import random
import os
import pandas as pd
import numpy as np


def run_random_sampling(max_structure_id=None, current_sample_id=None):
    """
    Execute random sampling function

    Args:
        max_structure_id (int): Maximum structure ID
        current_sample_id (int): Current sample ID for on-the-fly seeding

    Returns:
        int or dict: Selected structure ID, or on-the-fly candidate payload
    """
    if not os.path.exists("candidates.csv"):
        from .on_the_fly import sample_random_on_the_fly_structure

        return sample_random_on_the_fly_structure(current_sample_id)

    if max_structure_id is None:
        return None

    # Random sampling
    sampled_id = random.randint(0, max_structure_id)

    return sampled_id


def run_physbo_sampling(current_sample_id):
    """
    Execute PHYSBO sampling function

    Args:
        current_sample_id (int): Current sample ID to use as random seed

    Returns:
        int: Selected structure ID
    """
    from .utils import apx_print

    if not os.path.exists("candidates.csv"):
        apx_print("PHYSBO sampling requires candidates.csv; use OPTIMIZER = botorch for on-the-fly Bayesian sampling")
        return None

    try:
        import physbo
    except ImportError as e:
        apx_print(f"Error importing PHYSBO dependency: {e}")
        return None

    # Read PHYSBO settings from apx.in
    from .utils import read_physbo_setting, load_encoded_data_from_cache
    score, num_rand_basis = read_physbo_setting()

    # Load pre-encoded data from cache
    success, X = load_encoded_data_from_cache()
    if not success:
        return None

    # Center the data for PHYSBO
    X = physbo.misc.centering(X.astype(np.float32))

    # Load initial data from samples.csv
    samples_df = pd.read_csv("samples.csv")
    t_initial = -samples_df['energy'].values  # Convert minimization problem to maximization problem
    calculated_ids = samples_df['structure_id'].values

    # Create policy with initial data
    policy = physbo.search.discrete.policy(test_X=X, initial_data=[calculated_ids, t_initial])
    #policy.set_seed(current_sample_id)  # Use current_sample_id as random seed

    # Perform Bayesian search using settings from apx.in
    if score == "TS":
        # Use num_rand_basis only for Thompson Sampling
        actions = policy.bayes_search(
            max_num_probes=1,
            simulator=None,
            score=score,
            interval=0,
            num_rand_basis=num_rand_basis
        )
    else:
        # For EI, PI, etc., don't use num_rand_basis
        actions = policy.bayes_search(
            max_num_probes=1,
            simulator=None,
            score=score,
            interval=0
        )

    return actions[0]


def run_botorch_sampling(current_sample_id):
    """
    Execute BoTorch GP sampling.
    """
    from .botorch_sampler import run_botorch_sampling as _run_botorch_sampling
    return _run_botorch_sampling(current_sample_id)
