"""
Custom energy evaluator for H/GaN system using Ising model
"""

import pandas as pd
import os
import numpy as np

def read_atomic_config_from_candidates(structure_id):
    """
    Read atomic configuration from candidates.csv for given structure ID
    
    Args:
        structure_id (int): Structure ID to read from candidates.csv
    
    Returns:
        list: Atomic configuration data
    """
    # Read candidates.csv from current directory
    csv_file = "candidates.csv"
    
    try:
        # Read only the row with specific structure ID
        df = pd.read_csv(csv_file)
        candidate_row = df[df.iloc[:, 0] == structure_id]
        if candidate_row.empty:
            raise ValueError(f"Structure ID {structure_id} not found in candidates.csv")
        
        # Get atomic configuration data only, excluding structure ID
        atomic_config_row = candidate_row.iloc[0, 1:].values.tolist()
        return atomic_config_row
        
    except FileNotFoundError:
        raise ValueError(f"candidates.csv not found at {csv_file}")
    except Exception as e:
        raise ValueError(f"Error reading candidates.csv: {e}")

def calculate_ising_energy(atomic_config, J, h, c):
    """
    Calculate Ising model energy for H/E system
    
    Args:
        atomic_config (list): Atomic configuration from candidates.csv (H/E)
        J (float): Interaction parameter
        h (float): Field parameter
        c (float): Constant energy offset
    
    Returns:
        float: Energy (eV)
    """
    # Convert 23-site atomic configuration to digital values
    # H atoms -> 1, E (empty sites) -> 0
    sites = [0] * 23
    for i, atom in enumerate(atomic_config):
        if atom == "H":
            sites[i] = 1
        # E (empty sites) remain 0
    
    # Convert 0/1 to -1/1 for Ising model
    site = np.array(sites) * 2 - 1
    
    # Initialize interaction terms
    ss = np.zeros(27)
    
    # Define interaction terms for each site (H-H interactions)
    ss[0] = site[0] * (-1 + site[1] + site[4] + site[18] + site[19])
    ss[1] = site[1] * (site[0] + site[2] + site[19] + site[20])
    ss[2] = site[2] * (site[1] + site[3] + site[5] + site[20] + site[21])
    ss[3] = site[3] * (site[2] + site[5] + site[21] + site[22] + 1 + 1)
    ss[4] = site[4] * (site[0] + 1 - 1)
    ss[5] = site[5] * (site[2] + site[3] + site[6] + site[8] + 1)
    ss[6] = site[6] * (site[5] + site[7] + site[8] + site[9] + 1)
    ss[7] = site[7] * (site[6] + 1 + 1)
    ss[8] = site[8] * (site[5] + site[6] + site[9] + site[11])
    ss[9] = site[9] * (site[6] + site[8] + site[11] + site[12] + site[13])
    ss[10] = site[10] * (site[11] + site[15] + site[16])
    ss[11] = site[11] * (site[8] + site[9] + site[10] + site[12] + site[16])
    ss[12] = site[12] * (site[9] + site[11] + site[13] + site[16] + site[17] + site[18])
    ss[13] = site[13] * (site[9] + site[12] + site[18] + site[19])
    ss[14] = site[14] * (site[15] + site[20] + site[21])
    ss[15] = site[15] * (site[10] + site[14] + site[16] + site[21] + site[22])
    ss[16] = site[16] * (site[10] + site[11] + site[12] + site[15] + site[17] + site[22])
    ss[17] = site[17] * (site[12] + site[16] + site[18] + site[22] + 1 - 1)
    ss[18] = site[18] * (site[0] + site[12] + site[13] + site[17] + site[19] - 1)
    ss[19] = site[19] * (site[0] + site[1] + site[13] + site[18] + site[20])
    ss[20] = site[20] * (site[1] + site[2] + site[14] + site[19] + site[21])
    ss[21] = site[21] * (site[2] + site[3] + site[14] + site[15] + site[20] + site[22])
    ss[22] = site[22] * (site[3] + site[15] + site[16] + site[17] + site[21] + 1)
    ss[-1] = +1 * (site[3] + site[17] + site[22] - 1 + 1 + 1)
    ss[-2] = -1 * (site[0] + site[4] + site[17] + site[18] + 1 + 1)
    ss[-3] = +1 * (site[3] + site[5] + site[6] + site[7] + 1 + 1)
    ss[-4] = +1 * (site[4] + site[7] + 1 + 1 - 1)
    
    # Calculate energy: E = -J × Σ(interaction terms) - h × Σ(specific sites) + c
    e = -J * np.sum(ss) - h * (site[1] + site[4] + site[7] + site[8] + site[10] + site[13] + site[14]) + c
    
    return e

def run_ising_calculation(sample_id, structure_id):
    """
    Custom energy calculation function using Ising model for H/E system
    
    Args:
        sample_id (int): Sample ID
        structure_id (int): Structure ID
    
    Returns:
        tuple: (success, total_energy, atomic_config, error_message)
    """
    try:
        # Read atomic configuration from candidates.csv
        atomic_config = read_atomic_config_from_candidates(structure_id)
        
        # Ising model parameters
        J = -0.021863  # eV (H-H interaction parameter)
        h = -0.097453  # eV (field parameter for specific sites)
        c = 0.818264   # eV (constant energy offset)
        
        # Calculate Ising energy for H/E system
        energy = calculate_ising_energy(atomic_config, J, h, c)
        
        # Print energy calculation result
        print(f"Ising energy calculation - Sample {sample_id}, Structure {structure_id}: {energy:.6f} eV")
        print()
        
        return True, energy, atomic_config, None
        
    except Exception as e:
        return False, None, None, str(e)