#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os

def generate_candidates(num_sites=23, num_hydrogens=15, output_file='./candidates.csv'):
    """
    Generate atomic configurations for candidates.csv
    
    Parameters:
    -----------
    num_sites : int
        Number of sites in the atomic configuration (default: 23)
    num_hydrogens : int
        Number of hydrogen atoms (default: 15)
    output_file : str
        Output file path for candidates.csv (default: './candidates.csv')
    
    Returns:
    --------
    pandas.DataFrame
        Generated candidates with structure_id and site columns
    """
    
    # Initialize with single site configurations [0] and [1]
    X = np.array([[0], [1]])
    
    # Generate all possible configurations by adding sites
    for i in range(1, num_sites):
        h = np.shape(X)[0]
        # Add 0 (None) to all existing configurations
        X0 = np.hstack([X, np.zeros([h, 1])])
        # Add 1 (H) to all existing configurations
        X1 = np.hstack([X, np.ones([h, 1])])
        # Combine both sets
        X = np.vstack([X0, X1])
        
        # Filter configurations with too many hydrogens
        X = np.delete(X, np.sum(X, axis=1) > num_hydrogens, axis=0)
    
    # Filter configurations with too few hydrogens
    X = np.delete(X, np.sum(X, axis=1) < num_hydrogens, axis=0)
    
    # Create DataFrame with proper format
    # Convert 0/1 to H/E for atomic configuration
    atomic_configs = []
    for i, config in enumerate(X):
        # Convert 0->E (Empty), 1->H
        atomic_config = ['H' if site == 1 else 'E' for site in config]
        atomic_configs.append(atomic_config)
    
    # Create DataFrame with structure_id and site columns
    df_data = []
    for i, config in enumerate(atomic_configs):
        row = {'structure_id': i}  # 0-based indexing
        for j, site in enumerate(config, 1):
            row[f'site_{j}'] = site
        df_data.append(row)
    
    df = pd.DataFrame(df_data)
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    print(f"Candidates saved to {output_file}")
    
    return df

def main():
    """Main function to generate candidates"""
    # Check if output directory exists
    output_dir = os.path.dirname('./candidates.csv')
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Generate candidates
    candidates_df = generate_candidates()
    
    # Display first few rows
    print("\nFirst 5 candidates:")
    print(candidates_df.head().to_string(index=False))
    
    # Display last few rows
    print("\nLast 5 candidates:")
    print(candidates_df.tail().to_string(index=False))

if __name__ == "__main__":
    main() 