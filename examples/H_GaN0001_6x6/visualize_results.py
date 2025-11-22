#!/usr/bin/env python3
"""
Simple visualization script for PyAPX results
Creates a plot of sample_id vs energy from samples.csv
"""

import pandas as pd
import matplotlib.pyplot as plt

def visualize_results():
    """
    Create a simple plot of sample_id vs energy
    """
    try:
        # Read results from samples.csv
        df = pd.read_csv('samples.csv')
        
        # Create plot
        plt.figure(figsize=(10, 6))
        plt.plot(df.iloc[:, 0], df.iloc[:, -1], 'bo-', markersize=4)
        plt.xlabel('Sample ID')
        plt.ylabel('Energy')
        plt.title('PyAPX Results: Sample ID vs Energy')
        plt.grid(True, alpha=0.3)
        
        # Save plot
        plt.savefig('results_plot.png', dpi=150, bbox_inches='tight')
        print("Plot saved as 'results_plot.png'")
        
        # Show plot
        #plt.show()
        
    except FileNotFoundError:
        print("Error: samples.csv not found. Please run PyAPX first.")
    except Exception as e:
        print(f"Error creating plot: {e}")

if __name__ == "__main__":
    visualize_results()
