"""
PyAPX - Energy Evaluator Module
Manages energy calculations using QE and other DFT codes
"""

import os
import re
import time
import pandas as pd
from .utils import apx_print, read_parallel_command

def run_qe_calculation_and_extract_energy(input_file, output_file, candidate_id):
    """
    Execute Quantum ESPRESSO calculation, wait for completion, and extract total energy
    
    Args:
        input_file (str): Quantum ESPRESSO input filename
        output_file (str): Quantum ESPRESSO output filename
        candidate_id (int): Candidate structure ID
    
    Returns:
        tuple: (total_energy, error_message)
            - total_energy: Extracted total energy (in Ry), None if error
            - error_message: Error message, None if successful
    """
    
    # Execute calculation in dft_calc directory
    input_path = f"dft_calc/{input_file}"
    output_path = f"dft_calc/{output_file}"
    
    # Ensure dft_calc directory exists
    os.makedirs("dft_calc", exist_ok=True)
    
    # Grant execution permission to file
    os.chmod(input_path, 0o755)
    
    # Execute Quantum ESPRESSO calculation (explicitly specify file path)
    apx_print("Waiting for DFT calculation completion...")
    # Read parallel command prefix from apx.in and construct full command
    parallel_command_prefix = read_parallel_command()
    qe_command_template = f"{parallel_command_prefix} pw.x -npool 1 -inp {{input_path}} > {{output_path}}"
    qe_command = qe_command_template.format(input_path=input_path, output_path=output_path)
    os.system(qe_command)
    

    # Wait for calculation completion
    job_done = False
    convergence_failed = False
    fatal_error = False
    coordinates_missing = False
    
    while not job_done:
        with open(output_path, 'r') as f:
            contents = f.read()
            if "JOB DONE" in contents:
                job_done = True
            if re.search(r'convergence NOT achieved after (\d+) iterations: stopping', contents):
                convergence_failed = True
            if "dE0s is positive which should never happen" in contents:
                fatal_error = True
            # Check for missing coordinates
            if job_done and not re.search(r'Begin final coordinates.*?End final coordinates', contents, re.DOTALL):
                coordinates_missing = True
        time.sleep(1)  # Wait 1 second before checking again

    # Skip processing if error messages exist, even with JOB DONE
    if convergence_failed:
        return None, "Skipped due to DFT calculation not converging"
    if fatal_error:
        return None, "Skipped due to DFT calculation not completing normally"
    if coordinates_missing:
        return None, "Skipped due to DFT calculation not completing normally"
    
    # Process when JOB DONE exists
    if job_done and not (convergence_failed or fatal_error or coordinates_missing):
        apx_print(f"Confirmed completion of DFT calculation from {output_file}.")

    os.chmod(output_path, 0o755)

    # Extract energy
    total_energy = None
    try:
        with open(output_path, 'r') as f:
            content = f.read()
            # Extract total energy (get the last occurrence)
            energy_matches = re.findall(r'!.*total energy.*=\s*([-\d.]+)\s*Ry', content)
            if energy_matches:
                total_energy = float(energy_matches[-1])
            else:
                return None, "Skipped due to DFT calculation not completing normally"
    except Exception as e:
        return None, f"Skipped due to DFT calculation not completing normally: {e}"
    
    return total_energy, None


def create_qe_input(candidate_id, input_file):
    """
    Read specified structure ID row from candidates.csv and generate QE input file
    
    Args:
        candidate_id (int): Structure ID to read from candidates.csv
        input_file (str): QE input filename to generate
    
    Returns:
        list: Atomic configuration data
    """
    # Read files from execution directory
    csv_file = "candidates.csv"  # Input CSV file
    input_template = "qe_template.in"  # Existing Quantum ESPRESSO template file
    
    # Ensure dft_calc directory exists
    dft_calc_dir = "dft_calc"
    os.makedirs(dft_calc_dir, exist_ok=True)

    # Read only the row with specific structure ID
    df = pd.read_csv(csv_file)
    candidate_row = df[df.iloc[:, 0] == candidate_id]
    if candidate_row.empty:
        raise ValueError(f"Structure ID {candidate_id} not found")
    
    atomic_config_row = candidate_row.iloc[0, 1:].values  # Get atomic configuration data only, excluding structure ID

    with open(input_template, "r") as template_file:
        content = template_file.readlines()

    # Identify ATOMIC_POSITIONS section
    positions_start = None
    for i, line in enumerate(content):
        if line.strip().startswith("ATOMIC_POSITIONS"):
            positions_start = i + 1
            break

    if positions_start is None:
        raise ValueError("ATOMIC_POSITIONS section not found in the template file.")

    # Get data until next section
    positions = []
    template_atoms = []
    for line in content[positions_start:]:
        if line.strip() == "" or any(line.strip().startswith(keyword) for keyword in ["K_POINTS", "CELL_PARAMETERS", "ATOMIC_SPECIES"]):
            break
        line_parts = line.strip().split()
        if len(line_parts) >= 1:
            template_atoms.append(line_parts[0])  # First column element
            positions.append("   ".join(line_parts[1:]))  # Coordinate part

    # Count X's and check if they match the number of sites in candidates.csv
    x_count = sum(1 for atom in template_atoms if atom.upper() == "X")
    if x_count != len(atomic_config_row):
        raise ValueError(f"Number of X's in ATOMIC_POSITIONS ({x_count}) does not match number of sites in candidates.csv ({len(atomic_config_row)})")

    # Create output file in dft_calc directory
    output_path = f"dft_calc/{input_file}"
    with open(output_path, "w") as file:
        file.writelines(content[:positions_start])  # Template header part

        atomic_config_index = 0  # Index for atomic_config_row
        for template_atom, position in zip(template_atoms, positions):
            # Replace X (wildcard) elements in template with elements from candidates.csv sequentially
            if template_atom.upper() == "X":
                atom = atomic_config_row[atomic_config_index]
                atomic_config_index += 1
            else:
                atom = template_atom
            file.write(f"{atom}   {position}\n")
    
    apx_print(f"Generated new Quantum ESPRESSO input file: {output_path}")
    
    return atomic_config_row  # Return atomic configuration data


def run_qe_calculation(sample_id, structure_id):
    """
    Execute QE calculation for structure
    
    Args:
        sample_id (int): Sample ID
        structure_id (int): Structure ID
    
    Returns:
        tuple: (success, total_energy, atomic_config, error_message)
            - success: Whether calculation was successful
            - total_energy: Total energy in Ry (if successful)
            - atomic_config: Atomic configuration (if successful)
            - error_message: Error message (if failed)
    """
    # Create dft_calc directory
    dft_calc_dir = "dft_calc"
    os.makedirs(dft_calc_dir, exist_ok=True)
    
    # Check existence of qe_template.in file
    qe_template_path = "qe_template.in"
    if not os.path.exists(qe_template_path):
        error_msg = "qe_template.in file not found. Please place qe_template.in in the execution directory."
        return False, None, None, error_msg
    
    # Execute QE calculation
    input_file = f"qe_sample_{sample_id}.in"
    output_file = f"qe_sample_{sample_id}.out"
    
    try:
        # Generate QE input file for candidate structure
        atomic_config = create_qe_input(structure_id, input_file)
        
        # Execute QE calculation and extract energy
        total_energy, error_message = run_qe_calculation_and_extract_energy(input_file, output_file, structure_id)
        
        if total_energy is not None:
            return True, total_energy, atomic_config, None
        else:
            return False, None, None, error_message
            
    except Exception as e:
        return False, None, None, str(e)

# VASP calculation function to be implemented in the future
def run_vasp_calculation(sample_id, structure_id):
    """
    Execute VASP calculation for structure (to be implemented in the future)
    
    Args:
        sample_id (int): Sample ID
        structure_id (int): Structure ID
    
    Returns:
        tuple: (success, total_energy, atomic_config, error_message)
            - success: Whether calculation was successful
            - total_energy: Total energy (if successful)
            - atomic_config: Atomic configuration (if successful)
            - error_message: Error message (if failed)
    """
    apx_print(f"VASP calculation is not currently implemented")
    return sample_id + 1
 