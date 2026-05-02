import os
import time
from datetime import datetime
import torch
from ase.io import read
from chgnet.model.dynamics import StructOptimizer
from pyapx.energy_evaluator import create_qe_input

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

def run_chgnet_calculation(sample_id, structure_id):
    input_filename = f"chgnet_temp_{sample_id}.in"
    input_path = f"dft_calc/{input_filename}"
    relaxer = StructOptimizer(use_device=DEVICE)

    try:
        atomic_config = create_qe_input(structure_id, input_filename)
        atoms = read(input_path, format='espresso-in')
        
        #timer start
        start_time = time.perf_counter()
        start_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{start_dt}] Relax START for Sample {sample_id}...", flush=True)

        print(f"Relaxing structure for Sample {sample_id}...", flush=True)
        result = relaxer.relax(atoms, steps=200, fmax=0.02, verbose=False)
        final_energy = result["trajectory"].energies[-1]

        print(f"CHGNet Relaxed Energy - Sample {sample_id}: {final_energy:.6f} eV")

        # timer stop
        end_time = time.perf_counter()
        end_dt = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        elapsed_time = end_time - start_time
        
        print(f"[{end_dt}] Relax END for Sample {sample_id}", flush=True)
        print(f"--> Elapsed time: {elapsed_time:.2f} seconds", flush=True)
        del relaxer
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        return True, final_energy, atomic_config, None

    except Exception as e:
        print(f"Error in CHGNet calculation: {e}")
        return False, None, None, str(e)
        
    finally:
        if os.path.exists(input_path):
            os.remove(input_path)
