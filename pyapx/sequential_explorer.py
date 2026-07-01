"""
PyAPX - Sequential Exploration Module
Manages coordination between energy evaluator and sampler, executes exploration loops
"""

import inspect
import os

SAMPLES_CSV_PATH = "samples.csv"


def _unpack_sampled_structure(sampled_structure):
    if isinstance(sampled_structure, dict):
        return (
            sampled_structure.get("structure_id"),
            sampled_structure.get("atomic_config"),
            sampled_structure.get("mode"),
        )
    return sampled_structure, None, None


def _call_energy_calculation_func(energy_calculation_func, sample_id, structure_id, atomic_config):
    if atomic_config is not None:
        try:
            signature = inspect.signature(energy_calculation_func)
            parameters = list(signature.parameters.values())
            positional = [
                parameter
                for parameter in parameters
                if parameter.kind
                in (
                    inspect.Parameter.POSITIONAL_ONLY,
                    inspect.Parameter.POSITIONAL_OR_KEYWORD,
                )
            ]
            accepts_varargs = any(
                parameter.kind == inspect.Parameter.VAR_POSITIONAL
                for parameter in parameters
            )
            if accepts_varargs or len(positional) >= 3:
                return energy_calculation_func(sample_id, structure_id, atomic_config)
        except (TypeError, ValueError):
            pass
    return energy_calculation_func(sample_id, structure_id)

def run_encoding():
    """
    Perform encoding and dimension reduction
    
    Returns:
        bool: True if successful, False otherwise
    """
    from .utils import apx_print, read_encode_type_setting, read_dimension_reduction_setting
    
    print("")
    apx_print("=== Executing ENCODING ===")

    if not os.path.exists("candidates.csv"):
        from .on_the_fly import use_on_the_fly_if_no_candidates

        if use_on_the_fly_if_no_candidates():
            apx_print("candidates.csv not found; skipping enumerated candidate encoding for on-the-fly mode")
            return True
    
    # Read encoding and dimension reduction settings
    encode_type, weight = read_encode_type_setting()
    use_dimension_reduction, reduction_method, reduction_params = read_dimension_reduction_setting()
    
    apx_print(f"Encode type: {encode_type}")
    apx_print(f"Weight: {weight}")
    apx_print(f"Use dimension reduction: {use_dimension_reduction}")
    if use_dimension_reduction:
        apx_print(f"Reduction method: {reduction_method}")
        apx_print(f"Reduction params: {reduction_params}")
        print("")  # Blank line after reduction params
    
    # Perform encoding and dimension reduction
    from .encoder import encode_options
    try:
        # Force fresh calculation by removing cache if it exists
        import os
        cache_file = "encoded_candidates.pkl"
        if os.path.exists(cache_file):
            os.remove(cache_file)
            apx_print(f"Removed existing cache file: {cache_file}")
        
        result = encode_options(encode=encode_type, weight=weight, 
                             use_dimension_reduction=use_dimension_reduction, 
                             reduction_method=reduction_method, 
                             reduction_params=reduction_params)
        
        return True
        
    except Exception as e:
        apx_print(f"Error during encoding: {e}")
        return False

def run_sampling_loop(num_iterations, sampling_type, start_sample_id=0, energy_evaluator="qe"):
    """
    Execute sampling loop function (random or Bayesian)
    
    Args:
        num_iterations (int): Number of sampling iterations
        sampling_type (str): Sampling type ("random" or "bayes")
        start_sample_id (int): Starting sample ID (for continuation)
        energy_evaluator (str): Energy evaluator ("qe", "vasp", etc.)
    
    Note:
        For Bayesian sampling, the specific optimizer is determined by OPTIMIZER setting in apx.in
        Currently supports: "physbo", "botorch"
        Future optimizers can be added by extending the optimizer selection logic
    """
    # Create samples.csv header (only if file doesn't exist)
    from .utils import ensure_samples_csv_header, apx_print
    from .on_the_fly import ON_THE_FLY_MODE, get_candidate_mode, use_on_the_fly_if_no_candidates

    try:
        ensure_samples_csv_header()
    except Exception as e:
        apx_print(f"Failed to initialize samples.csv: {e}")
        return
    
    candidate_mode = get_candidate_mode()
    max_structure_id = None
    if candidate_mode == ON_THE_FLY_MODE:
        if not use_on_the_fly_if_no_candidates():
            apx_print("candidates.csv not found and on-the-fly mode is disabled")
            return
        apx_print("Candidate mode: on_the_fly (candidates.csv not found)")
    else:
        # Get maximum structure ID
        from .utils import get_max_structure_id
        max_structure_id = get_max_structure_id()
        if max_structure_id is None:
            apx_print("Failed to get maximum structure ID")
            return
    
    # Select energy evaluator
    if energy_evaluator.lower() == "qe":
        from .energy_evaluator import run_qe_calculation
        energy_calculation_func = run_qe_calculation
    elif energy_evaluator.lower() == "vasp":
        from .energy_evaluator import run_vasp_calculation
        energy_calculation_func = run_vasp_calculation
    elif energy_evaluator.lower().startswith("custom:"):
        # User-defined energy evaluator
        try:
            # Extract module and function name from "custom:module.function"
            custom_spec = energy_evaluator.lower().replace("custom:", "")
            if "." in custom_spec:
                module_name, function_name = custom_spec.split(".", 1)
                # Import the custom module
                custom_module = __import__(module_name, fromlist=[function_name])
                energy_calculation_func = getattr(custom_module, function_name)
            else:
                apx_print(f"Invalid custom energy evaluator format: {energy_evaluator}")
                apx_print("Expected format: custom:module.function")
                return
        except ImportError as e:
            apx_print(f"Failed to import custom energy evaluator: {e}")
            return
        except AttributeError as e:
            apx_print(f"Custom energy evaluator function not found: {e}")
            return
    else:
        apx_print(f"Unsupported energy evaluator: {energy_evaluator}")
        return
    
    # Sequential processing: sampling -> energy evaluation -> recording cycle
    current_sample_id = start_sample_id
    for i in range(num_iterations):
        print()
        
        # 1. Sampling (different processing for random vs bayes)
        if sampling_type.lower() == "random":
            from .sampler import run_random_sampling
            sampled_structure = run_random_sampling(max_structure_id, current_sample_id)
        elif sampling_type.lower() == "bayes":
            # Check optimizer setting for Bayesian optimization
            from .utils import read_optimizer
            optimizer = read_optimizer()
            if optimizer.lower() == "physbo":
                from .sampler import run_physbo_sampling
                sampled_structure = run_physbo_sampling(current_sample_id)
            elif optimizer.lower() == "botorch":
                from .sampler import run_botorch_sampling
                sampled_structure = run_botorch_sampling(current_sample_id)
            else:
                apx_print(f"Unsupported optimizer for Bayesian sampling: {optimizer}")
                return
        else:
            apx_print(f"Unsupported sampling type: {sampling_type}")
            return
        
        sampled_structure_id, sampled_atomic_config, sampled_mode = _unpack_sampled_structure(sampled_structure)

        if sampled_structure_id is None:
            apx_print(f"Sample {current_sample_id}: Sampling failed")
            current_sample_id += 1
            continue
        
        apx_print(f"Sample {current_sample_id}: Selected structure ID {sampled_structure_id}")
        
        # 2. Energy calculation
        success, total_energy, atomic_config, error_message = _call_energy_calculation_func(
            energy_calculation_func,
            current_sample_id,
            sampled_structure_id,
            sampled_atomic_config,
        )
        
        # 3. Record calculation results
        if success:
            with open(SAMPLES_CSV_PATH, "a") as samples_file:
                # Record in order: sample ID, structure ID, atomic configuration, energy
                atomic_config_to_record = atomic_config if atomic_config is not None else sampled_atomic_config
                atomic_config_str = ','.join(str(atom) for atom in atomic_config_to_record)
                samples_file.write(f"{current_sample_id},{sampled_structure_id},{atomic_config_str},{total_energy}\n")
        else:
            apx_print(f"Sample {current_sample_id}: Skipped - {error_message}")
        
        # 4. Update next sample ID
        current_sample_id += 1
