"""
PyAPX - Utility Module
Provides common utility functions
"""

import os
import pandas as pd

def apx_print(message):
    """
    Print message with apx prefix
    """
    print(f"(apx) {message}", flush=True)

def get_atom_types_and_num_sites_from_candidates():
    """
    Function to get atom types and number of sites from structure ID 0 in candidates.csv
    
    Returns:
        tuple: (atom_types, num_sites)
            - atom_types: List of atom types
            - num_sites: Number of sites
    """
    # Read only the row with structure ID 0
    df = pd.read_csv("candidates.csv", nrows=1)
    first_row = df.iloc[0, 1:].values  # Atomic configuration of structure ID 0 (excluding structure ID column)
    
    # Extract unique atom types and sort
    atom_types = list(set(first_row))
    atom_types.sort()
    
    # Number of sites is the length of atomic configuration
    num_sites = len(first_row)
    
    return atom_types, num_sites

def read_card(file_path, card_name):
    """
    Function to read a specific card from input file
    
    Args:
        file_path (str): Path to input file
        card_name (str): Card name
    
    Returns:
        list: List of card content lines
    """
    card_data = []
    with open(file_path, "r") as file:
        lines = file.readlines()
        inside_card = False

        for line in lines:
            stripped_line = line.strip()

            # Detect start of specified card
            if stripped_line.startswith(card_name):
                inside_card = True
                continue

            # Stop reading when encountering another card or empty line
            # Detect uppercase card names (but exclude lowercase data values)
            if inside_card and (not stripped_line or 
                              (stripped_line.isupper() and not stripped_line.startswith('#'))):
                break

            # Read data lines
            if inside_card:
                card_data.append(stripped_line)

    return card_data

def get_max_structure_id():
    """
    Function to get maximum structure ID from candidates.csv
    
    Returns:
        int: Maximum structure ID (None if failed)
    """
    try:
        candidates_path = "candidates.csv"
        with open(candidates_path, "r") as file:
            # Skip header
            next(file)
            # Get last line
            for line in file:
                pass
            max_structure_id = int(line.strip().split(',')[0])
        return max_structure_id
    except FileNotFoundError:
        apx_print("candidates.csv not found")
        return None
    except (ValueError, IndexError):
        apx_print("candidates.csv format is incorrect")
        return None

def read_card_value(file_path, card_name):
    """
    Read a KEY = VALUE style setting from an input file.

    This parser intentionally does exact key matching and strips comments, so
    settings such as

        SCORE = TS    # acquisition function
        BOTORCH_DEVICE=cuda

    are read as "TS" and "cuda" respectively.
    """
    target = str(card_name).strip().upper()

    try:
        with open(file_path, "r") as file:
            for raw_line in file:
                # Remove inline comments and surrounding whitespace.
                line = raw_line.split("#", 1)[0].strip()
                if not line or "=" not in line:
                    continue

                key, value = line.split("=", 1)
                if key.strip().upper() == target:
                    return value.strip()

        return None
    except Exception as e:
        apx_print(f"Error reading card value {card_name}: {e}")
        return None


def read_energy_evaluator():
    """
    Read energy evaluator settings from apx.in
    
    Returns:
        str: Energy evaluator name
            - "qe": Quantum ESPRESSO
            - "vasp": VASP
            - "custom:module.function": User-defined function
              Example: "custom:my_evaluator.run_calculation"
    """
    try:
        evaluator = read_card_value("apx.in", "ENERGY_EVALUATOR")
        if evaluator:
            return evaluator
        else:
            return "qe"  # Default
    except Exception as e:
        apx_print(f"Error reading energy evaluator settings: {e}")
        return "qe"  # Default

def read_sampling_sequence():
    """
    Read sampling settings from apx.in
    
    Returns:
        list: List of sampling tasks
    """
    try:
        tasks = []
        
        # Read random sampling settings
        random_sampling = read_card_value("apx.in", "RANDOM_SAMPLING")
        if random_sampling:
            try:
                num_iterations = int(random_sampling)
                tasks.append(("RANDOM_SAMPLING", num_iterations))
            except (ValueError, IndexError):
                apx_print("RANDOM_SAMPLING value is incorrect")
        
        # Read Bayesian sampling settings
        bayes_sampling = read_card_value("apx.in", "BAYES_SAMPLING")
        if bayes_sampling:
            try:
                num_iterations = int(bayes_sampling)
                tasks.append(("BAYES_SAMPLING", num_iterations))
            except (ValueError, IndexError):
                apx_print("BAYES_SAMPLING value is incorrect")
        
        return tasks
    
    except Exception as e:
        apx_print(f"Error reading sampling settings: {e}")
        return []

def get_current_sample_id():
    """
    Function to get current sample ID from samples.csv (for continuation)
    
    Returns:
        int: Next sample ID (0 if file doesn't exist)
    """
    samples_path = "samples.csv"
    
    if not os.path.exists(samples_path):
        return 0
    
    try:
        with open(samples_path, "r") as file:
            lines = file.readlines()
            if len(lines) > 1:  # If there is data other than header
                # Get sample ID from last line and add 1
                last_line = lines[-1].strip()
                if last_line:
                    last_sample_id = int(last_line.split(',')[0])
                    return last_sample_id + 1
        return 0
    except Exception as e:
        apx_print(f"Error reading existing samples.csv: {e}")
        return 0


def ensure_samples_csv_header():
    """
    Function to create samples.csv header (only if file doesn't exist)
    """
    samples_csv_path = "samples.csv"
    
    if not os.path.exists(samples_csv_path):
        if os.path.exists("candidates.csv"):
            # Get number of sites from structure ID 0
            _, num_sites = get_atom_types_and_num_sites_from_candidates()
        else:
            from .on_the_fly import load_on_the_fly_space, use_on_the_fly_if_no_candidates

            if not use_on_the_fly_if_no_candidates():
                raise FileNotFoundError("candidates.csv not found")
            num_sites = load_on_the_fly_space().num_sites
        
        header_parts = ["sample_id", "structure_id"]
        for i in range(num_sites):
            header_parts.append(f"site_{i+1}")
        header_parts.append("energy")
        
        with open(samples_csv_path, "w") as samples_file:
            samples_file.write(','.join(header_parts) + '\n')

def read_encode_type_setting():
    """
    Read encode type setting from apx.in.

    ENCODE_TYPE is the preferred key. ENCODE is accepted as a legacy alias so
    existing input files using "ENCODE = NAmod" still work.

    Returns:
        tuple: (encode_type, weight)
            - encode_type (str): Encoding method ("OH", "NA", or "NAmod")
            - weight (float): Weight for neighbor atom encoding
    """
    try:
        encode_type_raw = (
            read_card_value("apx.in", "ENCODE_TYPE")
            or read_card_value("apx.in", "ENCODE")
            or "OH"
        )

        encode_key = encode_type_raw.strip().lower()
        encode_map = {
            "oh": "OH",
            "na": "NA",
            "namod": "NAmod",
        }
        if encode_key not in encode_map:
            raise ValueError(
                f"Unsupported ENCODE_TYPE/ENCODE: {encode_type_raw}. "
                "Choose OH, NA, or NAmod."
            )
        encode_type = encode_map[encode_key]

        weight_str = read_card_value("apx.in", "WEIGHT")
        weight = float(weight_str) if weight_str else 0.0

        return encode_type, weight
    except Exception as e:
        apx_print(f"Error reading encode type setting: {e}")
        return "OH", 0.0


def read_encoding_setting():
    """
    Read encoding setting from apx.in
    
    Returns:
        bool: Whether to perform encoding preprocessing before sampling
    """
    try:
        encoding_str = read_card_value("apx.in", "ENCODING")
        encoding = encoding_str.lower() == "true" if encoding_str else False
        return encoding
    except Exception as e:
        apx_print(f"Error reading encoding setting: {e}")
        return False  # Default value

def read_dimension_reduction_setting():
    """
    Read dimension reduction settings from apx.in.

    For PCA_N_COMPONENTS, the following forms are supported:
      - integer, e.g. 30
      - float in (0, 1), e.g. 0.95, for sklearn's explained-variance target
      - none/default, to use sklearn's default

    Returns:
        tuple: (use_dimension_reduction, method, params)
            - use_dimension_reduction (bool): Whether to use dimension reduction
            - method (str): Method ("PCA" or "AUTOENCODER")
            - params (dict): Parameters for the selected method
    """

    def _parse_bool(value, default=False):
        if value is None:
            return default
        return str(value).strip().lower() in ["true", "t", "yes", "y", "1", "on"]

    def _parse_optional_int(value, name, default=None):
        if value is None or str(value).strip() == "":
            return default
        text = str(value).strip().lower()
        if text in ["none", "default"]:
            return None
        try:
            return int(text)
        except ValueError as exc:
            raise ValueError(f"{name} must be an integer, none, or omitted: {value}") from exc

    def _parse_pca_n_components(value):
        if value is None or str(value).strip() == "":
            return None
        text = str(value).strip().lower()
        if text in ["none", "default"]:
            return None
        if text == "mle":
            return "mle"

        try:
            number = float(text)
        except ValueError as exc:
            raise ValueError(
                f"PCA_N_COMPONENTS must be an integer, a float in (0, 1), "
                f"mle, none, or omitted: {value}"
            ) from exc

        if 0.0 < number < 1.0:
            return number
        if number >= 1.0 and abs(number - round(number)) < 1.0e-12:
            return int(round(number))
        raise ValueError(
            f"PCA_N_COMPONENTS must be an integer >= 1 or a float in (0, 1): {value}"
        )

    try:
        use_dimension_reduction_str = read_card_value("apx.in", "USE_DIMENSION_REDUCTION")
        use_dimension_reduction = _parse_bool(use_dimension_reduction_str, default=False)

        if use_dimension_reduction:
            method = read_card_value("apx.in", "DIMENSION_REDUCTION_METHOD") or "PCA"
            method = method.strip().upper()

            random_state_str = read_card_value("apx.in", "AUTOENCODER_RANDOM_STATE")
            random_state = _parse_optional_int(random_state_str, "AUTOENCODER_RANDOM_STATE", default=42)

            if method == "PCA":
                pca_n_components_str = read_card_value("apx.in", "PCA_N_COMPONENTS")
                pca_n_components = _parse_pca_n_components(pca_n_components_str)

                # PCA_RANDOM_STATE is only used by sklearn for randomized/arpack
                # solvers, but keeping it fixed is useful for reproducibility.
                pca_random_state_str = read_card_value("apx.in", "PCA_RANDOM_STATE")
                pca_random_state = _parse_optional_int(pca_random_state_str, "PCA_RANDOM_STATE", default=None)

                params = {
                    'n_components': pca_n_components,
                    'random_state': pca_random_state
                }
            elif method == "AUTOENCODER":
                latent_dim_str = read_card_value("apx.in", "AUTOENCODER_LATENT_DIM")
                latent_dim = int(latent_dim_str) if latent_dim_str else 20

                hidden_layers_str = read_card_value("apx.in", "AUTOENCODER_HIDDEN_LAYERS")
                if hidden_layers_str:
                    hidden_layers = tuple(int(x.strip()) for x in hidden_layers_str.split(','))
                else:
                    hidden_layers = (128, 64)

                learning_rate_str = read_card_value("apx.in", "AUTOENCODER_LEARNING_RATE")
                learning_rate = float(learning_rate_str) if learning_rate_str else 0.001

                epochs_str = read_card_value("apx.in", "AUTOENCODER_EPOCHS")
                epochs = int(epochs_str) if epochs_str else 100

                batch_size_str = read_card_value("apx.in", "AUTOENCODER_BATCH_SIZE")
                batch_size = int(batch_size_str) if batch_size_str else 1024

                params = {
                    'latent_dim': latent_dim,
                    'hidden_layers': hidden_layers,
                    'learning_rate': learning_rate,
                    'epochs': epochs,
                    'batch_size': batch_size,
                    'random_state': random_state
                }
            else:
                apx_print(f"Unknown dimension reduction method: {method}")
                return False, None, None
        else:
            method = None
            params = None

        return use_dimension_reduction, method, params
    except Exception as e:
        apx_print(f"Error reading dimension reduction settings: {e}")
        return False, None, None


def read_autoencoder_setting():
    """
    Read auto encoder settings from apx.in (for backward compatibility)
    
    Returns:
        tuple: (use_autoencoder, autoencoder_params)
            - use_autoencoder (bool): Whether to use auto encoder
            - autoencoder_params (dict): Auto encoder parameters
    """
    use_dimension_reduction, method, params = read_dimension_reduction_setting()
    
    if use_dimension_reduction and method and method.upper() == "AUTOENCODER":
        return True, params
    else:
        return False, None

def read_physbo_setting():
    """
    Read PHYSBO settings from apx.in.

    Returns:
        tuple: (score, num_rand_basis)
            - score (str): Score function ("TS", "EI", "PI")
            - num_rand_basis (int): Number of random basis (only used for "TS")
    """
    try:
        score = read_card_value("apx.in", "SCORE")
        score = score.strip().upper() if score else "TS"

        if score == "TS":
            num_rand_basis_str = read_card_value("apx.in", "NUM_RAND_BASIS")
            num_rand_basis = int(num_rand_basis_str) if num_rand_basis_str else 3000
        else:
            num_rand_basis = 3000

        return score, num_rand_basis
    except Exception as e:
        apx_print(f"Error reading PHYSBO settings: {e}")
        return "TS", 3000


def read_botorch_setting():
    """
    Read BoTorch settings from apx.in using a compact schema.
    """
    truthy = {"true", "t", "yes", "y", "1", "on"}
    bool_parser = lambda value: str(value).strip().lower() in truthy
    lower = lambda value: str(value).strip().lower()
    upper = lambda value: str(value).strip().upper()
    text = lambda value: str(value).strip()
    schema = [
        ("score", ("SCORE",), upper, "TS"),
        ("device", ("BOTORCH_DEVICE",), text, "auto"),
        ("dtype", ("BOTORCH_DTYPE",), lower, "float64"),
        ("batch_size", ("BOTORCH_BATCH_SIZE",), int, 65536),
        ("maxiter", ("BOTORCH_MAXITER",), int, 100),
        ("xi", ("BOTORCH_XI",), float, 0.0),
        ("input_transform", ("BOTORCH_INPUT_TRANSFORM", "BOTORCH_SCALING"), lower, "unit_cube"),
        ("gp_kernel", ("BOTORCH_GP_KERNEL",), lower, "default"),
        ("standardize_y", ("BOTORCH_STANDARDIZE_Y",), bool_parser, False),
        ("discrete_chunk_size", ("BOTORCH_DISCRETE_CHUNK_SIZE",), int, 4096),
        ("hamming_lengthscale", ("BOTORCH_HAMMING_LENGTHSCALE",), float, 0.2),
        ("ot_lengthscale", ("BOTORCH_OT_LENGTHSCALE",), float, 0.5),
        ("ot_atom_mismatch_penalty", ("BOTORCH_OT_ATOM_MISMATCH_PENALTY",), float, 1.0),
        ("ot_sinkhorn_epsilon", ("BOTORCH_OT_SINKHORN_EPSILON",), float, 0.05),
        ("ot_sinkhorn_iterations", ("BOTORCH_OT_SINKHORN_ITERATIONS",), int, 30),
        ("ot_chunk_size", ("BOTORCH_OT_CHUNK_SIZE",), int, 1024),
        ("use_local_env", ("BOTORCH_USE_LOCAL_ENV",), bool_parser, False),
        ("local_env_type", ("BOTORCH_LOCAL_ENV_TYPE",), text, "NA"),
        ("env_distance", ("BOTORCH_ENV_DISTANCE",), lower, "l1"),
        ("env_lengthscale", ("BOTORCH_ENV_LENGTHSCALE",), float, 0.2),
        ("ot_env_mismatch_penalty", ("BOTORCH_OT_ENV_MISMATCH_PENALTY",), float, 1.0),
        ("local_env_cache", ("BOTORCH_LOCAL_ENV_CACHE",), text, "local_env_candidates.pkl"),
        ("sa_screening", ("BOTORCH_SA_SCREENING",), bool_parser, False),
        ("sa_score", ("BOTORCH_SA_SCORE",), upper, ""),
        ("sa_initial_pool_size", ("BOTORCH_SA_INITIAL_POOL_SIZE",), int, 65536),
        ("sa_chains", ("BOTORCH_SA_CHAINS",), int, 1024),
        ("sa_steps", ("BOTORCH_SA_STEPS",), int, 500),
        ("sa_eval_batch_size", ("BOTORCH_SA_EVAL_BATCH_SIZE",), int, 4096),
        ("sa_initial_temperature", ("BOTORCH_SA_INITIAL_TEMPERATURE",), float, 1.0),
        ("sa_final_temperature", ("BOTORCH_SA_FINAL_TEMPERATURE",), float, 0.01),
        ("sa_random_fraction", ("BOTORCH_SA_RANDOM_FRACTION",), float, 0.05),
        ("sa_swap_neighbors", ("BOTORCH_SA_SWAP_NEIGHBORS",), bool_parser, True),
        ("on_the_fly_sa_restarts", ("ON_THE_FLY_SA_RESTARTS", "SA_NUM_RESTARTS"), int, None),
        ("on_the_fly_sa_steps", ("ON_THE_FLY_SA_STEPS", "SA_NUM_STEPS"), int, None),
        (
            "on_the_fly_sa_initial_temperature",
            ("ON_THE_FLY_SA_INITIAL_TEMPERATURE", "SA_INITIAL_TEMPERATURE"),
            float,
            None,
        ),
        (
            "on_the_fly_sa_cooling_rate",
            ("ON_THE_FLY_SA_COOLING_RATE", "SA_COOLING_RATE"),
            float,
            None,
        ),
    ]

    settings = {}
    for key, names, parser, default in schema:
        value = None
        source = names[0]
        for name in names:
            value = read_card_value("apx.in", name)
            if value is not None and str(value).strip() != "":
                source = name
                break
        if value is None or str(value).strip() == "":
            settings[key] = default
            continue
        try:
            settings[key] = parser(value)
        except (TypeError, ValueError):
            apx_print(f"Invalid {source}: {value}; using default {default}")
            settings[key] = default

    local_env_type_key = str(settings["local_env_type"]).strip().lower()
    local_env_type_map = {
        "na": "NA",
        "namod": "NAmod",
    }
    if local_env_type_key not in local_env_type_map:
        raise ValueError(
            f"Unsupported BOTORCH_LOCAL_ENV_TYPE: {settings['local_env_type']}. "
            "Choose NA or NAmod."
        )
    settings["local_env_type"] = local_env_type_map[local_env_type_key]

    env_distance_key = str(settings["env_distance"]).strip().lower()
    env_distance_map = {
        "l1": "l1",
        "manhattan": "l1",
        "l2": "l2",
        "euclidean": "l2",
    }
    if env_distance_key not in env_distance_map:
        raise ValueError(
            f"Unsupported BOTORCH_ENV_DISTANCE: {settings['env_distance']}. "
            "Choose l1 or l2."
        )
    settings["env_distance"] = env_distance_map[env_distance_key]

    if settings["env_lengthscale"] <= 0:
        raise ValueError("BOTORCH_ENV_LENGTHSCALE must be positive")
    if settings["ot_env_mismatch_penalty"] < 0:
        raise ValueError("BOTORCH_OT_ENV_MISMATCH_PENALTY must be non-negative")

    return settings


def read_optimizer():
    """
    Read optimizer settings from apx.in
    
    Returns:
        str: Optimizer name ("physbo", etc.)
    """
    try:
        optimizer = read_card_value("apx.in", "OPTIMIZER")
        if optimizer:
            return optimizer
        else:
            return "physbo"  # Default
    except Exception as e:
        apx_print(f"Error reading optimizer settings: {e}")
        return "physbo"  # Default

def read_parallel_command():
    """
    Read parallel execution command prefix from apx.in (part before DFT code executable)
    
    Returns:
        str: Parallel execution command prefix (e.g., "mpiexec", "mpirun -np 8")
             Default: "mpiexec"
    """
    try:
        parallel_command = read_card_value("apx.in", "PARALLEL_COMMAND")
        if parallel_command:
            return parallel_command.strip()
        else:
            # Default command prefix
            return "mpiexec"
    except Exception as e:
        apx_print(f"Error reading parallel command settings: {e}")
        # Default command prefix
        return "mpiexec"

def load_encoded_data_from_cache():
    """
    Load pre-encoded data from cache
    
    Returns:
        tuple: (success, data)
            - success (bool): True if successful, False otherwise
            - data (np.ndarray): Encoded data if successful, None otherwise
    """
    import os
    import pickle
    import numpy as np
    
    cache_file = "encoded_candidates.pkl"
    if not os.path.exists(cache_file):
        apx_print(f"Error: Cache file {cache_file} not found.")
        apx_print("Please run encoding first by setting ENCODING = True")
        return False, None
    
    try:
        with open(cache_file, 'rb') as f:
            data = pickle.load(f)
        return True, data
    except Exception as e:
        apx_print(f"Error loading cache file: {e}")
        return False, None

def read_use_initial_density():
    """
    Read USE_INITIAL_DENSITY setting from apx.in
    
    Returns:
        bool: Whether to use initial density (True if USE_INITIAL_DENSITY = True, False otherwise)
    """
    try:
        use_initial_density_str = read_card_value("apx.in", "USE_INITIAL_DENSITY")
        use_initial_density = use_initial_density_str.lower() == "true" if use_initial_density_str else False
        return use_initial_density
    except Exception as e:
        apx_print(f"Error reading USE_INITIAL_DENSITY setting: {e}")
        return False  # Default value

def read_qe_prefix_and_outdir():
    """
    Read prefix and outdir from qe_template.in
    
    Returns:
        tuple: (prefix, outdir)
            - prefix (str): Prefix value from qe_template.in (None if not found)
            - outdir (str): Outdir value from qe_template.in (None if not found)
    """
    try:
        qe_template_path = "qe_template.in"
        if not os.path.exists(qe_template_path):
            return None, None
        
        prefix = None
        outdir = None
        
        with open(qe_template_path, "r") as file:
            lines = file.readlines()
            inside_control = False
            
            for line in lines:
                stripped_line = line.strip()
                
                # Skip empty lines and comments
                if not stripped_line or stripped_line.startswith("#"):
                    continue
                
                # Detect start of &control section
                if stripped_line.startswith("&control") or stripped_line.startswith("&CONTROL"):
                    inside_control = True
                    continue
                
                # Stop reading when encountering another section or end of control
                if inside_control:
                    if stripped_line.startswith("/"):
                        break
                    if stripped_line.startswith("&"):
                        break
                    
                    # Parse prefix (handle both 'prefix = ...' and 'prefix=' formats)
                    if "prefix" in stripped_line.lower() and "=" in stripped_line:
                        # Remove comment if present
                        line_without_comment = stripped_line.split("#")[0].strip()
                        if "=" in line_without_comment:
                            value = line_without_comment.split("=", 1)[1].strip()
                            # Remove quotes if present
                            prefix = value.strip("'\"")
                    
                    # Parse outdir (handle both 'outdir = ...' and 'outdir=' formats)
                    if "outdir" in stripped_line.lower() and "=" in stripped_line:
                        # Remove comment if present
                        line_without_comment = stripped_line.split("#")[0].strip()
                        if "=" in line_without_comment:
                            value = line_without_comment.split("=", 1)[1].strip()
                            # Remove quotes if present
                            outdir = value.strip("'\"")
        
        return prefix, outdir
    except Exception as e:
        apx_print(f"Error reading prefix and outdir from qe_template.in: {e}")
        return None, None
