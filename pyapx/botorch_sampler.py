"""
BoTorch sampling workflow for PyAPX.
"""

import os
import pandas as pd
import numpy as np

from .botorch_kernels import (
    build_site_graph_distance,
    create_botorch_model,
    is_discrete_botorch_kernel,
)
from .botorch_sa import run_botorch_sa_screening


def _load_botorch_modules():
    """
    Import BoTorch-related modules only when the BoTorch optimizer is selected.
    """
    try:
        import torch
        import gpytorch
        from botorch.fit import fit_gpytorch_mll
        from botorch.models import SingleTaskGP
        from gpytorch.mlls import ExactMarginalLogLikelihood
    except ImportError as e:
        from .utils import apx_print
        apx_print(f"Error importing BoTorch dependencies: {e}")
        apx_print("Install torch, gpytorch, and botorch to use OPTIMIZER = botorch")
        return None

    return torch, gpytorch, fit_gpytorch_mll, SingleTaskGP, ExactMarginalLogLikelihood


def _select_botorch_device(torch, requested_device):
    """
    Select a torch device for BoTorch calculations.
    """
    requested_device = requested_device.lower()

    if requested_device == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        return torch.device("cpu")

    if requested_device.startswith("cuda"):
        if not torch.cuda.is_available():
            raise RuntimeError("BOTORCH_DEVICE requests CUDA, but torch.cuda.is_available() is False")
        return torch.device(requested_device)

    return torch.device(requested_device)


def _select_botorch_dtype(torch, requested_dtype):
    """
    Select a torch dtype for BoTorch calculations.
    """
    requested_dtype = requested_dtype.lower()

    if requested_dtype in ["float", "float32", "single"]:
        return torch.float32
    if requested_dtype in ["double", "float64"]:
        return torch.float64

    raise ValueError(f"Unsupported BOTORCH_DTYPE: {requested_dtype}")


def _scale_encoded_candidates_to_unit_cube(X, numpy_dtype=np.float32):
    """
    Scale encoded candidates to the [0, 1] unit cube.
    """
    X = np.asarray(X, dtype=numpy_dtype)
    if not X.flags.writeable:
        X = X.copy()

    feature_min = np.min(X, axis=0).astype(numpy_dtype)
    feature_max = np.max(X, axis=0).astype(numpy_dtype)
    feature_range = feature_max - feature_min
    feature_range[feature_range < 1.0e-12] = 1.0

    X -= feature_min
    X /= feature_range
    np.clip(X, 0.0, 1.0, out=X)

    return X, feature_min, feature_range


def _standardize_encoded_candidates(X, numpy_dtype=np.float32):
    """
    Standardize non-constant columns to zero mean and unit variance.
    """
    X = np.asarray(X, dtype=numpy_dtype)
    std = np.std(X, axis=0)
    keep = std > 0.0
    if not np.any(keep):
        return np.empty((X.shape[0], 0), dtype=numpy_dtype)
    transformed = (X[:, keep] - np.mean(X[:, keep], axis=0)) / std[keep]
    return np.ascontiguousarray(transformed, dtype=numpy_dtype)


def _transform_encoded_candidates_for_botorch(X, input_transform, numpy_dtype=np.float32):
    """
    Transform encoded candidates before BoTorch fitting/scoring.
    """
    input_transform = str(input_transform).strip().lower()

    if input_transform in ["standardize", "center", "centering", "standardize_x"]:
        return _standardize_encoded_candidates(X, numpy_dtype=numpy_dtype)

    if input_transform in ["unit_cube", "minmax", "normalize", "botorch"]:
        transformed, _, _ = _scale_encoded_candidates_to_unit_cube(X, numpy_dtype=numpy_dtype)
        return np.ascontiguousarray(transformed, dtype=numpy_dtype)

    if input_transform in ["none", "raw"]:
        return np.ascontiguousarray(np.asarray(X, dtype=numpy_dtype), dtype=numpy_dtype)

    raise ValueError(f"Unsupported BOTORCH_INPUT_TRANSFORM: {input_transform}")


def _load_discrete_candidate_features(numpy_dtype=np.float32):
    """
    Load raw atomic configurations from candidates.csv as integer site labels.
    """
    candidates_df = pd.read_csv("candidates.csv")
    if candidates_df.shape[1] < 2:
        raise ValueError("candidates.csv must contain structure_id and at least one site column")

    structure_ids = candidates_df.iloc[:, 0].to_numpy(dtype=int)
    expected_ids = np.arange(len(structure_ids), dtype=int)
    if not np.array_equal(structure_ids, expected_ids):
        sorted_order = np.argsort(structure_ids)
        sorted_ids = structure_ids[sorted_order]
        expected_sorted_ids = np.arange(int(sorted_ids[-1]) + 1, dtype=int)
        if not np.array_equal(sorted_ids, expected_sorted_ids):
            raise ValueError("Discrete BoTorch kernels require dense structure_id values starting at 0")
        candidates_df = candidates_df.iloc[sorted_order].reset_index(drop=True)

    site_values = candidates_df.iloc[:, 1:].astype(str)
    flat_values = site_values.to_numpy().reshape(-1)
    atom_types = sorted(pd.unique(flat_values).tolist())
    atom_to_code = {atom: idx for idx, atom in enumerate(atom_types)}
    X = site_values.apply(lambda column: column.map(atom_to_code)).to_numpy(dtype=numpy_dtype)

    return np.ascontiguousarray(X, dtype=numpy_dtype), atom_types


def _score_botorch_batch(torch, model, batch_X, score, best_f, xi, use_q_batch=True):
    """
    Score one candidate batch with BoTorch acquisition functions.
    """
    score = score.upper()
    threshold = best_f + torch.as_tensor(xi, device=batch_X.device, dtype=batch_X.dtype)

    if score not in ["TS", "EI"]:
        raise ValueError(f"Unsupported BoTorch SCORE: {score}")

    posterior = model.posterior(batch_X)
    mean = posterior.mean.reshape(-1)
    variance = posterior.variance.reshape(-1).clamp_min(0.0)
    std = variance.sqrt()

    if score == "TS":
        return mean + std * torch.randn_like(mean)

    improvement = mean - threshold
    normal = torch.distributions.Normal(
        torch.zeros((), device=batch_X.device, dtype=batch_X.dtype),
        torch.ones((), device=batch_X.device, dtype=batch_X.dtype),
    )
    min_std = torch.as_tensor(
        torch.finfo(batch_X.dtype).eps,
        device=batch_X.device,
        dtype=batch_X.dtype,
    )
    safe_std = std.clamp_min(min_std)
    u = improvement / safe_std
    ei = improvement * normal.cdf(u) + std * torch.exp(normal.log_prob(u))
    return torch.where(std > min_std, ei, improvement.clamp_min(0.0))


def _get_random_tie_batch_best(torch, acquisition, tie_atol=1.0e-12, tie_rtol=1.0e-12):
    """
    Return the best acquisition value and a random offset among tied maxima.
    """
    finite_mask = torch.isfinite(acquisition)
    if not finite_mask.any():
        return None, None, 0

    finite_values = acquisition[finite_mask]
    best_value_tensor = finite_values.max()
    scale = best_value_tensor.abs().clamp_min(1.0)
    atol = tie_atol * scale

    tie_mask = finite_mask & torch.isclose(
        acquisition,
        best_value_tensor,
        rtol=tie_rtol,
        atol=atol,
    )
    tie_offsets = torch.nonzero(tie_mask, as_tuple=False).reshape(-1)
    if tie_offsets.numel() == 0:
        return None, None, 0

    random_index = torch.randint(tie_offsets.numel(), (1,), device=acquisition.device)
    selected_offset = int(tie_offsets[random_index].item())
    return best_value_tensor.item(), selected_offset, int(tie_offsets.numel())


def _values_tied(a, b, tie_atol=1.0e-12, tie_rtol=1.0e-12):
    """
    Compare scalar acquisition values with a small tolerance for tie-breaking.
    """
    scale = max(abs(a), abs(b), 1.0)
    return abs(a - b) <= tie_atol * scale + tie_rtol * scale


def _score_on_the_fly_raw_codes(
    torch,
    model,
    raw_codes,
    atom_types,
    botorch_setting,
    score,
    best_f,
    xi,
    dtype,
    device,
    batch_size,
    numpy_dtype,
):
    from .on_the_fly import build_on_the_fly_kernel_features

    raw_codes = np.asarray(raw_codes, dtype=np.int64)
    values = np.empty(raw_codes.shape[0], dtype=np.float64)
    batch_size = max(1, int(batch_size))

    for start in range(0, raw_codes.shape[0], batch_size):
        end = min(start + batch_size, raw_codes.shape[0])
        batch_X_np, _ = build_on_the_fly_kernel_features(
            raw_codes[start:end],
            atom_types,
            use_local_env=botorch_setting["use_local_env"],
            local_env_type=botorch_setting["local_env_type"],
            numpy_dtype=numpy_dtype,
        )
        batch_X = torch.as_tensor(batch_X_np, device=device, dtype=dtype)
        acquisition = _score_botorch_batch(
            torch,
            model,
            batch_X,
            score,
            best_f,
            xi,
            use_q_batch=False,
        )
        acquisition = torch.nan_to_num(acquisition, nan=-torch.inf)
        values[start:end] = acquisition.detach().cpu().numpy()

    return values


def _raw_code_config_hash(raw_code_row, atom_types):
    from .on_the_fly import atomic_config_hash

    return atomic_config_hash(atom_types[int(code)] for code in raw_code_row)


def _penalize_existing_on_the_fly_configs(values, raw_codes, atom_types, existing_hashes):
    for index, raw_code_row in enumerate(raw_codes):
        if _raw_code_config_hash(raw_code_row, atom_types) in existing_hashes:
            values[index] = -np.inf
    return values


def _select_best_finite(values):
    finite_indices = np.flatnonzero(np.isfinite(values))
    if finite_indices.size == 0:
        return None
    local_index = int(np.argmax(values[finite_indices]))
    return int(finite_indices[local_index])


def _run_on_the_fly_sa_search(
    torch,
    gpytorch,
    model,
    botorch_setting,
    space,
    atom_types,
    existing_hashes,
    score,
    best_f,
    xi,
    dtype,
    device,
    numpy_dtype,
    current_sample_id,
):
    from .on_the_fly import (
        configs_to_raw_codes,
        generate_unique_random_config,
        propose_swap_neighbors,
        random_raw_code_configs,
    )

    rng = np.random.default_rng(int(current_sample_id))
    num_restarts = botorch_setting["on_the_fly_sa_restarts"]
    if num_restarts is None:
        num_restarts = botorch_setting["sa_chains"]
    num_steps = botorch_setting["on_the_fly_sa_steps"]
    if num_steps is None:
        num_steps = botorch_setting["sa_steps"]
    initial_temperature = botorch_setting["on_the_fly_sa_initial_temperature"]
    if initial_temperature is None:
        initial_temperature = botorch_setting["sa_initial_temperature"]
    cooling_rate = botorch_setting["on_the_fly_sa_cooling_rate"]
    final_temperature = botorch_setting["sa_final_temperature"]
    eval_batch_size = max(1, int(botorch_setting["sa_eval_batch_size"]))

    num_restarts = max(1, int(num_restarts))
    num_steps = max(0, int(num_steps))
    initial_temperature = max(float(initial_temperature), 1.0e-12)
    final_temperature = max(float(final_temperature), 1.0e-12)
    if cooling_rate is not None:
        cooling_rate = min(max(float(cooling_rate), 1.0e-12), 1.0)

    current_codes = random_raw_code_configs(
        space,
        atom_types,
        num_restarts,
        rng,
        existing_hashes=existing_hashes,
    )
    if current_codes.shape[0] == 0:
        fallback_config = generate_unique_random_config(space=space)
        return configs_to_raw_codes([fallback_config], atom_types, numpy_dtype=np.int64)[0], None

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        current_values = _score_on_the_fly_raw_codes(
            torch,
            model,
            current_codes,
            atom_types,
            botorch_setting,
            score,
            best_f,
            xi,
            dtype,
            device,
            eval_batch_size,
            numpy_dtype,
        )
        current_values = _penalize_existing_on_the_fly_configs(
            current_values,
            current_codes,
            atom_types,
            existing_hashes,
        )
        best_index = _select_best_finite(current_values)
        if best_index is None:
            fallback_config = generate_unique_random_config(space=space)
            return configs_to_raw_codes([fallback_config], atom_types, numpy_dtype=np.int64)[0], None

        best_codes = current_codes[best_index].copy()
        best_value = float(current_values[best_index])

        for step in range(num_steps):
            if cooling_rate is not None:
                temperature = initial_temperature * (cooling_rate ** step)
            elif num_steps <= 1:
                temperature = final_temperature
            else:
                fraction = step / float(num_steps - 1)
                temperature = initial_temperature * (final_temperature / initial_temperature) ** fraction
            temperature = max(float(temperature), 1.0e-12)

            proposal_codes = propose_swap_neighbors(rng, current_codes)
            proposal_values = _score_on_the_fly_raw_codes(
                torch,
                model,
                proposal_codes,
                atom_types,
                botorch_setting,
                score,
                best_f,
                xi,
                dtype,
                device,
                eval_batch_size,
                numpy_dtype,
            )
            proposal_values = _penalize_existing_on_the_fly_configs(
                proposal_values,
                proposal_codes,
                atom_types,
                existing_hashes,
            )

            finite_mask = np.isfinite(proposal_values)
            improvement = proposal_values - current_values
            accept_probability = np.exp(np.clip(improvement / temperature, -80.0, 0.0))
            accept_mask = finite_mask & (
                (improvement >= 0.0) | (rng.random(current_codes.shape[0]) < accept_probability)
            )
            current_codes[accept_mask] = proposal_codes[accept_mask]
            current_values[accept_mask] = proposal_values[accept_mask]

            proposal_best_index = _select_best_finite(proposal_values)
            if proposal_best_index is not None and proposal_values[proposal_best_index] > best_value:
                best_value = float(proposal_values[proposal_best_index])
                best_codes = proposal_codes[proposal_best_index].copy()

    if _raw_code_config_hash(best_codes, atom_types) in existing_hashes:
        fallback_config = generate_unique_random_config(space=space)
        return configs_to_raw_codes([fallback_config], atom_types, numpy_dtype=np.int64)[0], None

    return best_codes, best_value


def _run_on_the_fly_botorch_sampling(
    current_sample_id,
    torch,
    gpytorch,
    fit_gpytorch_mll,
    SingleTaskGP,
    ExactMarginalLogLikelihood,
    botorch_setting,
    device,
    dtype,
    numpy_dtype,
):
    from .utils import apx_print
    from .on_the_fly import (
        ON_THE_FLY_MODE,
        load_existing_atomic_config_hashes,
        load_on_the_fly_training_data,
        raw_codes_to_configs,
    )

    score = botorch_setting["score"].upper()
    gp_kernel = botorch_setting["gp_kernel"]
    maxiter = botorch_setting["maxiter"]
    xi = botorch_setting["xi"]
    standardize_y = botorch_setting["standardize_y"]
    use_discrete_kernel = is_discrete_botorch_kernel(gp_kernel)

    if not use_discrete_kernel:
        apx_print(
            "On-the-fly BoTorch sampling currently requires "
            "BOTORCH_GP_KERNEL = hamming, ot, or hamming_ot"
        )
        return None

    try:
        X, samples_df, metadata = load_on_the_fly_training_data(
            numpy_dtype=numpy_dtype,
            use_local_env=botorch_setting["use_local_env"],
            local_env_type=botorch_setting["local_env_type"],
        )
    except Exception as e:
        apx_print(f"Error loading on-the-fly BoTorch inputs: {e}")
        return None

    if len(samples_df) < 2:
        apx_print("On-the-fly BoTorch GP requires at least two completed samples in samples.csv")
        return None

    space = metadata["space"]
    atom_types = metadata["atom_types"]
    num_sites = int(metadata["num_sites"])
    kernel_context = dict(botorch_setting)
    kernel_context["use_local_env"] = bool(botorch_setting["use_local_env"])
    kernel_context["site_cost"] = build_site_graph_distance(num_sites)
    kernel_context["atom_types"] = atom_types
    if botorch_setting["use_local_env"]:
        kernel_context["local_env_type"] = metadata["local_env_type"]
        kernel_context["num_sites"] = num_sites
        kernel_context["env_dim"] = int(metadata["env_dim"])
        kernel_context["env_distance"] = botorch_setting["env_distance"]
        kernel_context["env_lengthscale"] = botorch_setting["env_lengthscale"]
        kernel_context["ot_env_mismatch_penalty"] = botorch_setting["ot_env_mismatch_penalty"]

    apx_print(
        f"On-the-fly BoTorch input: atoms={atom_types}, sites={num_sites}, "
        f"completed_samples={len(samples_df)}, kernel_dim={X.shape[1]}"
    )

    energy_numpy_dtype = np.float64 if dtype == torch.float64 else np.float32
    t_initial = -samples_df["energy"].to_numpy(dtype=energy_numpy_dtype)
    torch.manual_seed(int(current_sample_id))
    if device.type == "cuda":
        torch.cuda.manual_seed_all(int(current_sample_id))

    train_X = torch.as_tensor(X, device=device, dtype=dtype)
    train_Y_raw = torch.as_tensor(t_initial, device=device, dtype=dtype).unsqueeze(-1)
    train_Y_mean = train_Y_raw.mean()
    train_Y_std = train_Y_raw.std(unbiased=False).clamp_min(1.0e-12)
    train_Y = (train_Y_raw - train_Y_mean) / train_Y_std if standardize_y else train_Y_raw
    apx_print(
        f"On-the-fly BoTorch train_Y std: {float(train_Y_std):.6g}, "
        f"standardize_y: {standardize_y}"
    )

    try:
        model = create_botorch_model(
            torch,
            gpytorch,
            SingleTaskGP,
            train_X,
            train_Y,
            gp_kernel,
            kernel_context=kernel_context,
        )
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll, optimizer_kwargs={"options": {"maxiter": maxiter}})
        model.eval()
    except Exception as e:
        apx_print(f"Error fitting on-the-fly BoTorch GP: {e}")
        return None

    existing_hashes = load_existing_atomic_config_hashes(space.num_sites)
    best_f = train_Y.max()
    sa_score = str(botorch_setting["sa_score"] or score).strip().upper()
    if sa_score not in ["TS", "EI"]:
        apx_print(f"Unsupported BOTORCH_SA_SCORE: {sa_score}; using {score}")
        sa_score = score

    try:
        best_codes, best_value = _run_on_the_fly_sa_search(
            torch,
            gpytorch,
            model,
            botorch_setting,
            space,
            atom_types,
            existing_hashes,
            sa_score,
            best_f,
            xi,
            dtype,
            device,
            numpy_dtype,
            current_sample_id,
        )
    except Exception as e:
        apx_print(f"Error during on-the-fly BoTorch SA: {e}")
        return None

    atomic_config = raw_codes_to_configs(best_codes.reshape(1, -1), atom_types)[0]
    if best_value is None:
        apx_print("On-the-fly BoTorch SA fell back to a unique random configuration")
    else:
        apx_print(f"On-the-fly BoTorch SA score: {best_value:.6g}")

    if device.type == "cuda":
        torch.cuda.empty_cache()

    return {
        "mode": ON_THE_FLY_MODE,
        "sampling_method": "bayes_sa",
        "structure_id": int(current_sample_id),
        "atomic_config": atomic_config,
    }


def run_botorch_sampling(current_sample_id):
    """
    Execute BoTorch GP sampling function with torch device support.
    """
    from .utils import apx_print, read_botorch_setting, load_encoded_data_from_cache

    botorch_modules = _load_botorch_modules()
    if botorch_modules is None:
        return None

    torch, gpytorch, fit_gpytorch_mll, SingleTaskGP, ExactMarginalLogLikelihood = botorch_modules
    try:
        botorch_setting = read_botorch_setting()
    except Exception as e:
        apx_print(f"Error reading BoTorch settings: {e}")
        return None

    try:
        device = _select_botorch_device(torch, botorch_setting["device"])
        dtype = _select_botorch_dtype(torch, botorch_setting["dtype"])
    except Exception as e:
        apx_print(f"Error selecting BoTorch device/dtype: {e}")
        return None

    if device.type == "cuda":
        apx_print(f"BoTorch device: {device} ({torch.cuda.get_device_name(device)})")
    else:
        apx_print(f"BoTorch device: {device}")

    score = botorch_setting["score"].upper()
    batch_size = botorch_setting["batch_size"]
    maxiter = botorch_setting["maxiter"]
    xi = botorch_setting["xi"]
    input_transform = botorch_setting["input_transform"]
    gp_kernel = botorch_setting["gp_kernel"]
    standardize_y = botorch_setting["standardize_y"]
    use_discrete_kernel = is_discrete_botorch_kernel(gp_kernel)

    if score not in ["TS", "EI"]:
        apx_print(f"Unsupported BoTorch SCORE: {score}")
        return None

    if batch_size <= 0:
        apx_print("BOTORCH_BATCH_SIZE must be positive")
        return None

    numpy_dtype = np.float64 if dtype == torch.float64 else np.float32

    if not os.path.exists("candidates.csv"):
        return _run_on_the_fly_botorch_sampling(
            current_sample_id,
            torch,
            gpytorch,
            fit_gpytorch_mll,
            SingleTaskGP,
            ExactMarginalLogLikelihood,
            botorch_setting,
            device,
            dtype,
            numpy_dtype,
        )

    kernel_context = dict(botorch_setting)
    discrete_lookup_X = None

    if use_discrete_kernel:
        try:
            if botorch_setting["use_local_env"]:
                from .encoder import load_or_build_local_environment_features

                local_env_data = load_or_build_local_environment_features(
                    local_env_type=botorch_setting["local_env_type"],
                    cache_filename=botorch_setting["local_env_cache"],
                )
                X = np.ascontiguousarray(local_env_data["X_kernel"], dtype=numpy_dtype)
                discrete_lookup_X = np.ascontiguousarray(local_env_data["X_raw"], dtype=numpy_dtype)
                atom_types = local_env_data["atom_types"]
                num_sites = int(local_env_data["num_sites"])
                kernel_context["use_local_env"] = True
                kernel_context["local_env_type"] = local_env_data["local_env_type"]
                kernel_context["num_sites"] = num_sites
                kernel_context["env_dim"] = int(local_env_data["env_dim"])
                kernel_context["env_distance"] = botorch_setting["env_distance"]
                kernel_context["env_lengthscale"] = botorch_setting["env_lengthscale"]
                kernel_context["ot_env_mismatch_penalty"] = botorch_setting["ot_env_mismatch_penalty"]
            else:
                X, atom_types = _load_discrete_candidate_features(numpy_dtype=numpy_dtype)
                num_sites = X.shape[1]
                kernel_context["use_local_env"] = False

            kernel_context["site_cost"] = build_site_graph_distance(num_sites)
            kernel_context["atom_types"] = atom_types
        except Exception as e:
            apx_print(f"Error loading discrete BoTorch inputs: {e}")
            return None

        if botorch_setting["use_local_env"]:
            apx_print(
                f"BoTorch discrete kernel input: atoms={atom_types}, sites={num_sites}, "
                f"local_env={kernel_context['local_env_type']}, env_dim={kernel_context['env_dim']}, "
                f"kernel_dim={X.shape[1]}"
            )
        else:
            apx_print(
                f"BoTorch discrete kernel input: atoms={atom_types}, "
                f"dim/sites: {X.shape[1]}"
            )
    else:
        success, X = load_encoded_data_from_cache()
        if not success:
            return None

        try:
            X = _transform_encoded_candidates_for_botorch(
                X,
                input_transform=input_transform,
                numpy_dtype=numpy_dtype,
            )
        except Exception as e:
            apx_print(f"Error transforming BoTorch inputs: {e}")
            return None

        apx_print(
            f"BoTorch input transform: {input_transform}, "
            f"range: [{float(np.nanmin(X)):.6g}, {float(np.nanmax(X)):.6g}], "
            f"dim: {X.shape[1]}"
        )

    num_candidates = X.shape[0]
    if X.shape[1] == 0:
        apx_print("BoTorch input transform removed all feature columns")
        return None

    samples_df = pd.read_csv("samples.csv")
    if len(samples_df) < 2:
        apx_print("BoTorch GP requires at least two initial samples in samples.csv")
        return None

    calculated_ids = samples_df["structure_id"].values.astype(int)
    if np.any(calculated_ids < 0) or np.any(calculated_ids >= num_candidates):
        apx_print("samples.csv contains structure_id outside encoded candidate range")
        return None

    energy_numpy_dtype = np.float64 if dtype == torch.float64 else np.float32
    t_initial = -samples_df["energy"].to_numpy(dtype=energy_numpy_dtype)
    train_X_np = X[calculated_ids]

    torch.manual_seed(int(current_sample_id))
    if device.type == "cuda":
        torch.cuda.manual_seed_all(int(current_sample_id))

    train_X = torch.as_tensor(train_X_np, device=device, dtype=dtype)
    train_Y_raw = torch.as_tensor(t_initial, device=device, dtype=dtype).unsqueeze(-1)
    train_Y_mean = train_Y_raw.mean()
    train_Y_std = train_Y_raw.std(unbiased=False).clamp_min(1.0e-12)
    train_Y = (train_Y_raw - train_Y_mean) / train_Y_std if standardize_y else train_Y_raw
    apx_print(f"BoTorch train_X range: [{float(train_X_np.min()):.6g}, {float(train_X_np.max()):.6g}]")
    apx_print(
        f"BoTorch train_Y std: {float(train_Y_std):.6g}, "
        f"standardize_y: {standardize_y}"
    )

    try:
        model = create_botorch_model(
            torch,
            gpytorch,
            SingleTaskGP,
            train_X,
            train_Y,
            gp_kernel,
            kernel_context=kernel_context,
        )
        mll = ExactMarginalLogLikelihood(model.likelihood, model)
        fit_gpytorch_mll(mll, optimizer_kwargs={"options": {"maxiter": maxiter}})
        model.eval()
    except Exception as e:
        apx_print(f"Error fitting BoTorch GP: {e}")
        return None

    best_f = train_Y.max()
    calculated_id_set = set(int(idx) for idx in calculated_ids)
    best_value = None
    best_id = None
    best_tie_count = 0

    if botorch_setting["sa_screening"]:
        try:
            best_id, best_value = run_botorch_sa_screening(
                torch,
                gpytorch,
                model,
                X,
                calculated_id_set,
                botorch_setting,
                score,
                best_f,
                xi,
                dtype,
                device,
                use_discrete_kernel,
                current_sample_id,
                _score_botorch_batch,
                discrete_lookup_X=discrete_lookup_X,
            )
        except Exception as e:
            apx_print(f"Error during BoTorch SA screening: {e}")
            return None

        if best_id is None or best_value is None or not np.isfinite(best_value):
            apx_print("BoTorch SA screening failed to select a candidate")
            return None

        apx_print(
            f"BoTorch SA selected structure ID {best_id} "
            f"with screening score {best_value:.6g}"
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return best_id

    try:
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            for start in range(0, num_candidates, batch_size):
                end = min(start + batch_size, num_candidates)
                batch_np = np.ascontiguousarray(X[start:end])
                batch_X = torch.as_tensor(batch_np, device=device, dtype=dtype)
                offset_tensor = None
                acquisition = _score_botorch_batch(
                    torch,
                    model,
                    batch_X,
                    score,
                    best_f,
                    xi,
                    use_q_batch=not use_discrete_kernel,
                )
                acquisition = torch.nan_to_num(acquisition, nan=-torch.inf)

                sampled_offsets = [idx - start for idx in calculated_id_set if start <= idx < end]
                if sampled_offsets:
                    offset_tensor = torch.as_tensor(sampled_offsets, device=device, dtype=torch.long)
                    acquisition[offset_tensor] = -torch.inf

                batch_best_value, batch_best_offset, batch_tie_count = _get_random_tie_batch_best(torch, acquisition)
                if batch_best_value is not None:
                    if best_value is None:
                        best_value = batch_best_value
                        best_id = start + int(batch_best_offset)
                        best_tie_count = batch_tie_count
                    elif batch_best_value > best_value and not _values_tied(batch_best_value, best_value):
                        best_value = batch_best_value
                        best_id = start + int(batch_best_offset)
                        best_tie_count = batch_tie_count
                    elif _values_tied(batch_best_value, best_value):
                        total_tie_count = best_tie_count + batch_tie_count
                        choose_new = torch.randint(total_tie_count, (1,), device=device).item() >= best_tie_count
                        if choose_new:
                            best_value = batch_best_value
                            best_id = start + int(batch_best_offset)
                        best_tie_count = total_tie_count

                del acquisition, batch_X, batch_np
                if offset_tensor is not None:
                    del offset_tensor
                if device.type == "cuda":
                    torch.cuda.empty_cache()
    except Exception as e:
        apx_print(f"Error scoring candidates with BoTorch GP: {e}")
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return None

    if best_id is None or not np.isfinite(best_value):
        apx_print("BoTorch sampling failed to select a candidate")
        return None

    apx_print(
        f"BoTorch SCORE: {score}, batch size: {batch_size}, "
        f"maxiter: {maxiter}, kernel: {gp_kernel}"
    )

    if device.type == "cuda":
        torch.cuda.empty_cache()

    return best_id
