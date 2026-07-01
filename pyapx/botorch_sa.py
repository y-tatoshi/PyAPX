"""
Simulated annealing screening helpers for BoTorch sampling.
"""

import numpy as np


def _score_botorch_candidate_ids(
    torch,
    model,
    X,
    candidate_ids,
    score,
    best_f,
    xi,
    dtype,
    device,
    batch_size,
    score_batch_fn,
    use_q_batch=True,
):
    """
    Score selected candidate IDs with the fitted BoTorch GP.
    """
    candidate_ids = np.asarray(candidate_ids, dtype=int)
    values = np.empty(candidate_ids.shape[0], dtype=np.float64)
    batch_size = max(1, int(batch_size))

    for start in range(0, candidate_ids.shape[0], batch_size):
        end = min(start + batch_size, candidate_ids.shape[0])
        batch_np = np.ascontiguousarray(X[candidate_ids[start:end]])
        batch_X = torch.as_tensor(batch_np, device=device, dtype=dtype)
        acquisition = score_batch_fn(
            torch,
            model,
            batch_X,
            score,
            best_f,
            xi,
            use_q_batch=use_q_batch,
        )
        acquisition = torch.nan_to_num(acquisition, nan=-torch.inf)
        values[start:end] = acquisition.detach().cpu().numpy()

    return values


def _sample_uncalculated_ids(rng, num_candidates, calculated_id_set, count):
    """
    Sample structure IDs that are not already calculated.
    """
    count = max(0, int(count))
    if count == 0:
        return np.empty(0, dtype=int)

    sampled = []
    seen = set()
    max_draws = max(count * 4, 1024)

    while len(sampled) < count and max_draws > 0:
        draw_size = min(max(count - len(sampled), 1024), max_draws)
        draws = rng.integers(0, num_candidates, size=draw_size, dtype=np.int64)
        for candidate_id in draws:
            candidate_id = int(candidate_id)
            if candidate_id in calculated_id_set or candidate_id in seen:
                continue
            sampled.append(candidate_id)
            seen.add(candidate_id)
            if len(sampled) >= count:
                break
        max_draws -= draw_size

    if len(sampled) < count:
        remaining = np.setdiff1d(
            np.arange(num_candidates, dtype=int),
            np.fromiter(calculated_id_set.union(seen), dtype=int),
            assume_unique=False,
        )
        if remaining.size > 0:
            rng.shuffle(remaining)
            sampled.extend(int(value) for value in remaining[: count - len(sampled)])

    return np.asarray(sampled, dtype=int)


def _build_discrete_config_lookup(X):
    """
    Build a raw-configuration-to-structure-id lookup for swap-neighbor SA.
    """
    codes = np.ascontiguousarray(np.rint(X).astype(np.int16, copy=False))
    lookup = {codes[index].tobytes(): int(index) for index in range(codes.shape[0])}
    return codes, lookup


def _propose_sa_candidate_ids(
    rng,
    current_ids,
    num_candidates,
    calculated_id_set,
    random_fraction,
    discrete_codes=None,
    config_lookup=None,
):
    """
    Propose new SA candidate IDs using swap-neighbor moves when possible.
    """
    proposals = np.empty_like(current_ids)
    random_fraction = min(max(float(random_fraction), 0.0), 1.0)
    use_swaps = discrete_codes is not None and config_lookup is not None

    for index, current_id in enumerate(current_ids):
        use_random_jump = (not use_swaps) or (rng.random() < random_fraction)
        proposal_id = None

        if not use_random_jump:
            current_config = discrete_codes[int(current_id)]
            num_sites = current_config.shape[0]
            for _ in range(12):
                site_a, site_b = rng.choice(num_sites, size=2, replace=False)
                if current_config[site_a] == current_config[site_b]:
                    continue
                proposed_config = current_config.copy()
                proposed_config[site_a], proposed_config[site_b] = proposed_config[site_b], proposed_config[site_a]
                found_id = config_lookup.get(proposed_config.tobytes())
                if found_id is not None and found_id not in calculated_id_set:
                    proposal_id = int(found_id)
                    break

        if proposal_id is None:
            random_ids = _sample_uncalculated_ids(rng, num_candidates, calculated_id_set, 1)
            proposal_id = int(random_ids[0]) if random_ids.size > 0 else int(current_id)

        proposals[index] = proposal_id

    return proposals


def run_botorch_sa_screening(
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
    score_batch_fn,
    discrete_lookup_X=None,
):
    """
    Use multi-chain simulated annealing to pre-screen candidates by acquisition.
    """
    from .utils import apx_print

    num_candidates = X.shape[0]
    num_chains = max(1, int(botorch_setting["sa_chains"]))
    num_steps = max(0, int(botorch_setting["sa_steps"]))
    initial_pool_size = max(num_chains, int(botorch_setting["sa_initial_pool_size"]))
    eval_batch_size = max(1, int(botorch_setting["sa_eval_batch_size"]))
    initial_temperature = max(float(botorch_setting["sa_initial_temperature"]), 1.0e-12)
    final_temperature = max(float(botorch_setting["sa_final_temperature"]), 1.0e-12)
    random_fraction = float(botorch_setting["sa_random_fraction"])
    sa_score = str(botorch_setting["sa_score"] or score).strip().upper()

    if sa_score not in ["TS", "EI"]:
        apx_print(f"Unsupported BOTORCH_SA_SCORE: {sa_score}; using {score}")
        sa_score = score

    rng = np.random.default_rng(int(current_sample_id))
    initial_ids = _sample_uncalculated_ids(
        rng,
        num_candidates,
        calculated_id_set,
        initial_pool_size,
    )
    if initial_ids.size == 0:
        apx_print("SA screening found no uncalculated candidates")
        return None, None

    if initial_ids.size < num_chains:
        num_chains = int(initial_ids.size)

    discrete_codes = None
    config_lookup = None
    if use_discrete_kernel and botorch_setting["sa_swap_neighbors"]:
        try:
            lookup_X = X if discrete_lookup_X is None else discrete_lookup_X
            discrete_codes, config_lookup = _build_discrete_config_lookup(lookup_X)
        except MemoryError:
            apx_print("SA swap-neighbor lookup exceeded memory; falling back to random jumps")
        except Exception as e:
            apx_print(f"SA swap-neighbor lookup failed: {e}; falling back to random jumps")

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        initial_values = _score_botorch_candidate_ids(
            torch,
            model,
            X,
            initial_ids,
            sa_score,
            best_f,
            xi,
            dtype,
            device,
            eval_batch_size,
            score_batch_fn,
            use_q_batch=not use_discrete_kernel,
        )

        order = np.argsort(initial_values)[::-1]
        current_ids = initial_ids[order[:num_chains]].copy()
        current_values = initial_values[order[:num_chains]].copy()
        best_id = int(current_ids[0])
        best_value = float(current_values[0])

        for step in range(num_steps):
            if num_steps <= 1:
                temperature = final_temperature
            else:
                fraction = step / float(num_steps - 1)
                temperature = initial_temperature * (final_temperature / initial_temperature) ** fraction

            proposal_ids = _propose_sa_candidate_ids(
                rng,
                current_ids,
                num_candidates,
                calculated_id_set,
                random_fraction,
                discrete_codes=discrete_codes,
                config_lookup=config_lookup,
            )
            proposal_values = _score_botorch_candidate_ids(
                torch,
                model,
                X,
                proposal_ids,
                sa_score,
                best_f,
                xi,
                dtype,
                device,
                eval_batch_size,
                score_batch_fn,
                use_q_batch=not use_discrete_kernel,
            )

            finite_mask = np.isfinite(proposal_values)
            improvement = proposal_values - current_values
            accept_probability = np.exp(np.clip(improvement / temperature, -80.0, 0.0))
            accept_mask = finite_mask & ((improvement >= 0.0) | (rng.random(num_chains) < accept_probability))
            current_ids[accept_mask] = proposal_ids[accept_mask]
            current_values[accept_mask] = proposal_values[accept_mask]

            proposal_best_index = int(np.nanargmax(proposal_values))
            if np.isfinite(proposal_values[proposal_best_index]) and proposal_values[proposal_best_index] > best_value:
                best_value = float(proposal_values[proposal_best_index])
                best_id = int(proposal_ids[proposal_best_index])

    return best_id, best_value
