"""
On-the-fly candidate generation for PyAPX.

This module supports runs without candidates.csv by treating the
ATOMIC_POSITIONS entries in qe_template.in as configurable sites.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import os
import random
from typing import Optional

import numpy as np
import pandas as pd


ENUMERATED_CANDIDATES_MODE = "enumerated_candidates"
ON_THE_FLY_MODE = "on_the_fly"

_PLACEHOLDER_ATOMS = {"X"}
_QE_POSITION_STOP_CARDS = {
    "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS",
    "K_POINTS",
    "CELL_PARAMETERS",
    "CONSTRAINTS",
    "OCCUPATIONS",
    "ATOMIC_FORCES",
    "HUBBARD",
    "SOLVENTS",
}


@dataclass(frozen=True)
class AtomicPositionsTemplate:
    lines: list[str]
    positions_start: int
    positions_end: int
    atom_labels: list[str]
    position_tails: list[str]
    variable_indices: list[int]

    @property
    def num_variable_sites(self) -> int:
        return len(self.variable_indices)


@dataclass(frozen=True)
class OnTheFlySpace:
    template: AtomicPositionsTemplate
    composition: dict[str, int]
    atom_types: list[str]
    composition_labels: list[str]

    @property
    def num_sites(self) -> int:
        return self.template.num_variable_sites


def has_candidates_csv(path: str = "candidates.csv") -> bool:
    return os.path.exists(path)


def get_candidate_mode() -> str:
    return ENUMERATED_CANDIDATES_MODE if has_candidates_csv() else ON_THE_FLY_MODE


def _parse_bool_setting(value, default=False):
    if value is None:
        return default
    return str(value).strip().lower() in {"true", "t", "yes", "y", "1", "on"}


def use_on_the_fly_if_no_candidates() -> bool:
    from .utils import read_card_value

    value = read_card_value("apx.in", "USE_ON_THE_FLY_IF_NO_CANDIDATES")
    return _parse_bool_setting(value, default=True)


def read_max_duplicate_trials(default: int = 10000) -> int:
    from .utils import apx_print, read_card_value

    for key in ("MAX_DUPLICATE_TRIALS", "ON_THE_FLY_MAX_DUPLICATE_TRIALS"):
        value = read_card_value("apx.in", key)
        if value is None or str(value).strip() == "":
            continue
        try:
            return max(1, int(value))
        except ValueError:
            apx_print(f"Invalid {key}: {value}; using default {default}")
            return default
    return default


def read_random_seed(default=None):
    from .utils import apx_print, read_card_value

    value = read_card_value("apx.in", "RANDOM_SEED")
    if value is None or str(value).strip() == "":
        return default
    try:
        return int(value)
    except ValueError:
        apx_print(f"Invalid RANDOM_SEED: {value}; using default {default}")
        return default


def _is_qe_position_stop_line(stripped_line: str) -> bool:
    if not stripped_line:
        return True
    if stripped_line.startswith("&"):
        return True
    upper = stripped_line.split(None, 1)[0].upper()
    return upper in _QE_POSITION_STOP_CARDS


def read_atomic_positions_template(qe_template_path: str = "qe_template.in") -> AtomicPositionsTemplate:
    if not os.path.exists(qe_template_path):
        raise FileNotFoundError(
            "candidates.csv is not found and qe_template.in is also missing."
        )

    with open(qe_template_path, "r") as template_file:
        lines = template_file.readlines()

    positions_start = None
    for index, line in enumerate(lines):
        if line.strip().upper().startswith("ATOMIC_POSITIONS"):
            positions_start = index + 1
            break

    if positions_start is None:
        raise ValueError("ATOMIC_POSITIONS section was not found in qe_template.in.")

    atom_labels = []
    position_tails = []
    positions_end = positions_start
    for index, line in enumerate(lines[positions_start:], start=positions_start):
        stripped = line.strip()
        if _is_qe_position_stop_line(stripped):
            positions_end = index
            break
        parts = stripped.split()
        if len(parts) < 4:
            raise ValueError(f"Invalid ATOMIC_POSITIONS line: {line.rstrip()}")
        atom_labels.append(parts[0])
        position_tails.append("   ".join(parts[1:]))
    else:
        positions_end = len(lines)

    if not atom_labels:
        raise ValueError("ATOMIC_POSITIONS section does not contain any atoms.")

    x_indices = [
        index
        for index, atom in enumerate(atom_labels)
        if atom.strip().upper() in _PLACEHOLDER_ATOMS
    ]
    variable_indices = x_indices if x_indices else list(range(len(atom_labels)))

    return AtomicPositionsTemplate(
        lines=lines,
        positions_start=positions_start,
        positions_end=positions_end,
        atom_labels=atom_labels,
        position_tails=position_tails,
        variable_indices=variable_indices,
    )


def _parse_composition_text(raw_value: str) -> dict[str, int]:
    composition = {}
    normalized = raw_value.replace(";", ",")
    for item in normalized.split(","):
        text = item.strip()
        if not text:
            continue
        if ":" in text:
            atom, count = text.split(":", 1)
        elif "=" in text:
            atom, count = text.split("=", 1)
        else:
            parts = text.split()
            if len(parts) != 2:
                raise ValueError(
                    "ON_THE_FLY_COMPOSITION entries must look like B:18,C:18,N:18"
                )
            atom, count = parts
        atom = atom.strip()
        if not atom:
            raise ValueError("ON_THE_FLY_COMPOSITION contains an empty atom label")
        count_int = int(str(count).strip())
        if count_int <= 0:
            raise ValueError("ON_THE_FLY_COMPOSITION counts must be positive")
        composition[atom] = composition.get(atom, 0) + count_int
    if not composition:
        raise ValueError("ON_THE_FLY_COMPOSITION did not define any atoms")
    return composition


def _read_explicit_composition() -> Optional[dict[str, int]]:
    from .utils import read_card_value

    for key in ("ON_THE_FLY_COMPOSITION", "ATOMIC_COMPOSITION"):
        value = read_card_value("apx.in", key)
        if value is not None and str(value).strip() != "":
            return _parse_composition_text(value)
    return None


def infer_sample_site_columns(samples_df: pd.DataFrame) -> list[str]:
    if samples_df.empty:
        return []
    columns = list(samples_df.columns)
    excluded = {"sample_id", "structure_id", "energy"}
    return [column for column in columns if column not in excluded]


def _composition_from_samples(num_sites: int, samples_path: str = "samples.csv") -> Optional[dict[str, int]]:
    if not os.path.exists(samples_path):
        return None
    try:
        samples_df = pd.read_csv(samples_path)
    except Exception:
        return None
    if samples_df.empty:
        return None
    site_columns = infer_sample_site_columns(samples_df)
    if len(site_columns) != num_sites:
        return None
    first_config = samples_df.iloc[0][site_columns].astype(str).tolist()
    return dict(Counter(first_config))


def load_on_the_fly_space(qe_template_path: str = "qe_template.in") -> OnTheFlySpace:
    template = read_atomic_positions_template(qe_template_path)
    variable_atoms = [template.atom_labels[index] for index in template.variable_indices]
    placeholders = [
        atom for atom in variable_atoms if atom.strip().upper() in _PLACEHOLDER_ATOMS
    ]

    explicit_composition = _read_explicit_composition()
    if explicit_composition is not None:
        composition = explicit_composition
    elif placeholders:
        composition = _composition_from_samples(template.num_variable_sites)
        if composition is None:
            raise ValueError(
                "qe_template.in uses X placeholders, so on-the-fly composition "
                "cannot be inferred. Replace X with the desired initial composition "
                "or set ON_THE_FLY_COMPOSITION = B:18,C:18,N:18 in apx.in."
            )
    else:
        composition = dict(Counter(variable_atoms))

    total_atoms = sum(composition.values())
    if total_atoms != template.num_variable_sites:
        raise ValueError(
            f"On-the-fly composition contains {total_atoms} atoms, "
            f"but qe_template.in defines {template.num_variable_sites} variable sites"
        )

    atom_types = sorted(composition)
    composition_labels = []
    for atom in atom_types:
        composition_labels.extend([atom] * int(composition[atom]))

    return OnTheFlySpace(
        template=template,
        composition=composition,
        atom_types=atom_types,
        composition_labels=composition_labels,
    )


def atomic_config_hash(atomic_config) -> str:
    return "|".join(str(atom) for atom in atomic_config)


def load_existing_atomic_config_hashes(
    num_sites: Optional[int] = None,
    samples_path: str = "samples.csv",
) -> set[str]:
    if not os.path.exists(samples_path):
        return set()
    try:
        samples_df = pd.read_csv(samples_path)
    except Exception:
        return set()
    if samples_df.empty:
        return set()

    site_columns = infer_sample_site_columns(samples_df)
    if num_sites is not None and len(site_columns) != int(num_sites):
        return set()
    existing = set()
    for _, row in samples_df[site_columns].dropna().iterrows():
        config = row.astype(str).tolist()
        existing.add(atomic_config_hash(config))
    return existing


def generate_unique_random_config(
    space: Optional[OnTheFlySpace] = None,
    rng: Optional[random.Random] = None,
    max_duplicate_trials: Optional[int] = None,
    existing_hashes: Optional[set[str]] = None,
) -> list[str]:
    if space is None:
        space = load_on_the_fly_space()
    if rng is None:
        rng = random.Random()
    if max_duplicate_trials is None:
        max_duplicate_trials = read_max_duplicate_trials()
    if existing_hashes is None:
        existing_hashes = load_existing_atomic_config_hashes(space.num_sites)

    base_config = list(space.composition_labels)
    for _ in range(max_duplicate_trials):
        config = list(base_config)
        rng.shuffle(config)
        if atomic_config_hash(config) not in existing_hashes:
            return config

    raise RuntimeError("Failed to generate a new unique random configuration.")


def sample_random_on_the_fly_structure(current_sample_id: Optional[int] = None) -> dict:
    space = load_on_the_fly_space()
    base_seed = read_random_seed(default=None)
    if base_seed is None:
        seed = None if current_sample_id is None else int(current_sample_id)
    else:
        seed = int(base_seed) + (0 if current_sample_id is None else int(current_sample_id))
    rng = random.Random(seed)
    atomic_config = generate_unique_random_config(space=space, rng=rng)
    structure_id = int(current_sample_id) if current_sample_id is not None else 0
    return {
        "mode": ON_THE_FLY_MODE,
        "sampling_method": "random",
        "structure_id": structure_id,
        "atomic_config": atomic_config,
    }


def write_qe_input_from_atomic_config(
    atomic_config,
    input_file: str,
    qe_template_path: str = "qe_template.in",
    output_dir: str = "dft_calc",
) -> list[str]:
    template = read_atomic_positions_template(qe_template_path)
    atomic_config = [str(atom) for atom in atomic_config]
    if len(atomic_config) != template.num_variable_sites:
        raise ValueError(
            f"atomic_config has {len(atomic_config)} atoms, "
            f"but qe_template.in defines {template.num_variable_sites} variable sites"
        )

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, input_file)

    output_atoms = list(template.atom_labels)
    for atom, position_index in zip(atomic_config, template.variable_indices):
        output_atoms[position_index] = atom

    with open(output_path, "w") as output_file_obj:
        output_file_obj.writelines(template.lines[: template.positions_start])
        for atom, tail in zip(output_atoms, template.position_tails):
            output_file_obj.write(f"{atom}   {tail}\n")
        if template.positions_end < len(template.lines):
            output_file_obj.writelines(template.lines[template.positions_end :])

    return atomic_config


def configs_to_raw_codes(configs, atom_types: list[str], numpy_dtype=np.float32):
    atom_to_code = {atom: index for index, atom in enumerate(atom_types)}
    raw_codes = np.empty((len(configs), len(configs[0]) if configs else 0), dtype=np.int64)
    for row_index, config in enumerate(configs):
        for column_index, atom in enumerate(config):
            atom = str(atom)
            if atom not in atom_to_code:
                raise ValueError(f"Unknown atom label in on-the-fly config: {atom}")
            raw_codes[row_index, column_index] = atom_to_code[atom]
    return np.ascontiguousarray(raw_codes, dtype=numpy_dtype)


def raw_codes_to_configs(raw_codes, atom_types: list[str]) -> list[list[str]]:
    raw_codes = np.asarray(raw_codes, dtype=np.int64)
    return [[atom_types[int(code)] for code in row] for row in raw_codes]


def build_on_the_fly_kernel_features(
    raw_codes,
    atom_types: list[str],
    use_local_env: bool = False,
    local_env_type: str = "NA",
    numpy_dtype=np.float32,
):
    raw_codes_int = np.ascontiguousarray(np.rint(raw_codes).astype(np.int64, copy=False))
    X_raw = np.ascontiguousarray(raw_codes_int, dtype=numpy_dtype)
    metadata = {
        "atom_types": list(atom_types),
        "num_sites": int(raw_codes_int.shape[1]),
        "use_local_env": bool(use_local_env),
    }

    if not use_local_env:
        return X_raw, metadata

    from .encoder import (
        _append_namod_sigma,
        _build_na_site_descriptors,
        _canonical_local_env_type,
        _read_neighbor_sites_for_local_env,
    )

    local_env_type = _canonical_local_env_type(local_env_type)
    neighbor_sites = _read_neighbor_sites_for_local_env(raw_codes_int.shape[1])
    X_env = _build_na_site_descriptors(raw_codes_int, atom_types, neighbor_sites)
    if local_env_type == "NAmod":
        X_env = _append_namod_sigma(X_env, neighbor_sites)

    num_configs, num_sites = raw_codes_int.shape
    env_dim = int(X_env.shape[2])
    X_kernel = np.concatenate(
        [X_raw, X_env.reshape(num_configs, num_sites * env_dim)],
        axis=1,
    )
    metadata.update(
        {
            "local_env_type": local_env_type,
            "neighbor_sites": neighbor_sites,
            "env_dim": env_dim,
        }
    )
    return np.ascontiguousarray(X_kernel, dtype=numpy_dtype), metadata


def load_on_the_fly_training_data(numpy_dtype=np.float32, use_local_env=False, local_env_type="NA"):
    space = load_on_the_fly_space()
    if not os.path.exists("samples.csv"):
        raise FileNotFoundError("samples.csv is required for on-the-fly BoTorch sampling")

    samples_df = pd.read_csv("samples.csv")
    if samples_df.empty:
        raise ValueError("samples.csv does not contain any completed samples")
    site_columns = infer_sample_site_columns(samples_df)
    if len(site_columns) != space.num_sites:
        raise ValueError(
            f"samples.csv has {len(site_columns)} site columns, "
            f"but qe_template.in defines {space.num_sites} variable sites"
        )

    configs = samples_df[site_columns].astype(str).values.tolist()
    raw_codes = configs_to_raw_codes(configs, space.atom_types, numpy_dtype=numpy_dtype)
    X, metadata = build_on_the_fly_kernel_features(
        raw_codes,
        space.atom_types,
        use_local_env=use_local_env,
        local_env_type=local_env_type,
        numpy_dtype=numpy_dtype,
    )
    metadata["space"] = space
    metadata["site_columns"] = site_columns
    metadata["raw_codes"] = np.ascontiguousarray(raw_codes, dtype=numpy_dtype)
    return X, samples_df, metadata


def random_raw_code_configs(
    space: OnTheFlySpace,
    atom_types: list[str],
    count: int,
    rng: np.random.Generator,
    existing_hashes: Optional[set[str]] = None,
    max_draws: Optional[int] = None,
) -> np.ndarray:
    count = max(0, int(count))
    if count == 0:
        return np.empty((0, space.num_sites), dtype=np.int64)

    atom_to_code = {atom: index for index, atom in enumerate(atom_types)}
    base_codes = np.asarray(
        [atom_to_code[atom] for atom in space.composition_labels],
        dtype=np.int64,
    )
    if existing_hashes is None:
        existing_hashes = set()
    if max_draws is None:
        max_draws = max(count * 20, read_max_duplicate_trials())

    configs = []
    seen = set()
    attempts = 0
    while len(configs) < count and attempts < max_draws:
        attempts += 1
        codes = rng.permutation(base_codes)
        config = tuple(atom_types[int(code)] for code in codes)
        config_hash = atomic_config_hash(config)
        if config_hash in existing_hashes or config_hash in seen:
            continue
        configs.append(codes)
        seen.add(config_hash)

    return np.ascontiguousarray(configs, dtype=np.int64)


def propose_swap_neighbors(rng: np.random.Generator, current_codes: np.ndarray) -> np.ndarray:
    proposals = np.array(current_codes, copy=True)
    num_sites = proposals.shape[1]
    for row_index in range(proposals.shape[0]):
        for _ in range(24):
            site_a, site_b = rng.choice(num_sites, size=2, replace=False)
            if proposals[row_index, site_a] == proposals[row_index, site_b]:
                continue
            proposals[row_index, site_a], proposals[row_index, site_b] = (
                proposals[row_index, site_b],
                proposals[row_index, site_a],
            )
            break
    return proposals
