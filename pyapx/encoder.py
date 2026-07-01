"""
PyAPX - Encoding Module
Provides feature encoding functionality for atomic configurations
"""

import numpy as np
import pandas as pd
import os
import pickle
from .utils import read_card, read_card_value

LOCAL_ENV_SCHEMA_VERSION = 2

# Conditional imports for dimension reduction
SKLEARN_AVAILABLE = False
TORCH_AVAILABLE = False

try:
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# PyTorch imports for Auto Encoder
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import DataLoader, TensorDataset
    TORCH_AVAILABLE = True
    
    class AutoEncoder(nn.Module):
        """
        PyTorch Auto Encoder for dimension reduction
        """
        def __init__(self, input_dim, latent_dim, hidden_layers):
            super(AutoEncoder, self).__init__()
            
            # Build encoder
            encoder_layers = []
            prev_dim = input_dim
            for hidden_dim in hidden_layers:
                encoder_layers.extend([
                    nn.Linear(prev_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(0.1)
                ])
                prev_dim = hidden_dim
            
            # Final encoder layer (to latent space)
            encoder_layers.append(nn.Linear(prev_dim, latent_dim))
            
            # Build decoder
            decoder_layers = []
            prev_dim = latent_dim
            for hidden_dim in reversed(hidden_layers):
                decoder_layers.extend([
                    nn.Linear(prev_dim, hidden_dim),
                    nn.ReLU(),
                    nn.Dropout(0.1)
                ])
                prev_dim = hidden_dim
            
            # Final decoder layer (back to input dimension)
            decoder_layers.append(nn.Linear(prev_dim, input_dim))
            
            self.encoder = nn.Sequential(*encoder_layers[:-1])  # Exclude final linear layer
            self.latent_layer = encoder_layers[-1]  # Final encoder layer
            self.decoder = nn.Sequential(*decoder_layers)
            
        def forward(self, x):
            # Encode
            encoded = self.encoder(x)
            latent = self.latent_layer(encoded)
            
            # Decode
            decoded = self.decoder(latent)
            
            return decoded, latent
        
        def encode(self, x):
            """Extract latent representation"""
            encoded = self.encoder(x)
            latent = self.latent_layer(encoded)
            return latent
            
except ImportError:
    TORCH_AVAILABLE = False


def _canonical_local_env_type(local_env_type):
    local_env_key = str(local_env_type).strip().lower()
    local_env_map = {
        "na": "NA",
        "namod": "NAmod",
    }
    if local_env_key not in local_env_map:
        raise ValueError(
            f"Unsupported BOTORCH_LOCAL_ENV_TYPE: {local_env_type}. "
            "Choose NA or NAmod."
        )
    return local_env_map[local_env_key]


def _load_candidate_site_values():
    if not os.path.exists("candidates.csv"):
        raise FileNotFoundError("candidates.csv not found; local environment features require candidates.csv")

    df = pd.read_csv("candidates.csv")
    if df.shape[1] < 2:
        raise ValueError("candidates.csv must contain structure_id and at least one site column")

    candidate_ids = df.iloc[:, 0].to_numpy(dtype=int)
    expected_ids = np.arange(len(candidate_ids), dtype=int)
    if not np.array_equal(candidate_ids, expected_ids):
        sorted_order = np.argsort(candidate_ids)
        sorted_ids = candidate_ids[sorted_order]
        expected_sorted_ids = np.arange(int(sorted_ids[-1]) + 1, dtype=int)
        if not np.array_equal(sorted_ids, expected_sorted_ids):
            raise ValueError("Local environment features require dense structure_id values starting at 0")
        df = df.iloc[sorted_order].reset_index(drop=True)
        candidate_ids = df.iloc[:, 0].to_numpy(dtype=int)

    site_df = df.iloc[:, 1:].astype(str)
    site_columns = list(site_df.columns)
    site_values = site_df.to_numpy()
    flat_values = site_values.reshape(-1)
    atom_types = sorted(pd.unique(flat_values).tolist())
    atom_to_code = {atom_type: index for index, atom_type in enumerate(atom_types)}
    raw_codes = site_df.apply(lambda column: column.map(atom_to_code)).to_numpy(dtype=np.int64)

    return raw_codes, atom_types, candidate_ids, site_columns, site_values


def _read_local_env_weight_metadata():
    weight_value = read_card_value("apx.in", "WEIGHT")
    if weight_value is None or str(weight_value).strip() == "":
        return {"WEIGHT": 0.0}

    try:
        weight = float(weight_value)
    except ValueError as exc:
        raise ValueError(f"WEIGHT must be numeric when saved to local env metadata: {weight_value}") from exc

    return {"WEIGHT": weight}


def _read_neighbor_sites_for_local_env(num_sites):
    neighbor_sites_data = read_card("apx.in", "NEIGHBOR_SITES")
    if not neighbor_sites_data:
        raise ValueError("NEIGHBOR_SITES is required when BOTORCH_USE_LOCAL_ENV=True")
    if len(neighbor_sites_data) != num_sites:
        raise ValueError(
            f"NEIGHBOR_SITES must contain {num_sites} lines, "
            f"but {len(neighbor_sites_data)} lines were found"
        )

    neighbor_sites = []
    for site_index, raw_line in enumerate(neighbor_sites_data):
        line = raw_line.split("#", 1)[0].strip()
        try:
            neighbors = [int(value) for value in line.split()]
        except ValueError as exc:
            raise ValueError(
                f"Invalid NEIGHBOR_SITES entry for site {site_index + 1}: {raw_line}"
            ) from exc

        for neighbor_site in neighbors:
            if neighbor_site < 1 or neighbor_site > num_sites:
                raise ValueError(
                    f"Invalid neighbor site {neighbor_site} for site {site_index + 1}; "
                    f"valid range is 1 to {num_sites}"
                )
        neighbor_sites.append(neighbors)

    return neighbor_sites


def _build_na_site_descriptors(raw_codes, atom_types, neighbor_sites):
    num_candidates, num_sites = raw_codes.shape
    num_atom_types = len(atom_types)
    max_neighbors = max((len(neighbors) for neighbors in neighbor_sites), default=0)
    env_dim = (max_neighbors + 1) * num_atom_types

    atom_eye = np.eye(num_atom_types, dtype=np.float32)
    site_onehot = atom_eye[raw_codes]
    descriptors = np.zeros((num_candidates, num_sites, env_dim), dtype=np.float32)
    descriptors[:, :, :num_atom_types] = site_onehot

    for site_index, neighbors in enumerate(neighbor_sites):
        for neighbor_slot, neighbor_site in enumerate(neighbors):
            start = (neighbor_slot + 1) * num_atom_types
            end = start + num_atom_types
            descriptors[:, site_index, start:end] = site_onehot[:, neighbor_site - 1, :]

    return descriptors


def _append_namod_sigma(descriptors, neighbor_sites):
    num_candidates, num_sites, _ = descriptors.shape
    sigma = np.zeros((num_candidates, num_sites, 1), dtype=np.float32)

    for site_index, neighbors in enumerate(neighbor_sites):
        if not neighbors:
            continue
        neighbor_indices = [neighbor_site - 1 for neighbor_site in neighbors]
        neighbor_descriptors = descriptors[:, neighbor_indices, :]
        neighbor_mean = np.mean(neighbor_descriptors, axis=1, keepdims=True)
        squared_norm = np.sum((neighbor_descriptors - neighbor_mean) ** 2, axis=2)
        sigma[:, site_index, 0] = np.sqrt(np.mean(squared_norm, axis=1))

    return np.concatenate([descriptors, sigma], axis=2)


def _validate_local_env_cache(cache_data, local_env_type):
    required_keys = {
        "X_raw",
        "X_env",
        "X_kernel",
        "X_site_labels",
        "atom_types",
        "candidate_ids",
        "neighbor_sites",
        "site_columns",
        "weights",
        "num_sites",
        "env_dim",
        "local_env_type",
        "schema_version",
    }
    if not isinstance(cache_data, dict):
        raise ValueError("Invalid local env cache schema: cache payload must be a dictionary")
    if "schema_version" not in cache_data:
        raise ValueError("Invalid local env cache schema: missing schema_version")
    if cache_data["schema_version"] != LOCAL_ENV_SCHEMA_VERSION:
        return False
    missing_keys = required_keys - set(cache_data)
    if missing_keys:
        missing = ", ".join(sorted(missing_keys))
        raise ValueError(f"Invalid local env cache schema: missing keys: {missing}")
    if cache_data["local_env_type"] != local_env_type:
        return False

    X_raw = np.asarray(cache_data["X_raw"])
    X_env = np.asarray(cache_data["X_env"])
    X_kernel = np.asarray(cache_data["X_kernel"])
    X_site_labels = np.asarray(cache_data["X_site_labels"])
    candidate_ids = np.asarray(cache_data["candidate_ids"])
    num_sites = int(cache_data["num_sites"])
    env_dim = int(cache_data["env_dim"])
    expected_kernel_dim = num_sites + num_sites * env_dim

    if X_raw.ndim != 2 or X_raw.shape[1] != num_sites:
        raise ValueError("Invalid local env cache schema: X_raw shape does not match num_sites")
    if X_env.ndim != 3 or X_env.shape[1:] != (num_sites, env_dim):
        raise ValueError("Invalid local env cache schema: X_env shape does not match num_sites/env_dim")
    if X_kernel.ndim != 2 or X_kernel.shape != (X_raw.shape[0], expected_kernel_dim):
        raise ValueError("Invalid local env cache schema: X_kernel shape is inconsistent")
    if X_site_labels.ndim != 2 or X_site_labels.shape != X_raw.shape:
        raise ValueError("Invalid local env cache schema: X_site_labels shape does not match X_raw")
    if candidate_ids.ndim != 1 or candidate_ids.shape[0] != X_raw.shape[0]:
        raise ValueError("Invalid local env cache schema: candidate_ids length does not match X_raw")
    if len(cache_data["site_columns"]) != num_sites:
        raise ValueError("Invalid local env cache schema: site_columns length does not match num_sites")
    if not isinstance(cache_data["weights"], dict):
        raise ValueError("Invalid local env cache schema: weights must be a dictionary")

    return True


def build_local_environment_candidate_features(
    local_env_type="NA",
    cache_filename="local_env_candidates.pkl",
    force_rebuild=False,
):
    """
    Build or load site-wise local environment descriptors for discrete BoTorch kernels.
    """
    local_env_type = _canonical_local_env_type(local_env_type)

    if os.path.exists(cache_filename) and not force_rebuild:
        with open(cache_filename, "rb") as f:
            cached_result = pickle.load(f)
        if _validate_local_env_cache(cached_result, local_env_type):
            print(f"(apx) Loaded local environment data from cache: {cache_filename}", flush=True)
            return cached_result
        print(
            f"(apx) Local environment cache is stale or for "
            f"{cached_result.get('local_env_type', 'unknown')}; rebuilding for {local_env_type}.",
            flush=True,
        )

    raw_codes, atom_types, candidate_ids, site_columns, site_values = _load_candidate_site_values()
    num_candidates, num_sites = raw_codes.shape
    neighbor_sites = _read_neighbor_sites_for_local_env(num_sites)
    weights = _read_local_env_weight_metadata()

    X_env = _build_na_site_descriptors(raw_codes, atom_types, neighbor_sites)
    if local_env_type == "NAmod":
        X_env = _append_namod_sigma(X_env, neighbor_sites)

    env_dim = int(X_env.shape[2])
    X_raw = np.ascontiguousarray(raw_codes, dtype=np.float32)
    X_kernel = np.concatenate(
        [X_raw, X_env.reshape(num_candidates, num_sites * env_dim)],
        axis=1,
    )
    X_kernel = np.ascontiguousarray(X_kernel, dtype=np.float32)

    result = {
        "X_raw": X_raw,
        "X_env": np.ascontiguousarray(X_env, dtype=np.float32),
        "X_kernel": X_kernel,
        "X_site_labels": np.ascontiguousarray(site_values),
        "atom_types": list(atom_types),
        "candidate_ids": np.ascontiguousarray(candidate_ids, dtype=np.int64),
        "neighbor_sites": neighbor_sites,
        "site_columns": site_columns,
        "weights": weights,
        "num_sites": num_sites,
        "env_dim": env_dim,
        "local_env_type": local_env_type,
        "schema_version": LOCAL_ENV_SCHEMA_VERSION,
    }

    with open(cache_filename, "wb") as f:
        pickle.dump(result, f)
    print(
        f"(apx) Saved local environment data to cache: {cache_filename} "
        f"(type={local_env_type}, sites={num_sites}, env_dim={env_dim})",
        flush=True,
    )

    return result


def load_or_build_local_environment_features(*args, **kwargs):
    return build_local_environment_candidate_features(*args, **kwargs)


def encode_options(encode="OH", weight=0, use_dimension_reduction=False, reduction_method=None, reduction_params=None):
    """
    Perform feature encoding for atomic configurations
    
    Parameters:
    - encode (str): Encoding method ('OH', 'NA', or 'NAmod').
      - 'OH': One-hot encoding.
      - 'NA': Neighbor atom encoding.
      - 'NAmod': Modified neighbor atom encoding.
    - weight (float, optional): Weight to use in (modified) neighbor atom encoding (default is 0).
    - use_dimension_reduction (bool): Whether to apply dimension reduction after encoding (default is False).
    - reduction_method (str): Method for dimension reduction ('PCA' or 'AUTOENCODER').
    - reduction_params (dict): Parameters for dimension reduction (default is None).

    Returns:
    - np.ndarray: Encoded input data.
    """

    # Try to load from cache
    cache_filename = "encoded_candidates.pkl"
    
    # Use cache if available
    if os.path.exists(cache_filename):
        try:
            with open(cache_filename, 'rb') as f:
                cached_result = pickle.load(f)
            print(f"(apx) Loaded encoded data from cache: {cache_filename}")
            return cached_result
        except Exception as e:
            print(f"(apx) Cache loading failed: {e}, recalculating...")
    
    # Load candidates of atomic configuration
    # Read all data for encoding
    df = pd.read_csv("candidates.csv")
    if df.shape[1] < 2:
        raise ValueError("candidates.csv must contain structure_id and at least one site column")

    # Exclude first column (structure ID) and get only atomic configuration data.
    # Treat labels consistently as strings so per-site atom sets are stable even
    # when pandas infers mixed dtypes.
    site_df = df.iloc[:, 1:].astype(str)
    x_candidates = site_df.to_numpy()
    num_candidates = x_candidates.shape[0]
    num_sites = x_candidates.shape[1]

    # Detect atom types independently for each site. This avoids allocating
    # columns for elements that never appear at a given site.
    site_atom_types = [
        sorted(site_df.iloc[:, site_index].unique().tolist())
        for site_index in range(num_sites)
    ]
    atom_types = sorted({atom_type for atoms in site_atom_types for atom_type in atoms})
    num_atom_types = len(atom_types)
    atom_type_to_index = {atom_type: index for index, atom_type in enumerate(atom_types)}

    site_feature_offsets = np.zeros(num_sites + 1, dtype=int)
    for site_index, atoms in enumerate(site_atom_types):
        site_feature_offsets[site_index + 1] = site_feature_offsets[site_index] + len(atoms)
    num_site_features = int(site_feature_offsets[-1])

    site_atom_to_column = []
    site_atom_global_indices = []
    for site_index, atoms in enumerate(site_atom_types):
        offset = site_feature_offsets[site_index]
        site_atom_to_column.append({
            atom_type: offset + atom_index
            for atom_index, atom_type in enumerate(atoms)
        })
        site_atom_global_indices.append(
            np.array([atom_type_to_index[atom_type] for atom_type in atoms], dtype=int)
        )

    namod_feature_offsets = np.zeros(num_sites + 1, dtype=int)
    for site_index, atoms in enumerate(site_atom_types):
        namod_feature_offsets[site_index + 1] = (
            namod_feature_offsets[site_index] + len(atoms) + 1
        )
    num_namod_features = int(namod_feature_offsets[-1])

    full_oh_dim = num_sites * num_atom_types
    print(
        f"(apx) Detected site-specific atom types: OH dimension "
        f"{num_site_features} (full site x atom dimension would be {full_oh_dim}).",
        flush=True
    )
    
    # Define neighbor_sites only if needed
    if encode in ["NA", "NAmod"]:
        neighbor_sites_data = read_card("apx.in", "NEIGHBOR_SITES")
        neighbor_sites = [list(map(int, line.split())) for line in neighbor_sites_data]
        if len(neighbor_sites) != num_sites:
            raise ValueError(
                f"NEIGHBOR_SITES must contain {num_sites} lines, "
                f"but {len(neighbor_sites)} lines were found"
            )
        for site_index, neighbors in enumerate(neighbor_sites):
            for neighbor_site in neighbors:
                if neighbor_site < 1 or neighbor_site > num_sites:
                    raise ValueError(
                        f"Invalid neighbor site {neighbor_site} for site {site_index + 1}; "
                        f"valid range is 1 to {num_sites}"
                    )

    def onehot():
        """Perform one-hot encoding."""
        print(f"(apx) Performing one-hot encoding...", flush=True)
        phi_candidates = np.zeros((num_candidates, num_site_features), dtype=float)
        for site_index in range(num_sites):
            site_values = x_candidates[:, site_index]
            for atom_type, column_index in site_atom_to_column[site_index].items():
                phi_candidates[:, column_index] = (site_values == atom_type).astype(int)
        print(f"(apx) One-hot encoding completed.", flush=True)
        return phi_candidates

    def neighbor_atom():
        """Perform neighbor atom encoding."""
        print(f"(apx) Performing neighbor atom encoding...", flush=True)
        phi_candidates = onehot()
        if weight != 0:
            for site_index in range(num_sites):
                atom_to_column = site_atom_to_column[site_index]
                for neighbor_site in neighbor_sites[site_index]:
                    neighbor_values = x_candidates[:, neighbor_site - 1]
                    for atom_type, column_index in atom_to_column.items():
                        phi_candidates[:, column_index] += (
                            weight * (neighbor_values == atom_type).astype(float)
                        )
        print(f"(apx) Neighbor atom encoding completed.", flush=True)
        return phi_candidates

    def modified_neighbor_atom():
        """Perform modified neighbor atom encoding."""
        print(f"(apx) Performing modified neighbor atom encoding...", flush=True)
        phi_candidates = neighbor_atom()
        phi_candidates_mod = np.zeros((num_candidates, num_namod_features), dtype=float)

        for site_index in range(num_sites):
            src_start = site_feature_offsets[site_index]
            src_end = site_feature_offsets[site_index + 1]
            dst_start = namod_feature_offsets[site_index]
            dst_end = dst_start + (src_end - src_start)
            phi_candidates_mod[:, dst_start:dst_end] = phi_candidates[:, src_start:src_end]

        def compact_site_block_as_global(site_index):
            """Return one compact site block projected onto the global atom basis."""
            block = np.zeros((num_candidates, num_atom_types), dtype=float)
            src_start = site_feature_offsets[site_index]
            src_end = site_feature_offsets[site_index + 1]
            block[:, site_atom_global_indices[site_index]] = phi_candidates[:, src_start:src_end]
            return block

        for site_index in range(num_sites):
            neighbors = neighbor_sites[site_index]
            if not neighbors:
                continue

            neighbor_blocks = [
                compact_site_block_as_global(neighbor_site - 1)
                for neighbor_site in neighbors
            ]
            neighbor_stack = np.stack(neighbor_blocks, axis=0)
            neighbor_mean = np.mean(neighbor_stack, axis=0)
            sigma = np.sqrt(
                np.mean(np.sum((neighbor_stack - neighbor_mean) ** 2, axis=2), axis=0)
            )
            phi_candidates_mod[:, namod_feature_offsets[site_index + 1] - 1] = sigma
        print(f"(apx) Modified neighbor atom encoding completed.", flush=True)
        return phi_candidates_mod

    def apply_dimension_reduction(features, method, params):
        """Apply dimension reduction to features."""
        input_dim = features.shape[1]
        
        if method.upper() == "PCA":
            if not SKLEARN_AVAILABLE:
                raise ImportError("scikit-learn is required for PCA but not available")
            
            from sklearn.decomposition import PCA
            
            # Use PCA for dimension reduction
            n_components = params.get('n_components', None)
            random_state = params.get('random_state', None)
            
            pca = PCA(n_components=n_components, random_state=random_state)
            
            # Fit and transform
            print(f"(apx) Training PCA for dimension reduction...", flush=True)
            latent_features = pca.fit_transform(features)
            
            return latent_features, pca, None
            
        elif method.upper() == "AUTOENCODER":
            if not TORCH_AVAILABLE:
                raise ImportError("PyTorch is required for Auto Encoder but not available")
            
            # Prepare data for PyTorch
            features_tensor = torch.FloatTensor(features)
            
            # Create auto encoder
            latent_dim = params['latent_dim']
            hidden_layers = params['hidden_layers']
            learning_rate = params['learning_rate']
            epochs = params['epochs']
            batch_size = params.get('batch_size', 1024)
            
            model = AutoEncoder(input_dim, latent_dim, hidden_layers)
            criterion = nn.MSELoss()
            optimizer = optim.Adam(model.parameters(), lr=learning_rate)
            
            # Create DataLoader for batch processing
            dataset = TensorDataset(features_tensor, features_tensor)
            dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
            
            # Training with batch processing
            print(f"(apx) Training Auto Encoder for dimension reduction...", flush=True)
            model.train()
            
            for epoch in range(epochs):
                epoch_loss = 0.0
                num_batches = 0
                
                for batch_features, batch_targets in dataloader:
                    optimizer.zero_grad()
                    decoded, latent = model(batch_features)
                    loss = criterion(decoded, batch_targets)
                    loss.backward()
                    optimizer.step()
                    
                    epoch_loss += loss.item()
                    num_batches += 1
                
                avg_loss = epoch_loss / num_batches
                
                if (epoch + 1) % 20 == 0:
                    print(f"(apx) Epoch {epoch + 1}/{epochs}, Average Loss: {avg_loss:.6f}", flush=True)
            
            # Extract latent features
            model.eval()
            with torch.no_grad():
                # Process in batches to handle large datasets
                latent_features_list = []
                for batch_features, _ in dataloader:
                    batch_latent = model.encode(batch_features).numpy()
                    latent_features_list.append(batch_latent)
                
                latent_features = np.vstack(latent_features_list)
            
            return latent_features, model, None
        else:
            raise ValueError(f"Unknown dimension reduction method: {method}")

    # Perform encoding
    if encode == "OH":
        result = onehot()
    elif encode == "NA":
        result = neighbor_atom()
    elif encode == "NAmod":
        result = modified_neighbor_atom()
    else:
        raise ValueError("Invalid option. Choose 'OH', 'NA', or 'NAmod'.")
    
    # Apply dimension reduction if requested
    if use_dimension_reduction and reduction_method and reduction_params:
        try:
            # Store original dimension before reduction
            original_dim = result.shape[1]
            result, reduction_model, encoder_model = apply_dimension_reduction(result, reduction_method, reduction_params)
            
            # Print dimension reduction results
            output_dim = result.shape[1]
            print(f"(apx) Dimension reduction completed.", flush=True)
            print(f"(apx) Input dimension: {original_dim} -> Latent dimension: {output_dim}", flush=True)
            
            if reduction_method.upper() == "PCA" and reduction_model is not None:
                explained_variance_ratio = reduction_model.explained_variance_ratio_.sum()
                print(f"(apx) Explained variance ratio: {explained_variance_ratio:.6f}", flush=True)
                
        except ImportError as e:
            print(f"(apx) Error: {e}", flush=True)
            print(f"(apx) Program terminated due to missing required library.", flush=True)
            raise SystemExit(1)
        except Exception as e:
            print(f"(apx) Error during dimension reduction: {e}", flush=True)
            print(f"(apx) Program terminated due to dimension reduction error.", flush=True)
            raise SystemExit(1)
    
    # Save to cache
    try:
        print(f"(apx) Saving encoded data to cache: {cache_filename}", flush=True)
        
        with open(cache_filename, 'wb') as f:
            pickle.dump(result, f)
        
        print(f"(apx) Saved encoded data to cache: {cache_filename}", flush=True)
    except Exception as e:
        print(f"(apx) Cache saving failed: {e}", flush=True)
    
    return result
