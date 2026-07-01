"""
BoTorch kernel helpers for PyAPX.
"""

import warnings
import numpy as np


def is_discrete_botorch_kernel(gp_kernel):
    """
    Return True when the BoTorch kernel consumes raw atomic site labels.
    """
    gp_kernel = str(gp_kernel).strip().lower()
    return gp_kernel in ["hamming", "ot", "hamming_ot"]


def build_site_graph_distance(num_sites):
    """
    Build normalized shortest-path distances between sites from NEIGHBOR_SITES.
    """
    try:
        from .utils import read_card
        neighbor_sites_data = read_card("apx.in", "NEIGHBOR_SITES")
    except Exception:
        neighbor_sites_data = []

    distances = np.full((num_sites, num_sites), np.inf, dtype=np.float64)
    np.fill_diagonal(distances, 0.0)

    parsed_any_edge = False
    for site_index, raw_line in enumerate(neighbor_sites_data[:num_sites]):
        line = raw_line.split("#", 1)[0].strip()
        try:
            neighbors = [int(value) for value in line.split()]
        except ValueError:
            continue

        for neighbor in neighbors:
            neighbor_index = neighbor - 1
            if 0 <= neighbor_index < num_sites and neighbor_index != site_index:
                distances[site_index, neighbor_index] = 1.0
                distances[neighbor_index, site_index] = 1.0
                parsed_any_edge = True

    if parsed_any_edge:
        for k in range(num_sites):
            distances = np.minimum(distances, distances[:, [k]] + distances[[k], :])

    if not parsed_any_edge or not np.isfinite(distances).all():
        distances = np.ones((num_sites, num_sites), dtype=np.float64)
        np.fill_diagonal(distances, 0.0)

    max_distance = float(np.max(distances))
    if max_distance > 0.0:
        distances /= max_distance

    return distances


def _initialize_gp_output_hyperparameters(torch, model, train_Y):
    """
    Initialize mean, outputscale, and noise from objective statistics.
    """
    flat_y = train_Y.reshape(-1)
    y_std = flat_y.std(unbiased=False).clamp_min(1.0e-8)

    with torch.no_grad():
        if hasattr(model.mean_module, "constant"):
            model.mean_module.constant = flat_y.median()

        covar_module = model.covar_module
        if hasattr(covar_module, "outputscale"):
            covar_module.outputscale = y_std.pow(2).clamp_min(1.0e-8)

        if hasattr(model, "likelihood") and hasattr(model.likelihood, "noise"):
            noise = y_std.div(10.0).pow(2)
            noise_constraint = getattr(model.likelihood.noise_covar, "raw_noise_constraint", None)
            if noise_constraint is not None:
                lower_bound = getattr(noise_constraint, "lower_bound", None)
                if lower_bound is not None:
                    lower_bound = lower_bound.to(device=train_Y.device, dtype=train_Y.dtype).reshape(())
                    noise = noise.clamp_min(float(lower_bound.item()))
            model.likelihood.noise = noise


def _split_raw_and_env(X, num_sites, env_dim):
    """
    Split a local-env kernel input into raw site labels and site-wise descriptors.
    """
    num_sites = int(num_sites)
    env_dim = int(env_dim)
    expected_dim = num_sites + num_sites * env_dim
    if X.shape[-1] < expected_dim:
        raise ValueError(
            f"Local env kernel input has dimension {X.shape[-1]}, "
            f"but expected at least {expected_dim}"
        )

    raw = X[..., :num_sites]
    env_flat = X[..., num_sites:expected_dim]
    env = env_flat.reshape(X.shape[:-1] + (num_sites, env_dim))
    return raw, env


def _pairwise_env_distance(torch, env1, env2, distance_type):
    """
    Compute descriptor distances with broadcasting over all leading dimensions.
    """
    distance_type = str(distance_type).strip().lower()
    difference = env1 - env2

    if distance_type in ["l1", "manhattan"]:
        return difference.abs().sum(dim=-1)
    if distance_type in ["l2", "euclidean"]:
        return torch.linalg.vector_norm(difference, ord=2, dim=-1)

    raise ValueError(
        f"Unsupported BOTORCH_ENV_DISTANCE: {distance_type}. "
        "Choose l1 or l2."
    )


def _make_hamming_kernel(
    torch,
    gpytorch,
    chunk_size,
    lengthscale,
    use_local_env=False,
    num_sites=None,
    env_dim=None,
    env_distance="l1",
):
    """
    Create an exponential Hamming-distance kernel for site-label vectors.
    """
    class HammingDistanceKernel(gpytorch.kernels.Kernel):
        has_lengthscale = True

        def __init__(
            self,
            chunk_size,
            use_local_env,
            num_sites,
            env_dim,
            env_distance,
            **kwargs,
        ):
            super().__init__(**kwargs)
            self.chunk_size = max(1, int(chunk_size))
            self.use_local_env = bool(use_local_env)
            self.num_sites = None if num_sites is None else int(num_sites)
            self.env_dim = None if env_dim is None else int(env_dim)
            self.env_distance = env_distance

            if self.use_local_env and (self.num_sites is None or self.env_dim is None):
                raise ValueError("Local env Hamming kernel requires num_sites and env_dim")

        def _split_inputs(self, x1, x2):
            if not self.use_local_env:
                return x1.round(), x2.round(), None, None

            raw1, env1 = _split_raw_and_env(x1, self.num_sites, self.env_dim)
            raw2, env2 = _split_raw_and_env(x2, self.num_sites, self.env_dim)
            return raw1.round(), raw2.round(), env1, env2

        def _aligned_env_distance(self, env1, env2):
            distance = torch.zeros(env1.shape[0], device=env1.device, dtype=env1.dtype)
            for site_index in range(self.num_sites):
                distance = distance + _pairwise_env_distance(
                    torch,
                    env1[:, site_index, :],
                    env2[:, site_index, :],
                    self.env_distance,
                )
            return distance / float(self.num_sites)

        def _pairwise_aligned_env_distance(self, env1, env2):
            distance = torch.zeros(
                env1.shape[0],
                env2.shape[0],
                device=env1.device,
                dtype=env1.dtype,
            )
            for site_index in range(self.num_sites):
                distance = distance + _pairwise_env_distance(
                    torch,
                    env1[:, None, site_index, :],
                    env2[None, :, site_index, :],
                    self.env_distance,
                )
            return distance / float(self.num_sites)

        def forward(self, x1, x2, diag=False, **params):
            raw1, raw2, env1, env2 = self._split_inputs(x1, x2)
            lengthscale = self.lengthscale.to(device=x1.device, dtype=x1.dtype).reshape(()).clamp_min(1.0e-8)

            if diag:
                distance = (raw1 != raw2).to(dtype=x1.dtype).mean(dim=-1)
                if self.use_local_env:
                    distance = distance + self._aligned_env_distance(env1, env2)
                return torch.exp(-distance / lengthscale)

            n1 = raw1.shape[-2]
            n2 = raw2.shape[-2]
            result = torch.empty(n1, n2, device=x1.device, dtype=x1.dtype)
            for start in range(0, n1, self.chunk_size):
                end = min(start + self.chunk_size, n1)
                distance = (
                    raw1[start:end].unsqueeze(1) != raw2.unsqueeze(0)
                ).to(dtype=x1.dtype).mean(dim=-1)
                if self.use_local_env:
                    distance = distance + self._pairwise_aligned_env_distance(
                        env1[start:end],
                        env2,
                    )
                result[start:end] = torch.exp(-distance / lengthscale)

            return result

    kernel = HammingDistanceKernel(
        chunk_size=chunk_size,
        use_local_env=use_local_env,
        num_sites=num_sites,
        env_dim=env_dim,
        env_distance=env_distance,
    )
    kernel.lengthscale = torch.as_tensor(lengthscale, dtype=torch.float64).clamp_min(1.0e-8)
    return kernel


def _make_ot_kernel(
    torch,
    gpytorch,
    site_cost,
    atom_mismatch_penalty,
    sinkhorn_epsilon,
    sinkhorn_iterations,
    chunk_size,
    lengthscale,
    use_local_env=False,
    num_sites=None,
    env_dim=None,
    env_distance="l1",
    ot_env_mismatch_penalty=1.0,
):
    """
    Create a Sinkhorn-OT kernel over raw site-label configurations.
    """
    class OptimalTransportKernel(gpytorch.kernels.Kernel):
        has_lengthscale = True

        def __init__(
            self,
            site_cost_tensor,
            atom_mismatch_penalty,
            sinkhorn_epsilon,
            sinkhorn_iterations,
            chunk_size,
            use_local_env,
            num_sites,
            env_dim,
            env_distance,
            ot_env_mismatch_penalty,
            **kwargs,
        ):
            super().__init__(**kwargs)
            self.register_buffer("site_cost", site_cost_tensor)
            self.atom_mismatch_penalty = float(atom_mismatch_penalty)
            self.sinkhorn_epsilon = max(float(sinkhorn_epsilon), 1.0e-6)
            self.sinkhorn_iterations = max(1, int(sinkhorn_iterations))
            self.chunk_size = max(1, int(chunk_size))
            self.use_local_env = bool(use_local_env)
            self.num_sites = None if num_sites is None else int(num_sites)
            self.env_dim = None if env_dim is None else int(env_dim)
            self.env_distance = env_distance
            self.ot_env_mismatch_penalty = float(ot_env_mismatch_penalty)

            if self.use_local_env and (self.num_sites is None or self.env_dim is None):
                raise ValueError("Local env OT kernel requires num_sites and env_dim")

        def _sinkhorn_distance(self, cost):
            num_sites = cost.shape[-1]
            epsilon = torch.as_tensor(self.sinkhorn_epsilon, device=cost.device, dtype=cost.dtype)
            log_a = -torch.log(torch.as_tensor(float(num_sites), device=cost.device, dtype=cost.dtype))
            log_b = log_a
            log_kernel = -cost / epsilon
            u = torch.zeros(cost.shape[:-1], device=cost.device, dtype=cost.dtype)
            v = torch.zeros(cost.shape[:-2] + (num_sites,), device=cost.device, dtype=cost.dtype)

            for _ in range(self.sinkhorn_iterations):
                u = log_a - torch.logsumexp(log_kernel + v.unsqueeze(-2), dim=-1)
                v = log_b - torch.logsumexp(log_kernel.transpose(-2, -1) + u.unsqueeze(-2), dim=-1)

            log_plan = log_kernel + u.unsqueeze(-1) + v.unsqueeze(-2)
            plan = torch.exp(log_plan)
            return (plan * cost).sum(dim=(-2, -1))

        def _split_inputs(self, x1, x2):
            if not self.use_local_env:
                return x1.round(), x2.round(), None, None

            raw1, env1 = _split_raw_and_env(x1, self.num_sites, self.env_dim)
            raw2, env2 = _split_raw_and_env(x2, self.num_sites, self.env_dim)
            return raw1.round(), raw2.round(), env1, env2

        def _distance_to_one(self, x1_chunk, x2_one, env1_chunk=None, env2_one=None):
            site_cost = self.site_cost.to(device=x1_chunk.device, dtype=x1_chunk.dtype)
            mismatch = (
                x1_chunk.unsqueeze(-1) != x2_one.reshape(1, 1, -1)
            ).to(dtype=x1_chunk.dtype)
            cost = site_cost.unsqueeze(0) + self.atom_mismatch_penalty * mismatch
            if self.use_local_env:
                env_cost = _pairwise_env_distance(
                    torch,
                    env1_chunk[:, :, None, :],
                    env2_one.reshape(1, 1, self.num_sites, self.env_dim),
                    self.env_distance,
                )
                cost = cost + self.ot_env_mismatch_penalty * env_cost
            return self._sinkhorn_distance(cost)

        def forward(self, x1, x2, diag=False, **params):
            raw1, raw2, env1, env2 = self._split_inputs(x1, x2)
            lengthscale = self.lengthscale.to(device=x1.device, dtype=x1.dtype).reshape(()).clamp_min(1.0e-8)

            if diag:
                result = torch.empty(raw1.shape[-2], device=x1.device, dtype=x1.dtype)
                identical_pairs = None
                if self.use_local_env:
                    raw_equal = (raw1 == raw2).all(dim=-1)
                    env_equal = torch.isclose(
                        env1,
                        env2,
                        rtol=1.0e-12,
                        atol=1.0e-12,
                    ).reshape(env1.shape[0], -1).all(dim=-1)
                    identical_pairs = raw_equal & env_equal
                for start in range(0, raw1.shape[-2], self.chunk_size):
                    end = min(start + self.chunk_size, raw1.shape[-2])
                    distances = []
                    for offset in range(start, end):
                        if self.use_local_env:
                            distances.append(
                                self._distance_to_one(
                                    raw1[offset:offset + 1],
                                    raw2[offset],
                                    env1[offset:offset + 1],
                                    env2[offset],
                                ).reshape(())
                            )
                        else:
                            distances.append(
                                self._distance_to_one(
                                    raw1[offset:offset + 1],
                                    raw2[offset],
                                ).reshape(())
                            )
                    result[start:end] = torch.stack(distances)
                if identical_pairs is not None:
                    result = torch.where(identical_pairs, torch.zeros_like(result), result)
                return torch.exp(-result / lengthscale)

            n1 = raw1.shape[-2]
            n2 = raw2.shape[-2]
            result = torch.empty(n1, n2, device=x1.device, dtype=x1.dtype)
            for column in range(n2):
                for start in range(0, n1, self.chunk_size):
                    end = min(start + self.chunk_size, n1)
                    if self.use_local_env:
                        distance = self._distance_to_one(
                            raw1[start:end],
                            raw2[column],
                            env1[start:end],
                            env2[column],
                        )
                    else:
                        distance = self._distance_to_one(raw1[start:end], raw2[column])
                    result[start:end, column] = torch.exp(-distance / lengthscale)

            return result

    site_cost_tensor = torch.as_tensor(site_cost, dtype=torch.float64)
    kernel = OptimalTransportKernel(
        site_cost_tensor=site_cost_tensor,
        atom_mismatch_penalty=atom_mismatch_penalty,
        sinkhorn_epsilon=sinkhorn_epsilon,
        sinkhorn_iterations=sinkhorn_iterations,
        chunk_size=chunk_size,
        use_local_env=use_local_env,
        num_sites=num_sites,
        env_dim=env_dim,
        env_distance=env_distance,
        ot_env_mismatch_penalty=ot_env_mismatch_penalty,
    )
    kernel.lengthscale = torch.as_tensor(lengthscale, dtype=torch.float64).clamp_min(1.0e-8)
    return kernel


def _create_discrete_botorch_kernel(torch, gpytorch, gp_kernel, kernel_context):
    """
    Create a Hamming, OT, or mixed discrete kernel.
    """
    gp_kernel = str(gp_kernel).strip().lower()
    chunk_size = kernel_context["discrete_chunk_size"]
    use_local_env = bool(kernel_context.get("use_local_env", False))
    num_sites = kernel_context.get("num_sites")
    env_dim = kernel_context.get("env_dim")
    env_distance = kernel_context.get("env_distance", "l1")
    env_lengthscale = kernel_context.get("env_lengthscale", kernel_context["hamming_lengthscale"])
    hamming_kernel = None
    ot_kernel = None

    if gp_kernel in ["hamming", "hamming_ot"]:
        hamming_kernel = _make_hamming_kernel(
            torch,
            gpytorch,
            chunk_size=chunk_size,
            lengthscale=env_lengthscale if use_local_env else kernel_context["hamming_lengthscale"],
            use_local_env=use_local_env,
            num_sites=num_sites,
            env_dim=env_dim,
            env_distance=env_distance,
        )

    if gp_kernel in ["ot", "hamming_ot"]:
        ot_kernel = _make_ot_kernel(
            torch,
            gpytorch,
            site_cost=kernel_context["site_cost"],
            atom_mismatch_penalty=kernel_context["ot_atom_mismatch_penalty"],
            sinkhorn_epsilon=kernel_context["ot_sinkhorn_epsilon"],
            sinkhorn_iterations=kernel_context["ot_sinkhorn_iterations"],
            chunk_size=kernel_context["ot_chunk_size"],
            lengthscale=kernel_context["ot_lengthscale"],
            use_local_env=use_local_env,
            num_sites=num_sites,
            env_dim=env_dim,
            env_distance=env_distance,
            ot_env_mismatch_penalty=kernel_context.get("ot_env_mismatch_penalty", 1.0),
        )

    if hamming_kernel is not None and ot_kernel is not None:
        return hamming_kernel + ot_kernel
    if hamming_kernel is not None:
        return hamming_kernel
    if ot_kernel is not None:
        return ot_kernel

    raise ValueError(f"Unsupported discrete BOTORCH_GP_KERNEL: {gp_kernel}")


def create_botorch_model(torch, gpytorch, SingleTaskGP, train_X, train_Y, gp_kernel, kernel_context=None):
    """
    Create a BoTorch GP model.
    """
    gp_kernel = str(gp_kernel).strip().lower()

    if is_discrete_botorch_kernel(gp_kernel):
        if kernel_context is None:
            raise ValueError("Discrete BoTorch kernels require kernel_context")
        mean_module = gpytorch.means.ConstantMean()
        covar_module = gpytorch.kernels.ScaleKernel(
            _create_discrete_botorch_kernel(torch, gpytorch, gp_kernel, kernel_context)
        )
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Input data is not contained to the unit cube.*")
            warnings.filterwarnings("ignore", message="Input data is not standardized.*")
            model = SingleTaskGP(
                train_X,
                train_Y,
                mean_module=mean_module,
                covar_module=covar_module,
            )
        _initialize_gp_output_hyperparameters(torch, model, train_Y)
        return model

    if gp_kernel == "default":
        return SingleTaskGP(train_X, train_Y)

    raise ValueError(
        f"Unsupported BOTORCH_GP_KERNEL: {gp_kernel}. "
        "Choose default, hamming, ot, or hamming_ot."
    )
