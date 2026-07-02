</div>
<div align="center">
  <img src="logo/logo_banner.png" alt="PyAPX Logo" width="450">
</div>

[![arXiv](https://img.shields.io/badge/arXiv-2511.17972-b31b1b.svg)](http://arxiv.org/abs/2511.17972)

# PyAPX - Python Toolkit for Atomic Configuration Pattern Exploration

Unofficial experimental fork of PyAPX. For stable or citable use, please use the official PyAPX repository.

## Requirements

- Python 3.9.6 (tested)
- numpy == 1.26.4
- pandas == 2.2.3
- [physbo](https://github.com/issp-center-dev/PHYSBO) == 3.1.0
- [BoTorch](https://botorch.org/) == 0.10.0
- [GPyTorch](https://gpytorch.ai/) == 1.11
- [PyTorch](https://pytorch.org/) >= 1.13.1
- scikit-learn == 1.6.1
- scipy == 1.13.1

### DFT code

- [Quantum ESPRESSO](https://www.quantum-espresso.org/) v7.3.1 (tested)

## Installation

```bash
git clone https://github.com/a-ksb/PyAPX.git
cd PyAPX
pip install -r requirements.txt
```

## Quick Start

### Example Using an Ising Model (Instead of a DFT Code)

For details about this material system, please refer to articles <sup>[2,3]</sup>.

```bash
cd PyAPX/examples/H_GaN0001_6x6
python gen_candidates.py
PYTHONPATH=../.. python -m pyapx.cli
python visualize_results.py
```

## How to Use

### Example Using Quantum ESPRESSO

For details about this material system, please refer to the article <sup>[4]</sup>.

```bash
cd PyAPX/examples/h-BCN_3x3
```

Users need to prepare the following three files:

- `apx.in`: PyAPX input file
- `qe_template.in`: Quantum ESPRESSO input template
- `candidates.csv`: List of candidate atomic configurations

In `apx.in`, users can set parameters as follows:

```ini
ENCODING = True
RANDOM_SAMPLING = 2    # the number of random sampling
BAYES_SAMPLING = 2    # the number of Bayesian sampling
ENERGY_EVALUATOR = qe    # "qe" or custom:module.function
#PARALLEL_COMMAND = mpirun -np 8    # mpi setting for dft, Default: "mpiexec"
OPTIMIZER = physbo
```

When ENCODING = True, `encoded_candidates.pkl` is generated, which is required for Bayesian sampling.
Then, schedule the number of random sampling and Bayesian sampling iterations. At least two initial data points are required for Bayesian sampling.
PyAPX (v1.0.0) currently supports only Quantum ESPRESSO as a DFT code.
As an energy evaluator, you can also specify user-defined functions instead of DFT codes (see Quick Start example).

```ini
# settings for encoding
ENCODE = NAmod    # encoding method: "OH", "NA" or "NAmod"
WEIGHT = 0.3    # parameter for "NA" and "NAmod"
```

You can choose from the following encoding methods: one-hot (OH) encoding, neighbor-atom (NA) encoding, and modified neighbor-atom (NAmod) encoding. For details, please refer to the article <sup>[1]</sup>.

```ini
# settings for physbo
SCORE = TS    # acquisition function: "TS", "EI" or "PI"
NUM_RAND_BASIS = 3000    # the number of basis functions
```

You can choose from the following acquisition functions: Thompson Sampling (TS), Expected Improvement (EI), and Probability of Improvement (PI). For details, please refer to the [PHYSBO documentation](https://issp-center-dev.github.io/PHYSBO/manual/master/en/index.html).

```ini
# settings for botorch
OPTIMIZER = botorch
SCORE = EI    # acquisition function: "TS" or "EI"
BOTORCH_DEVICE = auto    # "auto", "cuda", "cuda:0" or "cpu"
BOTORCH_DTYPE = float32    # "float32" or "float64"
BOTORCH_BATCH_SIZE = 65536    # candidates scored per GPU/CPU batch
BOTORCH_MAXITER = 500    # GP hyperparameter fitting iterations
BOTORCH_INPUT_TRANSFORM = unit_cube    # "unit_cube", "standardize", or "none"
BOTORCH_GP_KERNEL = default    # "default", "hamming", "ot", or "hamming_ot"
BOTORCH_STANDARDIZE_Y = False    # False keeps raw objective values
BOTORCH_DISCRETE_CHUNK_SIZE = 4096    # row chunk size for Hamming / mixed discrete kernels
BOTORCH_HAMMING_LENGTHSCALE = 0.2    # initial Hamming-kernel lengthscale
BOTORCH_OT_LENGTHSCALE = 0.5    # initial optimal-transport-kernel lengthscale
BOTORCH_OT_ATOM_MISMATCH_PENALTY = 1.0    # extra OT cost for matching different atom types
BOTORCH_OT_SINKHORN_EPSILON = 0.05    # entropic regularization for Sinkhorn OT
BOTORCH_OT_SINKHORN_ITERATIONS = 30    # Sinkhorn update count
BOTORCH_OT_CHUNK_SIZE = 1024    # row chunk size for OT distance evaluation
BOTORCH_USE_LOCAL_ENV = False    # add NA / NAmod site-wise descriptors to discrete kernels
BOTORCH_LOCAL_ENV_TYPE = NA    # "NA" or "NAmod"
BOTORCH_ENV_DISTANCE = l1    # "l1" or "l2" descriptor distance
BOTORCH_ENV_LENGTHSCALE = 0.2    # Hamming local-env kernel lengthscale
BOTORCH_OT_ENV_MISMATCH_PENALTY = 1.0    # OT cost weight for descriptor mismatch
BOTORCH_LOCAL_ENV_CACHE = local_env_candidates.pkl
BOTORCH_SA_SCREENING = False    # use simulated annealing to pre-screen candidates
BOTORCH_SA_SCORE = EI    # "EI" or "TS" for screening
BOTORCH_SA_INITIAL_POOL_SIZE = 65536    # random candidates scored before SA
BOTORCH_SA_CHAINS = 1024    # parallel SA chains
BOTORCH_SA_STEPS = 500    # annealing steps per chain
BOTORCH_SA_EVAL_BATCH_SIZE = 4096    # scoring batch size during SA
BOTORCH_SA_INITIAL_TEMPERATURE = 1.0
BOTORCH_SA_FINAL_TEMPERATURE = 0.01
BOTORCH_SA_RANDOM_FRACTION = 0.05    # probability of a global random proposal
BOTORCH_SA_SWAP_NEIGHBORS = True    # use composition-preserving site swaps
#BOTORCH_XI = 0.0    # improvement margin for EI
```

When `OPTIMIZER = botorch`, PyAPX fits a BoTorch/GPyTorch Gaussian process and scores the discrete candidate list in batches on the selected PyTorch device. Use a CUDA-enabled PyTorch installation if `BOTORCH_DEVICE = auto` or `cuda` should run on an NVIDIA GPU. `BOTORCH_GP_KERNEL = default` uses BoTorch's `SingleTaskGP` with encoded/PCA features. `BOTORCH_GP_KERNEL = hamming`, `ot`, and `hamming_ot` consume the raw atomic labels from `candidates.csv` instead of encoded/PCA features by default. Set `BOTORCH_USE_LOCAL_ENV = True` to augment those discrete kernels with site-wise `NA` or `NAmod` local environment descriptors; this requires `NEIGHBOR_SITES` and creates `local_env_candidates.pkl` automatically. The cache also stores inspection metadata such as `candidate_ids`, `site_columns`, `weights`, and `X_site_labels`; `WEIGHT` is recorded there but is not included in the kernel input. The optimal-transport kernel uses normalized shortest-path distances from `NEIGHBOR_SITES` and a Sinkhorn approximation; `hamming_ot` adds both discrete kernels. `BOTORCH_SA_SCREENING = True` evaluates a large random pool and then runs multi-chain simulated annealing over uncalculated candidates, using the selected acquisition score as a cheap surrogate before returning one structure ID for expensive evaluation. The original `OPTIMIZER = physbo` path is unchanged.

```ini
# settings for "NA", "NAmod"
NEIGHBOR_SITES
10   12   18    # site_1's neighbors are site_10, site_12 and site_18
10   11   16    # site_2's neighbors are site_10, site_11 and site_16
11   12   17    # site_3's neighbors are site_11, site_12 and site_17
      .
      .
      .
```

When using NA or NAmod encoding, or `BOTORCH_USE_LOCAL_ENV = True`, list the neighboring sites in order from site_1 to define the site network.

In `qe_template.in`, users need to include atomic coordinates in the [Quantum ESPRESSO input file](https://www.quantum-espresso.org/Doc/INPUT_PW.html) as follows:

```ini
&control
calculation = 'vc-relax'    # DFT calc settings
      .
      .
      .

ATOMIC_SPECIES
   C   12.01060   C_paw.UPF
   B   10.81350   B_paw.UPF
   N   14.00650   N_paw.UPF
ATOMIC_POSITIONS {crystal}
X   0.0000   0.0000   0.0000    # site_1
X   0.0000   0.3333   0.0000    # site_2
X   0.0000   0.6667   0.0000    # site_3
      .
      .
      .
```

Here, sites that are subject to atomic rearrangement should be set with the element as a wildcard "X".

In this example, `candidates.csv` is generated using the following enumeration code, which contains as many as 5,717,712 candidate configurations.

```bash
julia gen_candidates.jl    # julia 1.10.3 (tested)
```

```ini
structure_id,site_1,site_2,site_3,...,site_18
0,C,B,B,B,B,B,B,C,C,C,C,C,N,N,N,N,N,N
1,C,B,B,B,B,B,B,C,C,C,C,N,C,N,N,N,N,N
2,C,B,B,B,B,B,B,C,C,C,C,N,N,C,N,N,N,N
      .
      .
      .
5717711,C,N,N,N,N,N,N,C,C,C,C,C,B,B,B,B,B,B
```

Users are required to prepare `candidates.csv` in the format shown above.
The sites in `candidates.csv` correspond to the sites specified by "X" in `qe_template.in` in order.

Finally, ensure that the three files `apx.in`, `qe_template.in`, and `candidates.csv` are present in the current directory, and run the following command with the PYTHONPATH set to the project root directory where the `pyapx` directory is located.

```bash
PYTHONPATH=../.. python -m pyapx.cli
```

Output files:

- `samples.csv`: Sampling results and energy values
- `dft_calc/qe_sample_[ID].in`: Generated QE input files
- `dft_calc/qe_sample_[ID].out`: QE calculation output files

## Fork notice

This repository is an unofficial experimental fork of PyAPX. For stable, reproducible, or citable use, please use the official PyAPX repository.

Modifications in this fork:
Copyright (c) 2026 Tatoshi YONEMORI
