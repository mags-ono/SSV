## Code for paper "Data-Driven Estimation of Structured Singular Values"

This repository contains MATLAB code associated with the paper:

**“Data-Driven Estimation of Structured Singular Values”**  
Margarita A. Guerrero, Braghadeesh Lakshminarayanan, Cristian R. Rojas  
📄 [arXiv link](https://arxiv.org/abs/2503.13410)  
✅ Accepted for publication in *IEEE Control Systems Letters (LCSS), 2025 Edition*.

### Overview

This code implements a **data-driven method** to compute a **lower bound on the structured singular value (SSV)** of a discrete-time linear time-invariant system. The method relies exclusively on **input-output time-domain experiments**, requiring **no model of the plant**, and assumes prior knowledge of the **structure of the uncertainty block**.

The approach is based on a frequency-wise **power iteration scheme**, where time-domain experiments are used to approximate frequency responses via the Discrete Fourier Transform (DFT), in both forward and backward directions. This procedure estimates two bounds, denoted as \( \tilde{\mu} \) and \( \bar{\mu} \), which converge under mild assumptions to a valid lower bound of the structured singular value.

As a baseline comparison, the method also computes the model-based lower bound via MATLAB's `mussv` function from the Robust Control Toolbox, using the exact same block structure.

---

### Contents

The repository provides a single, **self-contained MATLAB file**:

```
SSV.m
```

This file contains all necessary routines to simulate the proposed method. It includes:

- Utility functions for evaluating transfer functions and simulating the plant
- The main algorithm that computes \( \tilde{\mu} \) and \( \bar{\mu} \) using data only
- An example corresponding to **Figure 3** of the paper

---

### Main Functions

Below is a list of the key functions defined in `SSV.m`, along with their descriptions:

#### `power_method_convergence(G0, B_init, W_init, num_iters, r, mm, tol, tol2, N, sigma)`

Performs the iterative data-driven power method to estimate the structured singular value.  
- **Inputs**:  
  - `G0`: Plant as a cell array of discrete transfer functions  
  - `B_init`, `W_init`: Initial frequency-domain vectors  
  - `num_iters`: Maximum number of iterations  
  - `r`, `mm`: Block structure definitions  
    - `r = [s, d₁, d₂, ..., dₛ]` where `s` is the number of scalar complex blocks and `dᵢ` their dimensions  
    - `mm = [f, d₁, d₂, ..., d_f]` where `f` is the number of full complex blocks and `dᵢ` their dimensions  
    - Example structures used in the paper:

```
cases = { 
 [1, 3], [0, 0];    % Case 1: Scalar block 3x3 
 [0, 0], [1, 3];    % Case 2: Full block 3x3 
 [1, 2], [1, 1];    % Case 3: Scalar 2x2, Full 1x1 
 [1, 1], [1, 2];    % Case 4: Scalar 1x1, Full 2x2 (with noise) 
 [2, 1, 1], [1, 1]; % Case 5: Scalars 1x1, 1x1 and Full 2x2 
 [1, 1], [2, 1, 1]; % Case 6: Scalar 1x1 and Fulls 1x1, 1x1 
 [3, 1, 1, 1], [0, 0];  % Case 7: 3 separate scalar blocks 
 [0, 0], [3, 1, 1, 1]   % Case 8: 3 separate full blocks 
};
```

- **Other Inputs**:  
  - `tol`, `tol2`: Convergence tolerances  
  - `N`: Number of frequency points  
  - `sigma`: Noise standard deviation for simulation  
- **Outputs**:  
  - `mu_tilde`, `mu_bar`: Sequences of bounds over iterations  
  - `B_final`, `W_final`: Final frequency-domain vectors

> For how these block structures are used in `mussv`, refer to the MATLAB documentation:  
> https://se.mathworks.com/help/robust/ref/mussv.html

#### `generate_random_G0(n_outputs, n_inputs, n_states, seed)`

Generates a stable discrete-time MIMO system in transfer function format.  
- **Inputs**: Output and input dimensions, number of states, random seed  
- **Output**: `G0` in cell structure with fields `b` (numerator) and `a` (denominator)

#### `simulate_plant_general_discrete(G0, u_time, sigma)`

Simulates the time-domain response of a discrete LTI system with additive Gaussian noise.  
- **Inputs**:  
  - `G0`: Plant in cell format  
  - `u_time`: Input time-domain signal  
  - `sigma`: Noise standard deviation  
- **Output**: Output signal over time

#### `evaluate_discrete_tf(G0, omega)`

Evaluates the frequency response of `G0` at a given frequency \( \omega \).  
- **Inputs**:  
  - `G0`: Cell-based plant  
  - `omega`: Frequency value (in radians/sample)  
- **Output**: Frequency response matrix \( G(e^{j\omega}) \)

---

### Example Included: Frequency Response (Figure 3)

The file `SSV.m` includes an example that reproduces **Figure 3** from the paper.

This example simulates a system with:
- 3 outputs and inputs (`n = 3`)
- Uncertainty structure:  
  - One scalar block of size 1×1 → `r = [1, 1]`  
  - One full block of size 2×2 → `m = [1, 2]`

The following steps are executed:
- A random system `G0` is generated using `generate_random_G0`.
- Input/output data is created from experiments using `simulate_plant_general_discrete`.
- The power method is run using `power_method_convergence`.
- The `mussv` lower bound is computed using:

```matlab
[bounds, ~] = mussv(G_mussv, mussv_block);
mussv_lower_values = bounds(1).ResponseData(:);
[max_mussv_value, max_mussv_idx] = max(mussv_lower_values);
max_mussv_freq = 2*pi - freqs(max_mussv_idx);
```

This example is freely executable and provides a minimal working demonstration of the proposed method. It allows direct visual comparison between the data-driven bounds and the model-based `mussv` output.

---

### Citation

If you use this code in your work, please cite:

> M. A. Guerrero, B. Lakshminarayanan, and C. R. Rojas,  
> “Data-Driven Estimation of Structured Singular Values,”  
> *IEEE Control Systems Letters*, 2025.  
> [https://arxiv.org/abs/2503.13410](https://arxiv.org/abs/2503.13410)