% =========================================================================
% Code for: "Data-Driven Estimation of Structured Singular Values"
%
% Author: Margarita Andrea Guerrero Salazar (aka MAGS)
% Contact: mags3@kth.se
%
% If you use this code, please cite the associated paper:
%   M. A. Guerrero, B. Lakshminarayanan, C. R. Rojas,
%   "Data-Driven Estimation of Structured Singular Values",
%   https://arxiv.org/abs/2503.13410
%
% This code is provided for academic and research purposes only.
% Any reproduction, modification, or distribution should properly cite the paper above.
%
% © 2025 MAGS. All rights reserved.
% =========================================================================

seed=5; %For Paper reproducibility
rng(seed);

function G_eval = evaluate_discrete_tf(G0, omega)

% Evaluate frequency response of discrete-time transfer function matrix
%
% Args:
%   G0: Cell array (n_outputs x n_inputs) of discrete transfer functions (structs with fields 'b' and 'a')
%   omega: Frequency in radians/sample
%
% Returns:
%   G_eval: Frequency response matrix G(e^{jω}) of size (n_outputs x n_inputs)

    [n_outputs, n_inputs] = size(G0);
    G_eval = zeros(n_outputs, n_inputs);
    z = exp(1j * omega);
    for i = 1:n_outputs
        for j = 1:n_inputs
            b = G0{i, j}.b;
            a = G0{i, j}.a;
            num = sum(b .* (z .^ -(0:length(b)-1)));
            den = sum(a .* (z .^ -(0:length(a)-1)));
            G_eval(i, j) = num / den;
        end
    end
end

function y_time = simulate_plant_general_discrete(G0, u_time,sigma)

% Simulate the discrete-time response of a MIMO plant G0 with additive Gaussian noise
%
% Args:
%   G0: Cell array of transfer functions (structs with fields 'b' and 'a')
%   u_time: Input signal matrix of size (n_inputs x time_steps)
%   sigma: Standard deviation of the output noise
%
% Returns:
%   y_time: Output signal matrix of size (n_outputs x time_steps)

    [n_outputs, n_inputs] = size(G0);
    n_time_steps = size(u_time, 2);
    y_time = zeros(n_outputs, n_time_steps);

    for i = 1:n_outputs
        for j = 1:n_inputs
            if ~isempty(G0{i, j})
                b = G0{i, j}.b;
                a = G0{i, j}.a;
                response = filter(b, a, u_time(j, :)); %an~adir ruido acá
                y_time(i, :) = y_time(i, :) + response + sigma * randn(size(response));
            end
        end
    end
end

function G0 = generate_random_G0(n_outputs, n_inputs, n_states, seed)

% Generate a random, stable discrete-time MIMO system G0 as structured TFs
%
% Args:
%   n_outputs: Number of output channels
%   n_inputs: Number of input channels
%   n_states: Number of internal states of the system
%   seed: Random seed for reproducibility
%
% Returns:
%   G0: Cell array (n_outputs x n_inputs) of TF structs with fields:
%       - b: Numerator coefficients
%       - a: Denominator coefficients

    rng(seed); % Set seed for reproducibility

    % Generate a random discrete-time state-space system
    sysd = drss(n_states, n_outputs, n_inputs); % Already stable by design

    % Balance the realization to improve numerical properties
    sysd = balreal(sysd);

    % Convert to transfer function
    Gd_tf = tf(sysd);

    % Initialize the G0 cell structure
    G0 = cell(n_outputs, n_inputs);

    % Convert each transfer function entry to structured format
    for i = 1:n_outputs
        for j = 1:n_inputs
            % Extract numerator and denominator coefficients
            [b, a] = tfdata(Gd_tf(i, j), 'v');

            % Ensure denominator is nonzero
            if abs(a(1)) < 1e-12
                a(1) = 1e-12;
            end

            % Store in structured format
            G0{i, j} = struct('b', b, 'a', a);
        end
    end
end


function [mu_tilde, mu_bar, B_final, W_final] = power_method_convergence(G0, B_init, W_init, num_iters, r, mm, tol, tol2, N,sigma)
% Power method to estimate lower bound of structured singular value (μ)
%
% Args:
%   G0: Discrete-time plant in cell array TF format
%   B_init: Initial frequency-domain perturbation signal matrix (n x N)
%   W_init: Initial frequency-domain weighting signal matrix (n x N)
%   num_iters: Number of power method iterations
%   r: Vector [s, d1, ..., ds] defining scalar uncertainty blocks (s: number of blocks, d_i: sizes)
%   mm: Vector [f, d1, ..., df] defining full uncertainty blocks (f: number of blocks, d_i: sizes)
%   tol: Tolerance for convergence in μ̃ and μ̄
%   tol2: Secondary tolerance (optional, not used)
%   N: Number of frequency points (length of FFT)
%   sigma: Standard deviation of Gaussian output noise
%
% Returns:
%   mu_tilde: Matrix (N x num_iters) with μ̃ estimates per iteration and frequency
%   mu_bar: Matrix (N x num_iters) with μ̄ estimates per iteration and frequency
%   B_final: Final B matrix (n x N)
%   W_final: Final W matrix (n x N)
  
    freqs = 2 * pi * (0:N-1) / N;
    
    %[B_init, W_init] = initialize_bw_with_Dopt(G0, freqs, s, f, N);
    % for m = 1:N
    %     B_init(:, m) = B_init(:, m) / norm(B_init(:, m));
    %     W_init(:, m) = W_init(:, m) / norm(W_init(:, m));
    % end
    B = B_init;
    W = W_init;

    % Compute total dimension n
    s = r(1);  % Number of scalar blocks
    f = mm(1); % Number of full blocks
    
    if s==0
        n = sum(mm(2:end));
    elseif f==0
        n=sum(r(2:end));
    else
        n = sum(r(2:end)) + sum(mm(2:end)); % Total size of the uncertainty
    end

    num_freqs = N;
    [n_outputs, n_inputs] = size(G0);
    
    if n ~= n_outputs
        error('Mismatch: The total uncertainty dimension (%d) does not match the plant size (%d x %d)', n, n_outputs, n_inputs);
    end
    Z = zeros(size(B));

    mu_tilde = zeros(num_freqs, num_iters);
    mu_bar = zeros(num_freqs, num_iters);
    convergence_flags = false(1, num_freqs);
    

    for l = 1:num_iters
%% <<<<<<<<<  =====  First Experiment with G_0  =====    >>>>>>>>>>

        
        b_time=ifft(B,[],2,'symmetric');
        p_time = simulate_plant_general_discrete(G0, b_time,sigma);
        P_freq = fft(p_time,[],2);
%% <<<<<<<<< =====  End of First Experiment with G_0  =====    >>>>>>>>>

        for m = 1:num_freqs

            M_tilde = norm(P_freq(:, m));
           
            A = (1 / M_tilde) * P_freq(:, m);

            % Initialize structured uncertainty Z
            Z_rj = zeros(sum(r(2:end)), 1);
            Z_mk = zeros(sum(mm(2:end)), 1);

            % Process each scalar uncertainty block
            idx_r = 1;
            for j = 1:s
                r_j = r(j+1);
                W_rj = W(idx_r:idx_r+r_j-1, m);
                A_rj = A(idx_r:idx_r+r_j-1);
                Z_rj(idx_r:idx_r+r_j-1) = ((W_rj' * A_rj) / abs(W_rj' * A_rj)) * W_rj;
                idx_r = idx_r + r_j;
            end

            % Process each full uncertainty block
            idx_mm = sum(r(2:end)) + 1;
            for k = 1:f
                mm_k = mm(k+1);
                W_mk = W(idx_mm:idx_mm+mm_k-1, m);
                A_mk = A(idx_mm:idx_mm+mm_k-1);
                Z_mk(idx_mm-sum(r(2:end)):idx_mm-sum(r(2:end))+mm_k-1) = (norm(W_mk) / norm(A_mk)) * A_mk;
                idx_mm = idx_mm + mm_k;
            end

            
            if s>0 && f>0
                Z(:, m) = [Z_rj; Z_mk];
            elseif s==0
                Z(:, m) = Z_mk;
            elseif f==0
                Z(:, m) = Z_rj;
            end
          

            mu_tilde(m, l) = M_tilde;

        end

%% <<<<<<<<<<<  ===== Second Experiment with G_0^T  =====   >>>>>>>>> 

        Z_time=ifft(Z,[],2,'symmetric');
        Z_time=flip(Z_time,2);

        
        r_time = zeros(n_outputs, N);

        for alpha = 1:n_inputs
            for beta = 1:n_outputs
                e_alpha = zeros(n_inputs, 1);
                e_beta = zeros(n_outputs, 1);
                e_alpha(alpha) = 1;
                e_beta(beta) = 1;

                u_time = e_alpha * (e_beta' * Z_time);
                response = simulate_plant_general_discrete(G0, u_time,sigma);

                response = e_alpha * (e_beta' * response);
                r_time = r_time + response;
            end
        end

        r_time=flip(r_time,2);
        R_freq=fft(r_time,[],2);
%% <<<<<<<<<  ===== End of Second Experiment with G_0^T  =====   >>>>>>>>> 
        
        for m = 1:num_freqs     
            
            M_bar = norm(R_freq(:, m));
            W(:, m) = (1 / M_bar) * R_freq(:, m);

            % 
            % Initialize B blocks
            B_rj = zeros(sum(r(2:end)), 1);
            B_mk = zeros(sum(mm(2:end)), 1);

            % Process each scalar block
            idx_r = 1;
            for j = 1:s
                r_j = r(j+1);
                A_rj = A(idx_r:idx_r+r_j-1);
                W_rj = W(idx_r:idx_r+r_j-1, m);
                B_rj(idx_r:idx_r+r_j-1) = ((A_rj' * W_rj) / abs(A_rj' * W_rj)) * A_rj;
                idx_r = idx_r + r_j;
            end

            % Process each full block
            idx_mm = sum(r(2:end)) + 1;
            for k = 1:f
                mm_k = mm(k+1);
                A_mk = A(idx_mm:idx_mm+mm_k-1);
                W_mk = W(idx_mm:idx_mm+mm_k-1, m);
                B_mk(idx_mm-sum(r(2:end)):idx_mm-sum(r(2:end))+mm_k-1) = (norm(A_mk) / norm(W_mk)) * W_mk;
                idx_mm = idx_mm + mm_k;
            end
        
            % 
            if s > 0 && f > 0
                B(:, m) = [B_rj; B_mk];
            elseif s == 0
                B(:, m) = B_mk;
            elseif f == 0
                B(:, m) = B_rj;
            end
        
            mu_bar(m, l) = M_bar;
        end

        for m = 1:N
            B(:, m) = B(:, m) / norm(B(:, m));
            W(:, m) = W(:, m) / norm(W(:, m));
        end
        % norm_factor = sqrt(sum(vecnorm(B).^2)); % Sum of squared norms for all frequencies
        % B = B / norm_factor; % Normalize B
        for m = 1:num_freqs
            if l > 1
                convergence_flags(m) = abs(mu_tilde(m, l) - mu_tilde(m, l - 1)) < tol && abs(mu_bar(m, l) - mu_bar(m, l - 1)) < tol;
            end
        end
        
        if all(convergence_flags)
            mu_tilde = mu_tilde(:, 1:l);
            mu_bar = mu_bar(:, 1:l);
            break;
        end
    end

    B_final = B;
    W_final = W;
end


%% Frequency Response Case 4 - WITHOUT NOISE
% ===== Case 4: Scalar 1x1, Full 2x2 =====

n=3;
n_outputs = n; % Number of outputs
n_inputs = n;  % Number of inputs
n_states = 2;  % Number of states (adjust as needed)
seed=5;

% Generate random G0
G0 = generate_random_G0(n_outputs, n_inputs, n_states,seed);

clc;

% ===== Load System (Example Case 4: Scalar 1x1, Full 2x2) =====
[n_outputs, n_inputs] = size(G0);

% Number of states (fixed at 3)
n = 3; 

% Define the case
r = [1, 1];  % 1x1 scalar
m = [1, 2];  % 2x2 full block
mussv_block = [1, 0; 2, 2];  % MUSSV block structure

% Define the number of frequencies (N)
N = 10000;

% Generate frequency vector
freqs = 2 * pi * (0:N-1) / N; % Positive frequencies only

% ===== Generate Random B and W Matrices =====
B_freq = randn(n, N) + 1j * randn(n, N);
W_freq = randn(n, N) + 1j * randn(n, N);

% ===== Apply Symmetry on B and W for Tom Oomen Algorithm =====
for i = 1:n
    B_freq(i, [1, N/2+1]) = real(B_freq(i, [1, N/2+1]));
    W_freq(i, [1, N/2+1]) = real(W_freq(i, [1, N/2+1]));
    B_freq(i, 2:N/2) = B_freq(i, 2:N/2);
    B_freq(i, N:-1:N/2+2) = conj(B_freq(i, 2:N/2));
    W_freq(i, 2:N/2) = W_freq(i, 2:N/2);
    W_freq(i, N:-1:N/2+2) = conj(W_freq(i, 2:N/2));
end

% ===== Run the Power Method =====
[mu_tilde, mu_bar, ~, ~] = power_method_convergence(G0, B_freq, W_freq, 100, r, m, 1e-6, 1e-1, N,0);

% ===== Compute MUSSV Lower Bound =====
G_frd = zeros(n, n, N);  % Storage for frequency responses
for m_idx = 1:N
    Gf = evaluate_discrete_tf(G0, freqs(m_idx)); % Evaluate plant at each frequency
    G_frd(:, :, m_idx) = Gf;
end

% Create frequency response model
G_mussv = frd(G_frd, freqs);

% Compute MUSSV lower bound
[bounds, ~] = mussv(G_mussv, mussv_block);
mussv_lower_values = bounds(1).ResponseData(:); % Extract lower bound as vector
[max_mussv_value, max_mussv_idx] = max(mussv_lower_values);
max_mussv_freq =2*pi - freqs(max_mussv_idx);


% ===== Extract Last Iteration Values for Power Method =====
mu_tilde_final = mu_tilde(:, end);  % Last iteration of mu_tilde
mu_bar_final = mu_bar(:, end);      % Last iteration of mu_bar

% ===== Plot Results (Only Positive Frequencies) =====
% Generate mirrored frequencies from -pi to pi
freqs_mirrored =  [-flip(freqs(2:end)), freqs]; % Exclude zero to avoid duplication

% Mirror the data for symmetric plotting
mu_tilde_mirrored = [flip(mu_tilde_final(2:end)); mu_tilde_final];  
mu_bar_mirrored = [flip(mu_bar_final(2:end)); mu_bar_final];  
mussv_lower_mirrored = [flip(mussv_lower_values(2:end)); mussv_lower_values];

% Plot results
figure;
hold on;
set(gca,'Box','on')


% Plot Power Method Results with symmetry
plot(freqs_mirrored, mu_tilde_mirrored, 'b-', 'LineWidth', 0.2, 'DisplayName', '$\tilde{\mu}$');
plot(freqs_mirrored, mu_bar_mirrored, 'g-', 'LineWidth', 0.2, 'DisplayName', '$\bar{\mu}$');

% Plot MUSSV Lower Bound with symmetry
plot(freqs_mirrored, mussv_lower_mirrored, 'r:', 'LineWidth', 1, 'DisplayName', 'MUSSV');
% Set legend with LaTeX
lgd = legend('show', 'Location', 'best');
set(lgd, 'Interpreter', 'latex', 'FontSize', 11);

yline(max_mussv_value, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); % Horizontal line
xline(max_mussv_freq, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off');  % Vertical line at positive peak frequency
xline(-max_mussv_freq, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off'); % Vertical line at negative peak frequency

xlabel('$\omega$ (rad/s)', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$\mu_{\Delta}^{\prime} (G_0(e^{i\omega}))$', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'FontSize', 10);

xlim([-pi+0.005, pi-0.005]); % Ensure x-axis is from -pi to pi

hold off;
print('-depsc2', 'mu_plot2.eps');