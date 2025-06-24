function G_eval = evaluate_discrete_tf(G0, omega)
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
    [n_outputs, n_inputs] = size(G0);
    n_time_steps = size(u_time, 2);
    y_time = zeros(n_outputs, n_time_steps);

    for i = 1:n_outputs
        for j = 1:n_inputs
            if ~isempty(G0{i, j})
                b = G0{i, j}.b;
                a = G0{i, j}.a;
                response = filter(b, a, u_time(j, :));
                y_time(i, :) = y_time(i, :) + response + sigma * randn(size(response));
            end
        end
    end
end

function G0 = generate_random_G0(n_outputs, n_inputs, n_states, seed)
    rng(seed); 
    sysd = drss(n_states, n_outputs, n_inputs); 
    sysd = balreal(sysd);

    Gd_tf = tf(sysd);

    G0 = cell(n_outputs, n_inputs);

    for i = 1:n_outputs
        for j = 1:n_inputs
            [b, a] = tfdata(Gd_tf(i, j), 'v');

            if abs(a(1)) < 1e-12
                a(1) = 1e-12;
            end

            G0{i, j} = struct('b', b, 'a', a);
        end
    end
end


function [mu_tilde, mu_bar, B_final, W_final] = power_method_convergence(G0, B_init, W_init, num_iters, r, mm, tol, tol2, N,sigma)
  
    freqs = 2 * pi * (0:N-1) / N;
    

    B = B_init;
    W = W_init;

    s = r(1);
    f = mm(1);
    
    if s==0
        n = sum(mm(2:end));
    elseif f==0
        n=sum(r(2:end));
    else
        n = sum(r(2:end)) + sum(mm(2:end)); 
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

            Z_rj = zeros(sum(r(2:end)), 1);
            Z_mk = zeros(sum(mm(2:end)), 1);

            idx_r = 1;
            for j = 1:s
                r_j = r(j+1);
                W_rj = W(idx_r:idx_r+r_j-1, m);
                A_rj = A(idx_r:idx_r+r_j-1);
                Z_rj(idx_r:idx_r+r_j-1) = ((W_rj' * A_rj) / abs(W_rj' * A_rj)) * W_rj;
                idx_r = idx_r + r_j;
            end

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


            B_rj = zeros(sum(r(2:end)), 1);
            B_mk = zeros(sum(mm(2:end)), 1);

            idx_r = 1;
            for j = 1:s
                r_j = r(j+1);
                A_rj = A(idx_r:idx_r+r_j-1);
                W_rj = W(idx_r:idx_r+r_j-1, m);
                B_rj(idx_r:idx_r+r_j-1) = ((A_rj' * W_rj) / abs(A_rj' * W_rj)) * A_rj;
                idx_r = idx_r + r_j;
            end

            idx_mm = sum(r(2:end)) + 1;
            for k = 1:f
                mm_k = mm(k+1);
                A_mk = A(idx_mm:idx_mm+mm_k-1);
                W_mk = W(idx_mm:idx_mm+mm_k-1, m);
                B_mk(idx_mm-sum(r(2:end)):idx_mm-sum(r(2:end))+mm_k-1) = (norm(A_mk) / norm(W_mk)) * W_mk;
                idx_mm = idx_mm + mm_k;
            end
        
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