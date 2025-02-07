function [Xoutput, Xhatoutput, noise] = LEO_3_Particle(std_n, noise_switch)
    rng('default')
    %% SYSTEM
    step = 0.0001;   % sampling frequency
    final = 5; % total running time
    t = 0:step:(final - 1);
    
    d = 0.5;
    lambda = 1;
    M = 40;
    N = 40;
    
    std_w = 0.1;
    if noise_switch == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Gassian noise     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        w = normrnd(0, std_w, 4, (final - 1) / step + 1);
        noise = normrnd(0, std_n, 2, (final - 1)/step + 1);
    elseif noise_switch == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Composite wave   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define sine waves with different frequencies
        f1 = 0.5;  % 5Hz
        f2 = 1; % 15Hz
        f3 = 3; % 30Hz
        
        % composition of sine waves
        w = normrnd(0, std_w, 4, (final - 1) / step + 1);
        noise = normrnd(0, std_n, 2, (final - 1)/step + 1) + (0.05*sin(2*pi*f1*t) + 0.03*sin(2*pi*f2*t)) .* [1; 1];
    elseif noise_switch == 2   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Pulse noise    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate Bernoulli-Gaussian impulse noise
        p = 0.3;
        b = rand(1, (final - 1) / step + 1) < p; % Bernoulli random variable
        g = std_n * randn(2, (final - 1) / step + 1); % Gaussian noise
        w = normrnd(0, std_w, 4, (final - 1) / step + 1);
        noise = b .* g; % Combine Bernoulli and Gaussian
    elseif noise_switch == 3 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Rician noise    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r = makedist('Rician', 's', 1, 'sigma', std_n);
        w = normrnd(0, std_w, 4, (final - 1) / step + 1);
        noise = random(r, 2, (final - 1) / step + 1);
    end 
    
    % c
    function C = calcC(x, xhat, d, lambda, M, N)
        alphar = x(1);
        alphai = x(2);
        theta = x(3);
        phi = x(4);
        thetahat = xhat(3);
        phihat = xhat(4);
        psix = 2 * pi * d / lambda * cos(theta) * cos(phi);
        psiy = 2 * pi * d / lambda * cos(theta) * sin(phi);
        psihatx = 2 * pi * d / lambda * cos(thetahat) * cos(phihat);
        psihaty = 2 * pi * d / lambda * cos(thetahat) * sin(phihat);
        
        dyr_alphar = 0;
        for n = 1 : N
            for m = 1 : M
                dyr_alphar = dyr_alphar + cos((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
    
        dyr_alphai = 0;
        for n = 1 : N
            for m = 1 : M
                if (n == 1) && (m == 1)
                    continue;
                end
                dyr_alphai = dyr_alphai - sin((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
    
        dyi_alphar = 0;
        for n = 1 : N
            for m = 1 : M
                if (n == 1) && (m == 1)
                    continue;
                end
                dyi_alphar = dyi_alphar + sin((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
    
        dyi_alphai = 0;
        for n = 1 : N
            for m = 1 : M
                dyi_alphai = dyi_alphai + cos((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
    
        function dyr_ = calcdyr_(alphar, alphai, psix, psiy, ...
                psihatx, psihaty, dpsix_, dpsiy_, M, N)
            dyr__r = 0;
            for n = 1 : N
                for m = 1 : M
                    if (n == 1) && (m == 1)
                        continue;
                    end
                    dyr__r = dyr__r - sin((n - 1) * (psiy - psihaty) + ...
                        (m - 1) * (psix - psihatx)) * ((n - 1) * dpsiy_ + (m - 1) * dpsix_);
                end
            end
            dyr__i = 0;
            for n = 1 : N
                for m = 1 : M
                    if (n == 1) && (m == 1)
                        continue;
                    end
                    dyr__i = dyr__i + cos((n - 1) * (psiy - psihaty) + ...
                        (m - 1) * (psix - psihatx)) * ((n - 1) * dpsiy_ + (m - 1) * dpsix_);
                end
            end
            dyr_ = alphar * dyr__r - alphai * dyr__i;
        end
        
        function dyi_ = calcdyi_(alphar, alphai, psix, psiy, ...
                psihatx, psihaty, dpsix_, dpsiy_, M, N)
            dyi__r = 0;
            for n = 1 : N
                for m = 1 : M
                    if (n == 1) && (m == 1)
                        continue;
                    end
                    dyi__r = dyi__r + cos((n - 1) * (psiy - psihaty) + ...
                        (m - 1) * (psix - psihatx)) * ((n - 1) * dpsiy_ + (m - 1) * dpsix_);
                end
            end
            dyi__i = 0;
            for n = 1 : N
                for m = 1 : M
                    if (n == 1) && (m == 1)
                        continue;
                    end
                    dyi__i = dyi__i - sin((n - 1) * (psiy - psihaty) + ...
                        (m - 1) * (psix - psihatx)) * ((n - 1) * dpsiy_ + (m - 1) * dpsix_);
                end
            end
            dyi_ = alphar * dyi__r + alphai * dyi__i;
        end
        
        dpsix_theta = - 2 * pi * d / lambda * sin(theta) * cos(phi);
        dpsix_phi = - 2 * pi * d / lambda * cos(theta) * sin(phi);
        dpsiy_theta = - 2 * pi * d / lambda * sin(theta) * sin(phi);
        dpsiy_phi = 2 * pi * d / lambda * cos(theta) * cos(phi);
    
        dyr_theta = calcdyr_(alphar, alphai, psix, psiy, psihatx, psihaty, ...
            dpsix_theta, dpsiy_theta, M, N);
        dyr_phi = calcdyr_(alphar, alphai, psix, psiy, psihatx, psihaty, ...
            dpsix_phi, dpsiy_phi, M, N);
        dyi_theta = calcdyi_(alphar, alphai, psix, psiy, psihatx, psihaty, ...
            dpsix_theta, dpsiy_theta, M, N);
        dyi_phi = calcdyi_(alphar, alphai, psix, psiy, psihatx, psihaty, ...
            dpsix_phi, dpsiy_phi, M, N);
    
        C = [dyr_alphar, dyr_alphai, dyr_theta, dyr_phi;
             dyi_alphar, dyi_alphai, dyi_theta, dyi_phi];
    end
    
    function yk1 = received_signal(xk1, xhatk, d, lambda, M, N)
        alphar = xk1(1);
        alphai = xk1(2);
        theta = xk1(3);
        phi = xk1(4);
        thetahat = xhatk(3);
        phihat = xhatk(4);
        psix = 2 * pi * d / lambda * cos(theta) * cos(phi);
        psiy = 2 * pi * d / lambda * cos(theta) * sin(phi);
        psihatx = 2 * pi * d / lambda * cos(thetahat) * cos(phihat);
        psihaty = 2 * pi * d / lambda * cos(thetahat) * sin(phihat);
        
        yr_r = 0;
        for n = 1 : N
            for m = 1 : M
                yr_r = yr_r + cos((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
        yr_i = 0;
        for n = 1 : N
            for m = 1 : M
                yr_i = yr_i + sin((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
        yr = alphar * yr_r - alphai * yr_i;
        
        yi_r = 0;
        for n = 1 : N
            for m = 1 : M
                yi_r = yi_r + sin((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
        yi_i = 0;
        for n = 1 : N
            for m = 1 : M
                yi_i = yi_i + cos((n - 1) * (psiy - psihaty) + (m - 1) * (psix - psihatx));
            end
        end
        yi = alphar * yi_r + alphai * yi_i;
        yk1 = [yr; yi];
    end
    
    A1_alphar_k1 = 0.61 ; A1_alphar_k = 1 - A1_alphar_k1;
    A1_alphai_k1 = 0.59 ; A1_alphai_k = 1 - A1_alphai_k1;
    A1_theta_k1 = 0.5     ; A1_theta_k = 1 - A1_theta_k1;
    A1_phi_k1 = 0.5       ; A1_phi_k = 1 - A1_phi_k1;
    A = [A1_alphar_k1, 0, 0, 0, A1_alphar_k, 0, 0, 0;
         0, A1_alphai_k1, 0, 0, 0, A1_alphai_k, 0, 0;
         0, 0, A1_theta_k1, 0, 0, 0, A1_theta_k, 0;
         0, 0, 0, A1_phi_k1, 0, 0, 0, A1_phi_k;
         0.97, 0, 0, 0, 0.02, 0, 0, 0;
         0, 0.97, 0, 0, 0, 0.03, 0, 0;
         0, 0, 0.999, 0, 0, 0, 0.001, 0;
         0, 0, 0, 0.9997, 0, 0, 0, 0.0002];
    X0 = [1/(2^(0.5)); 1/(2^(0.5)); 1; 1; zeros(4, 1)];
    N = 10; % number of particle
    wk = 1 / N .* ones(8, N); % weight
    std_n = 0.1;
    sigma = std_n^2 .* eye(2);
    
    X = zeros(8, (final - 1) / step + 1);
    Xhat = zeros(8, (final - 1) / step + 1);
    X_ptcl = zeros(8, N, (final - 1) / step + 1);
    X_ptcl(:, :, 1) = X0 .* ones(1, N);
    weight = zeros(N, (final - 1) / step + 1);
    weight_norm = zeros(N, (final - 1) / step + 1);
    yk1 = zeros(2, (final - 1) / step + 1);

    wbar = [w; zeros(4, (final - 1) / step + 1)];

    for index = 1 : length(t)
        
        if index ~= (final - 1) / step + 1
            X(:, index + 1) = A * X(:, index) + wbar(:, index);
            
            Ck1 = calcC(X(:, index + 1), Xhat(:, index), d, lambda, M, N);
            Cbark1 = [Ck1, zeros(2, 4)];
            % yk1 = received_signal(X(:, index + 1), Xhat(:, index), d, lambda, M, N) + noise(:, index);
            yk1 = Cbark1 * X(:, index + 1) + noise(:, index);
            
            for i = 1 : N
                % generate particles
                X_ptcl(:, i, index + 1) = normrnd(A * X_ptcl(:, i, index), std_w, 8, 1);
                pyk1_xk1(i) = -0.5 * ((yk1 - Cbark1 * X_ptcl(:, i, index + 1)).' * inv(sigma) * (yk1 - Cbark1 * X_ptcl(:, i, index + 1)));
            end
            for i = 1 : N
                weight(i, index + 1) = exp(pyk1_xk1(i) - max(pyk1_xk1));
            end
            
            % normalize weights
            wk1_sum = sum(weight(:, index + 1), 1);
            for i = 1 : N
                weight_norm(i, index + 1) = weight(i, index + 1) / wk1_sum;
            end
    
            temp_ptcl = zeros(8, N);
            % state estimation
            for i = 1 : N
                Xhat(:, index + 1) = Xhat(:, index + 1) + weight_norm(i, index + 1) .* X_ptcl(:, i, index + 1);
                temp_ptcl(:, i) = X_ptcl(:, i, index + 1);
            end
    
            % resampling (residual resampling)
            Nres = N;
            for i = 1 : N
                Nres = Nres - floor(N * weight_norm(i, index + 1));
            end
            U_tilde = ones(1, Nres + 1);
            for i = 1 : Nres
                U = unifrnd(0, 1);
                U_tilde(i + 1) = U ^ (1 / (N - i + 1)) * U_tilde(i); 
                Nres_m = 0;
                for m = 1 : N
                    Nres_m = Nres_m + (N * weight_norm(m, index + 1) - floor(N * weight_norm(m, index + 1))) / Nres;
                    if U_tilde(i + 1) > Nres_m
                        X_ptcl(:, i, index + 1) = temp_ptcl(:, m + 1);
                    end
                end
            end
            for i = 1 : N
                if floor(N * weight_norm(i, index + 1)) >= 1
                    for m = Nres + 1 : Nres + floor(N * weight_norm(i, index + 1))
                        X_ptcl(:, m, index + 1) = temp_ptcl(:, i);
                    end
                    Nres = Nres + floor(N * weight_norm(i, index + 1));
                end
            end
        end
    end
    Xoutput = X(3 : 4, :);
    Xhatoutput = Xhat(3 : 4, :);
end
