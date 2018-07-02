function [ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, t ] = q1_pf_nonnonlinear_obs(mu,sigma, xtrue, ztrue, t, g, Ts, rho_0, k_rho, Q, R, N)

x_P = [];

NEES = [];
for i = 1:N
    x_P(:,i) = mu + sqrt(sigma) * randn(3,1);
end

Xpred = [mu]; % the vector of particle filter estimates.
sigma_xx = [sigma(1, 1)];
sigma_yy = [sigma(2, 2)];
sigma_tt = [sigma(3, 3)];
for j = 1:length(t)
    sigma_xx = [sigma_xx sigma(1,1)];
    sigma_yy = [sigma_yy sigma(2,2)];
    sigma_tt = [sigma_tt sigma(3,3)];
    
    %Here, we do the particle filter
    for i = 1:N

        x_P_update(:,i) = syst(x_P(1,i), x_P(2,i), x_P(3,i), g, Ts, rho_0, k_rho) + Ts*sqrt(R)*randn(3,1);

        z_update(:,i) = nonlinear_measurement_model(x_P_update(:,i));

        P_w(:,i) = (1/sqrt(det(2*pi*Q))) * exp(-(ztrue(:,j+1) - z_update(i))^2/(2*Q));
    end
    
    % Normalize to form a probability distribution (i.e. sum to 1).
    P_w = P_w./sum(P_w);

    %Resampling: From this new distribution, now we randomly sample from it to generate our new estimate particles
    

    for i = 1 : N
        x_P(:,i) = x_P_update(:, find(rand <= cumsum(P_w),1));
    end
    
    x_est = mean(x_P,2);
    sigma = cov(x_P');
    NEES = [NEES, (xtrue(:,j+1) - x_est)'*pinv(sigma)*(xtrue(:,j+1) - x_est)];
    Xpred = [Xpred x_est];
end
t = [t, t(end)+Ts];
end