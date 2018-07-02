function [ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, t ] = q1_enkf_linear_obs(mu,sigma, xtrue, ztrue, t, g, T, rho_0, k_rho, Q, R, H, n_ensembles)

X = zeros(3, n_ensembles);
for i = 1:n_ensembles
   X(:,i)= mu + [sqrt(sigma(1,1))*randn(1,1); sqrt(sigma(2,2))*randn(1,1); sqrt(sigma(3,3))*randn(1,1)];
end

sigma_xx = [sigma(1, 1)];
sigma_yy = [sigma(2, 2)];
sigma_tt = [sigma(3, 3)];
Xpred =[mu];
NEES = [];
for i = 1:length(t)
    
sigma_xx = [sigma_xx, sigma(1, 1)];
sigma_yy = [sigma_yy, sigma(2, 2)];
sigma_tt = [sigma_tt, sigma(3, 3)];

for j = 1:n_ensembles
    X(:,j) = syst(X(1,j), X(2,j), X(3,j), g, T, rho_0, k_rho) + T*sqrt(R)*randn(3,1);
end

Pf = cov(X');
K = Pf*H'/(H*Pf*H'+Q);
for j = 1:n_ensembles
    X(:,j) = X(:,j) + K*(ztrue(:,i+1) - linear_measurement_model(X(:,j)));
end
xa = mean(X,2);
Pa = cov(X');
sigma = Pa;
NEES = [NEES, (xtrue(:,i+1) - xa)'*pinv(sigma)*(xtrue(:,i+1) - xa)];
Xpred = [Xpred xa];
end
t = [t, t(end)+T];
end

