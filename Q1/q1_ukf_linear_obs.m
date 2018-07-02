function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, t] = q1_ukf_linear_obs( mu,sigma, xtrue, ztrue, t, g, T, rho_0, k_rho, Q, R)
alpha = 0.75;
beta = 2;
k = 3;
n = 3;
lambda = (alpha^2)*(n+k) - n;
sgma_pts = zeros(n, 2*n+1);

w0m = lambda/(n + lambda);
w0c = w0m + (1-alpha^2 + beta);

wim = 1/(2*(n+lambda));
wic = wim;

weights_m = [w0m, wim*ones(1, 2*n)];
weights_c = [w0c, wic*ones(1, 2*n)];

sigma_xx = [sigma(1, 1)];
sigma_yy = [sigma(2, 2)];
sigma_tt = [sigma(3, 3)];
Xpred =[mu];
NEES = [];
for i = 1:length(t)

sigma_xx = [sigma_xx, sigma(1, 1)];
sigma_yy = [sigma_yy, sigma(2, 2)];
sigma_tt = [sigma_tt, sigma(3, 3)];
L = chol(sigma, 'lower');

sgma_pts(:,1) = mu;
sgma_star_bar(:,1) = syst(sgma_pts(1,1), sgma_pts(2,1), sgma_pts(3,1), g, T, rho_0, k_rho);
for j = 2:n+1
    sgma_pts(:,j) = mu + sqrt(n+lambda)*L(:, j-1);
    sgma_star_bar(:,j) = syst(sgma_pts(1,j), sgma_pts(2,j), sgma_pts(3,j), g, T, rho_0, k_rho);
end

for j = n+2:2*n+1
    sgma_pts(:,j) = mu - sqrt(n+lambda)*L(:, j-4);
    sgma_star_bar(:,j) = syst(sgma_pts(1,j), sgma_pts(2,j), sgma_pts(3,j), g, T, rho_0, k_rho);
end

mu_bar = sum(weights_m.*sgma_star_bar,2);

sigma_bar = zeros(n, n);

for j=1:(2*n+1)
    sigma_bar = sigma_bar + weights_c(j)*(sgma_star_bar(:,j)-mu_bar)*(sgma_star_bar(:,j)-mu_bar)';
end

sigma_bar = sigma_bar + R;

L = chol(sigma_bar, 'lower');

Z_bar = zeros(1, 2*n+1);

sgma_pts_bar(:,1) = mu_bar;
Z_bar(1) = linear_measurement_model(sgma_pts_bar(:,1));
for j = 2:n+1
    sgma_pts_bar(:,j) = mu_bar + sqrt(n+lambda)*L(:, j-1);
    Z_bar(j) = linear_measurement_model(sgma_pts_bar(:,j));
end

for j = n+2:2*n+1
    sgma_pts_bar(:,j) = mu_bar - sqrt(n+lambda)*L(:, j-4);
    Z_bar(j) = linear_measurement_model(sgma_pts_bar(:,j));
end


z_hat = sum(weights_m.* Z_bar);

S = 0;
for j=1:(2*n+1)
    S = S + weights_c(j)*(Z_bar(j)-z_hat)*(Z_bar(j)-z_hat)';
end
S = S + Q;

sigma_bar_xz = zeros(n, 1);
for j=1:(2*n+1)
    sigma_bar_xz = sigma_bar_xz + weights_c(j)*(sgma_pts_bar(:,j) - mu_bar)*(Z_bar(j)-z_hat)';
end

K = sigma_bar_xz / S;

mu = mu_bar + K*(ztrue(:,i+1) - z_hat);
sigma = sigma_bar - K*S*K';
NEES = [NEES, ((xtrue(:,i+1) - mu)'/(sigma))*(xtrue(:,i+1) - mu)];
Xpred = [Xpred, mu];

end
t = [t, t(end)+T];
end

