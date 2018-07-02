function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_UKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend)

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
sigma_xx = [sigma(1,1)];
sigma_yy = [sigma(2,2)];
sigma_tt = [sigma(3,3)];
count = 1;
Xpred = [mu];
NEES = [];
for i=0:Ts:Tend

    vk = sin(i);
if i <= 50
   wk = 0.1;
elseif i > 50 && i <= 80
   wk = -0.2;
else
   wk = -0.1;
end


sigma_xx = [sigma_xx sigma(1,1)];
sigma_yy = [sigma_yy sigma(2,2)];
sigma_tt = [sigma_tt sigma(3,3)];



L = chol(sigma, 'lower');

sgma_pts(:,1) = mu;
sgma_star_bar(:,1) = syst(sgma_pts(1,1), sgma_pts(2,1), sgma_pts(3,1), vk, wk, Ts);
for j = 2:n+1
    sgma_pts(:,j) = mu + sqrt(n+lambda)*L(:, j-1);
    sgma_star_bar(:,j) = syst(sgma_pts(1,j), sgma_pts(2,j), sgma_pts(3,j),vk, wk, Ts);
end

for j = n+2:2*n+1
    sgma_pts(:,j) = mu - sqrt(n+lambda)*L(:, j-4);
    sgma_star_bar(:,j) = syst(sgma_pts(1,j), sgma_pts(2,j), sgma_pts(3,j), vk, wk, Ts);
end

mu = sum(weights_m.*sgma_star_bar,2);

sigma = zeros(n, n);

for j=1:(2*n+1)
    sigma = sigma + weights_c(j)*(sgma_star_bar(:,j)-mu)*(sgma_star_bar(:,j)-mu)';
end

sigma = sigma + R;
if mod(i, 5*Ts) == 0 && i ~= 0
    % disp('measurement taken in UKF');
%correction begins
L = chol(sigma, 'lower');

Z_bar = zeros(2, 2*n+1);

sgma_pts_bar(:,1) = mu;
Z_bar(:,1) = H*sgma_pts_bar(:,1);
for j = 2:n+1
    sgma_pts_bar(:,j) = mu + sqrt(n+lambda)*L(:, j-1);
    Z_bar(:,j) = H * sgma_pts_bar(:,j);
end

for j = n+2:2*n+1
    sgma_pts_bar(:,j) = mu - sqrt(n+lambda)*L(:, j-4);
    Z_bar(:,j) = H * sgma_pts_bar(:,j);
end


z_hat = sum(weights_m.* Z_bar, 2);

S = zeros(size(Q));
for j=1:(2*n+1)
    S = S + weights_c(j)*(Z_bar(j)-z_hat)*(Z_bar(j)-z_hat)';
end
S = S + Q;

sigma_bar_xz = zeros(n, 1);
for j=1:(2*n+1)
    sigma_bar_xz = sigma_bar_xz + weights_c(j)*(sgma_pts_bar(:,j) - mu)*(Z_bar(:,j)-z_hat)';
end

K = sigma_bar_xz / S;


mu = mu + K*(ztrue(:,count+1) - z_hat);
sigma = sigma - K*S*K';
% Correction Ends
end
NEES = [NEES, ((xtrue(:,count+1) - mu)'/sigma)*(xtrue(:,count+1) - mu)];
Xpred = [Xpred, mu];
count = count + 1;
end
time = [0:Ts:Tend, Tend+Ts];
end