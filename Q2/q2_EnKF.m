function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_EnKF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, n_ensembles)

sigma_xx = [sigma(1,1)];
sigma_yy = [sigma(2,2)];
sigma_tt = [sigma(3,3)];

n = 3;
m = 2;
% intial step of generating ensembles.
X = zeros(3, n_ensembles);

for i = 1:n_ensembles
   X(:,i) = mu + sqrt(sigma)*randn(n, 1);
end
Xpred = [mu];
NEES = [];
count = 1;
for i=0:Ts:Tend
sigma_xx = [sigma_xx sigma(1,1)];
sigma_yy = [sigma_yy sigma(2,2)];
sigma_tt = [sigma_tt sigma(3,3)];

vk = sin(i);
if i <= 50
   wk = 0.1;
elseif i > 50 && i <= 80
   wk = -0.2;
else
   wk = -0.1;
end
if ((mod(i, 5*Ts) ~= 0) || i == 0)
 
for j = 1:n_ensembles
   X(:,j) = syst(X(1,j), X(2,j), X(3,j), vk, wk, Ts) + sqrt(R)*randn(n, 1);
end
xa = mean(X,2);
Exx = X - xa;
Pxx = cov(X');
sigma = Pxx;
else
%disp('measurement taken in EnKF');
for j = 1:n_ensembles
   X(:,j) = syst(X(1,j), X(2,j), X(3,j), vk, wk, Ts) + sqrt(R)*randn(n, 1);
   Y(:,j) = measurement_model(X(:,j)) + sqrt(Q)*rand(m,1);
end
xa = mean(X,2);
Exx = X - xa;
ybar = mean(Y,2);
Eyy = Y - ybar;
Pxy = (1/(n_ensembles-1))*Exx*Eyy';
Pxx = cov(X');
Pyy = cov(Y');
K = Pxy/Pyy;
for j = 1:n_ensembles
    X(:,j) = X(:,j) + K*(ztrue(:,count+1) - measurement_model(X(:,j)));
end
xa = mean(X,2);
sigma = cov(X');
%sigma = Pxx - (Pxy/Pyy)*Pxy';
end

NEES = [NEES, ((xtrue(:,count+1) - xa)'/sigma)*(xtrue(:,count+1) - xa)];
Xpred = [Xpred xa];
count = count + 1;
end
time = [0:Ts:Tend, Tend+Ts];
end