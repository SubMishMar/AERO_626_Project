clear all
close all

rho_0 = 3.4 * 10^(-3);
g = 32.2; 
k_rho = 2.2 * 10^(4);

x_0_mean = 10^(5);
x_0_var = 5 * 10^(2);

y_0_mean = -6 * 10^(3);
y_0_var = 2 * 10^(4);

z_0_mean = 2.0 * 10^(3);
z_0_var = 2.5 * 10^(5);

T = 0.1;
Tend = 17;
mu = [x_0_mean; y_0_mean; z_0_mean];
sigma = diag([x_0_var, y_0_var, z_0_var]);

R = (T^2)*diag([ 0, 2, 0]); %% Process Noise
Q_linear = 100; %% Measurement Noise
Q_nonlinear = 200; %% Measurement Noise for nonlinear case
Q_nonlinear_PF = 2 * 10^(12);%% Measurement Noise for PF with nonlinear measurement
H = [1, 0, 0];

t = 0:T:Tend;
n = 3;
m = 1;
xtrue = zeros(n,length(t)+1);
ztrue_linear = zeros(m,length(t)+1);
ztrue_nonlinear = zeros(m,length(t)+1);
xtrue(:,1) = [x_0_mean; y_0_mean; z_0_mean];

n_ensembles = 175;
N = 1000; 
%% Truth Simulation
for i = 1:length(t)
   xtrue(:,i+1) = syst(xtrue(1,i), xtrue(2, i), xtrue(3, i), g, T, rho_0, k_rho) + T*sqrt(R)*randn(n, 1);
   ztrue_linear(:,i+1) = H*xtrue(:,i+1) + sqrt(Q_linear)*randn(m,1);
   ztrue_nonlinear(:,i+1) = xtrue(1,i+1)^2 + xtrue(2,i+1)^2 + sqrt(Q_nonlinear)*randn(m,1);
end


%% Q1 EKF Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ekf_linear_obs(mu, sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H);
plotname = "Q1 EKF Linear Measurment";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 EKF Non Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ekf_nonlinear_obs(mu, sigma, xtrue, ztrue_nonlinear, t, g, T, rho_0, k_rho, Q_nonlinear, R);
plotname = "Q1 EKF Non Linear Measurment";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 UKF Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ukf_linear_obs( mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R);
plotname = "Q1 UKF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 UKF Non Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ukf_nonlinear_obs( mu,sigma, xtrue, ztrue_nonlinear, t, g, T, rho_0, k_rho, Q_nonlinear, R);
plotname = "Q1 UKF Non Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 EnKF Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_enkf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H, n_ensembles);
plotname = "Q1 EnKF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 EnKF Non Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_enkf_nonlinear_obs(mu,sigma, xtrue, ztrue_nonlinear, t, g, T, rho_0, k_rho, Q_nonlinear, R, n_ensembles);
plotname = "Q1 EnKF Non Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 PF Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_pf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, N);
plotname = "Q1 PF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 PF Non Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_pf_nonlinear_obs(mu,sigma, xtrue, ztrue_nonlinear, t, g, T, rho_0, k_rho, Q_nonlinear_PF, R, N);
plotname = "Q1 PF Non Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

