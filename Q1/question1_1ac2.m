%% Description
% This script generates plots of estimates of the states, the ground truth, errors and the
% estimation error standard deviations for first 5 seconds

%% System Initialization
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
Tend = input('Enter the execution time(in seconds): ');
mu = [x_0_mean; y_0_mean; z_0_mean];
sigma = diag([x_0_var, y_0_var, z_0_var]);

R = (T^2)*diag([ 0, 2, 0]); %% Process Noise
Q_linear = 100; %% Measurement Noise
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
end


%% Q1 EKF Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ekf_linear_obs(mu, sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H);
plotname = "Q1 EKF Linear Measurment";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 UKF Linear Measurement
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q1_ukf_linear_obs( mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R);
plotname = "Q1 UKF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 EnKF Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_enkf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H, n_ensembles);
plotname = "Q1 EnKF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);

%% Q1 PF Linear Measurement
[ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time ] = q1_pf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, N);
plotname = "Q1 PF Linear Measurement";
plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname);


