%% Description
% This script generates plots of average estimation errors over a number of
% simulations(50)

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

NEES_KF = zeros(1, length(t));
NEES_UKF = zeros(1, length(t));
NEES_EnKF = zeros(1, length(t));
NEES_PF = zeros(1, length(t));

%% Q1 EKF Linear Measurement
tic
errors = zeros(n, length(t)+1, 50);
for count = 1:50
    % Truth Simulation
    for i = 1:length(t)
        xtrue(:,i+1) = syst(xtrue(1,i), xtrue(2, i), xtrue(3, i), g, T, rho_0, k_rho) + T*sqrt(R)*randn(n, 1);
        ztrue_linear(:,i+1) = H*xtrue(:,i+1) + sqrt(Q_linear)*randn(m,1);
    end
    % Filter
    [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_kf, time] = q1_ekf_linear_obs(mu, sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H);
    NEES_KF = NEES_KF + NEES_kf;
    errors(:, :, count) = abs(Xpred - xtrue);
end
mean_error_ekf = mean(mean(errors, 2),3);
timervalue_ekf = toc/count;
%% Q1 UKF Linear Measurement
tic
errors = zeros(n, length(t)+1, 50);
for count = 1:50
    % Truth Simulation
    for i = 1:length(t)
        xtrue(:,i+1) = syst(xtrue(1,i), xtrue(2, i), xtrue(3, i), g, T, rho_0, k_rho) + T*sqrt(R)*randn(n, 1);
        ztrue_linear(:,i+1) = H*xtrue(:,i+1) + sqrt(Q_linear)*randn(m,1);
    end
    % Filter
    [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_ukf, time] = q1_ukf_linear_obs( mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R);
    NEES_UKF = NEES_UKF + NEES_ukf; 
    errors(:, :, count) = abs(Xpred - xtrue);
end
mean_error_ukf = mean(mean(errors, 2),3);
timervalue_ukf = toc/count;
NEES_UKF = NEES_UKF/count;

%% Q1 EnKF Linear Measurement
tic
errors = zeros(n, length(t)+1, 50);
for count = 1:50
    % Truth Simulation
    for i = 1:length(t)
        xtrue(:,i+1) = syst(xtrue(1,i), xtrue(2, i), xtrue(3, i), g, T, rho_0, k_rho) + T*sqrt(R)*randn(n, 1);
        ztrue_linear(:,i+1) = H*xtrue(:,i+1) + sqrt(Q_linear)*randn(m,1);
    end
    % Filters
    [ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_enkf, time ] = q1_enkf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, H, n_ensembles);
    NEES_EnKF = NEES_EnKF + NEES_enkf; 
    errors(:, :, count) = abs(Xpred - xtrue);
end
mean_error_enkf = mean(mean(errors, 2),3);
timervalue_enkf = toc/count;
NEES_EnKF = NEES_EnKF/count;
%% Q1 PF Linear Measurement
tic
errors = [];
for count = 1:50
    % Truth Simulation
    for i = 1:length(t)
        xtrue(:,i+1) = syst(xtrue(1,i), xtrue(2, i), xtrue(3, i), g, T, rho_0, k_rho) + T*sqrt(R)*randn(n, 1);
        ztrue_linear(:,i+1) = H*xtrue(:,i+1) + sqrt(Q_linear)*randn(m,1);
    end
    % Filters
    [ Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_pf, time ] = q1_pf_linear_obs(mu,sigma, xtrue, ztrue_linear, t, g, T, rho_0, k_rho, Q_linear, R, N);
    if length(Xpred) == length(xtrue)
        NEES_PF = NEES_PF + NEES_pf;
        errors(:, :, count) = abs(Xpred - xtrue);
    else
        break;
    end
end
mean_error_pf = mean(mean(errors, 2),3);
timervalue_pf = toc/count;
NEES_PF = NEES_PF/count;
%% bar plots for avg error and execution time comparison


figure(1);
Y = [mean_error_ekf'; mean_error_ukf'; mean_error_enkf'; mean_error_pf'];
bar(Y);
grid;
text([0.7825, 1, 1.225], mean_error_ekf, num2str(mean_error_ekf), 'vert','bottom','horiz','center');
text([1.7825, 2, 2.225], mean_error_ukf, num2str(mean_error_ukf), 'vert','bottom','horiz','center');
text([2.7825, 3, 3.225], mean_error_enkf, num2str(mean_error_enkf), 'vert','bottom','horiz','center');
text([3.7825, 4, 4.225], mean_error_pf, num2str(mean_error_pf), 'vert','bottom','horiz','center');
u = legend('$x_{1}$','$x_{2}$', '$x_{3}$','Location','northwest');
set(u,'Interpreter','latex');
title('Absolute Average State Estimation Errors for 1: EKF, 2: UKF, 3: EnKF, 4: PF');

figure(2)
Y = [timervalue_ekf, timervalue_ukf, timervalue_enkf, timervalue_pf];
bar(Y);
grid;
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center'); 
box off
title('Average Execution time for 1 loop in seconds for 1: EKF, 2: UKF, 3: EnKF, 4: PF');

%% NEES plots
upperbound = (chi2inv(0.95,150)/50)*ones(1, length(NEES_KF));

figure(3)
subplot(221)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_KF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('NEES for EKF');

subplot(222)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_UKF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('NEES for UKF');

subplot(223)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_EnKF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('NEES for EnKF');

subplot(224)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_PF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('NEES for PF');