close all
clear all
clc

x_0_mean = 0;
x_0_var = 0.1;

y_0_mean = 0;
y_0_var = 0.1;

theta_0_mean = pi/2;
theta_0_var = 0.4;

Ts = 0.1;
Tend = 100;
mu = [x_0_mean; y_0_mean; theta_0_mean];
sigma = diag([x_0_var, y_0_var, theta_0_var]);
R = (Ts^2)*diag([0.01, 0.01, 0.2]);
Q = diag([0.2, 0.2]);
H = [1, 0, 0; 0, 1, 0];

t = 0:Ts:Tend;
n = 3;
m = 2;
xtrue = zeros(n,length(t)+1);

xtrue(:,1) = [x_0_mean; y_0_mean; theta_0_mean];

n_ensembles = 75;
N = 100; 
NEES_KF = zeros(size(t));
NEES_UKF = zeros(size(t));
NEES_EnKF = zeros(size(t));
NEES_PF = zeros(size(t));
%% For EKF
tic;
for ct = 1:50
% Truth Simulation
count = 1;
for i=0:Ts:Tend

   vk = sin(i);
   if i <= 50
       wk = 0.1;
   elseif i > 50 && i <= 80
       wk = -0.2;
   else
       wk = -0.1;
   end
   
   xtrue(:,count+1) = syst(xtrue(1,count), xtrue(2,count), xtrue(3,count), vk, wk, Ts) + Ts*sqrt(R)*randn(3,1);
   ztrue(:,count+1) = H*xtrue(:,count+1) + sqrt(Q)*randn(2,1);
   count = count + 1;
   
end
% Filter
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_ekf, time] = q2_EKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend);
NEES_KF = NEES_KF + NEES_ekf;
errors(:,:,ct) = (xtrue - Xpred).^2;
end
NEES_KF = NEES_KF/ct;
avg_RMSE_EKF = mean(sqrt(mean(errors,2)),3);
timervalue_ekf = toc/count;
%% For UKF
tic
for ct = 1:50
% Truth Simulation
count = 1;
for i=0:Ts:Tend

   vk = sin(i);
   if i <= 50
       wk = 0.1;
   elseif i > 50 && i <= 80
       wk = -0.2;
   else
       wk = -0.1;
   end
   
   xtrue(:,count+1) = syst(xtrue(1,count), xtrue(2,count), xtrue(3,count), vk, wk, Ts) + Ts*sqrt(R)*randn(3,1);
   ztrue(:,count+1) = H*xtrue(:,count+1) + sqrt(Q)*randn(2,1);
   count = count + 1;
   
end
% Filter
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_ukf, time] = q2_UKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend);
NEES_UKF = NEES_UKF + NEES_ukf;
errors(:,:,ct) = (xtrue - Xpred).^2;
end
NEES_UKF = NEES_UKF/ct;
avg_RMSE_UKF = mean(sqrt(mean(errors,2)),3);
timervalue_ukf = toc/count;
%% For EnKF
tic;
for ct = 1:50
% Truth Simulation
count = 1;
for i=0:Ts:Tend

   vk = sin(i);
   if i <= 50
       wk = 0.1;
   elseif i > 50 && i <= 80
       wk = -0.2;
   else
       wk = -0.1;
   end
   
   xtrue(:,count+1) = syst(xtrue(1,count), xtrue(2,count), xtrue(3,count), vk, wk, Ts) + Ts*sqrt(R)*randn(3,1);
   ztrue(:,count+1) = H*xtrue(:,count+1) + sqrt(Q)*randn(2,1);
   count = count + 1;
   
end
% Filter
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_enkf, time] = q2_EnKF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, n_ensembles);
NEES_EnKF = NEES_EnKF + NEES_enkf;
errors(:,:,ct) = (xtrue - Xpred).^2;
end
NEES_EnKF = NEES_EnKF/ct;
avg_RMSE_EnKF = mean(sqrt(mean(errors,2)),3);
timervalue_enkf = toc/count;
%% For PF
tic;
for ct = 1:50
% Truth Simulation
count = 1;
for i=0:Ts:Tend

   vk = sin(i);
   if i <= 50
       wk = 0.1;
   elseif i > 50 && i <= 80
       wk = -0.2;
   else
       wk = -0.1;
   end
   
   xtrue(:,count+1) = syst(xtrue(1,count), xtrue(2,count), xtrue(3,count), vk, wk, Ts) + Ts*sqrt(R)*randn(3,1);
   ztrue(:,count+1) = H*xtrue(:,count+1) + sqrt(Q)*randn(2,1);
   count = count + 1;
   
end
% Filter
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES_pf, time] = q2_PF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, N);
NEES_PF = NEES_PF + NEES_pf;
errors(:,:,ct) = (xtrue - Xpred).^2;
end
NEES_PF = NEES_PF/ct;
avg_RMSE_PF = mean(sqrt(mean(errors,2)),3);
timervalue_pf = toc/count;
%% bar plots for avg error and execution time comparison

figure(1);
Y = [avg_RMSE_EKF'; avg_RMSE_UKF'; avg_RMSE_EnKF'; avg_RMSE_PF'];
bar(Y);
grid;
text([0.7825, 1, 1.225],avg_RMSE_EKF, num2str(avg_RMSE_EKF), 'vert','bottom','horiz','center');
text([1.7825, 2, 2.225],avg_RMSE_UKF, num2str(avg_RMSE_UKF), 'vert','bottom','horiz','center');
text([2.7825, 3, 3.225],avg_RMSE_EnKF, num2str(avg_RMSE_EnKF), 'vert','bottom','horiz','center');
text([3.7825, 4, 4.225],avg_RMSE_PF, num2str(avg_RMSE_PF), 'vert','bottom','horiz','center');
u = legend('$x_{1}$','$x_{2}$', '$x_{3}$','Location','northwest');
set(u,'Interpreter','latex');
title('Average RMSE for 1: EKF, 2: UKF, 3: EnKF, 4: PF');

figure(2)
Y = [timervalue_ekf, timervalue_ukf, timervalue_enkf, timervalue_pf];
bar(Y);
grid;
text(1:length(Y),Y,num2str(Y'),'vert','bottom','horiz','center'); 
box off
title('Average Execution time for 1 iteration in seconds for 1: EKF, 2: UKF, 3: EnKF, 4: PF');

%% NEES plots
upperbound =(chi2inv(0.95, 150)/50)*ones(size(NEES_PF));


figure(3)
subplot(221)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_KF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('Average NEES value for EKF');

subplot(222)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_UKF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('Average NEES value for UKF');

subplot(223)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_EnKF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('Average NEES value for EnKF');

subplot(224)
plot(t, upperbound, '--r');
hold on;
plot(t, NEES_PF);
hold off;
grid;
xlabel('time(s)')
ylabel('Average NEES value');
title('Average NEES value for PF');