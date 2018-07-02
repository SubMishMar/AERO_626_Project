close all
clear all
clc

x_0_mean = 0;
x_0_var = 0.1;

y_0_mean = 0;
y_0_var = 0.1;

theta_0_mean = pi/2;
theta_0_var = 0.4;

b4orAftr = input('Enter 0 for before 1st measurement and 1 for after 1st measurement: ');

if b4orAftr == 0
    n = 4;
elseif b4orAftr == 1
    n = 5;
else
    disp('Wrong choice!');
    return;
end
Ts = 0.01;
Tend = n*Ts;
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
N = 1000; 
%% Truth Simulation
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


%% Q2 EKF 
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_EKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend);
figure(1)
p = plot_ellipse(Xpred, sigma_xx, sigma_yy, 1, 'r');
%% Q2 UKF
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_UKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend);
q = plot_ellipse(Xpred, sigma_xx, sigma_yy, 1, 'g');

%% Q2 EnKF
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_EnKF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, n_ensembles);
r = plot_ellipse(Xpred, sigma_xx, sigma_yy, 1, 'b');

%% Q2 PF
[Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time, x_P] = q2e_PF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, N);
s = plot_ellipse(Xpred, sigma_xx, sigma_yy, 1, 'c');
plot(x_P(1,:), x_P(2,:), 'c.');
hold off;
grid;

u = legend([p q r s ],'EKF','UKF','EnKF', 'PF');
set(u,'Interpreter','latex');
axis square
xlabel('x[m]');
ylabel('y[m]');
title('Probability Distribution of Estimates predicted by various Estimators');
