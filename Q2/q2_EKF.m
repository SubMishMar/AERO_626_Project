function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time] = q2_EKF(mu, sigma, xtrue, ztrue, Ts, Q, R, H, Tend)

sigma_xx = [sigma(1,1)];
sigma_yy = [sigma(2,2)];
sigma_tt = [sigma(3,3)];
count = 1;
Xpred = [mu];
NEES = [];
for i=0:Ts:Tend
   
   sigma_xx = [sigma_xx, sigma(1, 1)];
   sigma_yy = [sigma_yy, sigma(2, 2)];
   sigma_tt = [sigma_tt, sigma(3, 3)];
   
   vk = sin(i);
   if i <= 50
       wk = 0.1;
   elseif i > 50 && i <= 80
       wk = -0.2;
   else
       wk = -0.1;
   end
 
   mu = syst(mu(1), mu(2), mu(3), vk, wk, Ts);
   G = syst_jacobian(mu(1), mu(2), mu(3), vk, wk, Ts);
   sigma = G*sigma*G' + R;
   
   if mod(i, 5*Ts) == 0 && i ~= 0
    %disp('measurement taken in EKF');
    K = sigma*H'/(H*sigma*H' + Q);
    mu = mu + K*(ztrue(:,count+1) - H*mu);
    sigma = (eye(3) - K*H)*sigma;
   end
   
   NEES = [NEES, ((xtrue(:,count+1) - mu)'/sigma)*(xtrue(:,count+1) - mu)];
   Xpred = [Xpred, mu];
   count = count + 1;
end
   time = [0:Ts:Tend, Tend+Ts];
end
