function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, t] = q1_ekf_linear_obs(mu,sigma, xtrue, ztrue, t, g, T, rho_0, k_rho, Q, R, H)

sigma_xx = [sigma(1, 1)];
sigma_yy = [sigma(2, 2)];
sigma_tt = [sigma(3, 3)];
Xpred =[mu];
mu_posterior = mu;
NEES = [];
for i = 1:length(t)
    
   sigma_xx = [sigma_xx, sigma(1, 1)];
   sigma_yy = [sigma_yy, sigma(2, 2)];
   sigma_tt = [sigma_tt, sigma(3, 3)];
   
   mu = syst(mu(1), mu(2), mu(3), g, T, rho_0, k_rho);
   G = syst_jacobian(mu_posterior(1), mu_posterior(2), mu_posterior(3), T, rho_0, k_rho);
   sigma = G*sigma*G' + R;
   K = sigma*H'/(H*sigma*H' + Q);
 
   mu = mu + K*(ztrue(:,i+1) - H*mu);
   sigma = (eye(3) - K*H)*sigma;
   NEES = [NEES, ((xtrue(:,i+1) - mu)'/(sigma))*(xtrue(:,i+1) - mu)];
   Xpred = [Xpred, mu];
   mu_posterior = mu;
   
end
t = [t, t(end)+T];
end