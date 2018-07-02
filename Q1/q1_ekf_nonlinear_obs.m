function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, t] = q1_ekf_nonlinear_obs(mu,sigma, xtrue, ztrue, t, g, T, rho_0, k_rho, Q, R)
NEES = [];
sigma_xx = [sigma(1, 1)];
sigma_yy = [sigma(2, 2)];
sigma_tt = [sigma(3, 3)];
Xpred =[mu];
mu_posterior = mu;
for i = 1:length(t)
    
   sigma_xx = [sigma_xx, sigma(1, 1)];
   sigma_yy = [sigma_yy, sigma(2, 2)];
   sigma_tt = [sigma_tt, sigma(3, 3)];
      
   mu = syst(mu(1), mu(2), mu(3), g, T, rho_0, k_rho);
   G = syst_jacobian(mu_posterior (1), mu_posterior(2), mu_posterior(3), T, rho_0, k_rho);
   sigma = G*sigma*G' + R;
   H = [2*mu(1), 2*mu(2), 0];
   K = sigma*H'/(H*sigma*H' + Q);
    
   z_bar = (mu(1)^2+mu(2)^2);
   mu = mu + K*(ztrue(:,i+1) - z_bar);
   sigma = (eye(3) - K*H)*sigma;
   NEES = [NEES, ((xtrue(:,i+1) - mu)'/sigma)*(xtrue(:,i+1) - mu)];
   mu_posterior = mu;
   Xpred = [Xpred, mu];
   
end
t = [t, t(end)+T];
end