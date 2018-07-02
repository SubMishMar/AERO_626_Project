function [Xpred, sigma_xx, sigma_yy, sigma_tt, NEES, time, x_P] = q2_PF(mu, sigma, xtrue, ztrue, Ts, Q, R, Tend, N)

sigma_xx = [sigma(1,1)];
sigma_yy = [sigma(2,2)];
sigma_tt = [sigma(3,3)];
NEES = [];
%initilize our initial, prior particle distribution as a gaussian around
%the true initial value

x_P = [];
n = 3;
m = 2;
% make the randomly generated particles from the initial prior gaussian distribution
for i = 1:N
    x_P(:,i) = mu + sqrt(sigma) * randn(n,1);
end

Xpred = [mu]; % the vector of particle filter estimates.


%%
count = 1;
for j = 0:Ts:Tend
    sigma_xx = [sigma_xx sigma(1,1)];
    sigma_yy = [sigma_yy sigma(2,2)];
    sigma_tt = [sigma_tt sigma(3,3)];
    vk = sin(j);
    if j <= 50
       wk = 0.1;
    elseif j > 50 && j <= 80
       wk = -0.2;
    else
       wk = -0.1;
    end
 
    if mod(j, 5*Ts) == 0 && j ~= 0
        %disp('measurement taken in PF');
     for i = 1:N
        x_P_update(:,i) = syst(x_P(1,i), x_P(2,i), x_P(3,i), vk, wk, Ts) + sqrt(R)*randn(n,1);
        z_update(:,i) = measurement_model(x_P_update(:,i));
        P_w(:,i) = (1/sqrt(det(2*pi*Q))) * exp(-0.5*(ztrue(:,count+1) - z_update(:,i))'*inv(Q)*(ztrue(:,count+1) - z_update(:,i)));
     end
    % Normalize to form a probability distribution (i.e. sum to 1).
     P_w = P_w./sum(P_w);
    % Resampling: From this new distribution, now we randomly sample from it to generate our new estimate particles
     for i = 1 : N
        x_P(:,i) = x_P_update(:, find(rand <= cumsum(P_w),1));
     end
     x_est = mean(x_P,2);
     Pa = cov(x_P');
     sigma = Pa;
    else
     for i = 1:N
        x_P_update(:,i) = syst(x_P(1,i), x_P(2,i), x_P(3,i), vk, wk, Ts) + sqrt(R)*randn(n,1);
     end
     x_est = mean(x_P_update,2);
     Pa = cov(x_P_update');
     sigma = Pa; 
     x_P = x_P_update;
    end
    NEES = [NEES, ((xtrue(:,count+1) - x_est)'/sigma)*(xtrue(:,count+1) - x_est)];
    Xpred = [Xpred x_est];
    count = count + 1;
end
time = [0:Ts:Tend, Tend+Ts];
end
