function plot_file(xtrue, Xpred, sigma_xx, sigma_yy, sigma_tt, time, plotname)
figure('Name',plotname, 'NumberTitle','off');
subplot(221)
p = plot(time, xtrue(1,:),'--r', 'LineWidth', 3);
hold on;
q = plot(time, Xpred(1,:), 'g', 'LineWidth', 2);
hold on;
r = plot(time, xtrue(1,:) - Xpred(1,:), 'b', 'LineWidth', 2);
hold on;
s = plot(time, Xpred(1,:) + 3*sqrt(sigma_xx), '--k', 'LineWidth', 2);
hold on;
t = plot(time, Xpred(1,:) - 3*sqrt(sigma_xx), '--k', 'LineWidth', 2);
hold off;
title('True and Predicted x1[ft] vs t[s]','Interpreter','latex');
grid;
u = legend([p q r s ],'$x1_{true}$','$x1_{pred}$', '$e_{x1}$', '$\pm3\sigma$');
set(u,'Interpreter','latex');
% xlabel('t[s]','Interpreter','latex')
% ylabel('$x_{true}$ $\&$ $x_{pred}$ [m]','Interpreter','latex')


subplot(222)
p = plot(time, xtrue(2,:),'--r', 'LineWidth', 3);
hold on;
q = plot(time, Xpred(2,:),'g', 'LineWidth', 2);
hold on;
r = plot(time, xtrue(2,:) - Xpred(2,:), 'b', 'LineWidth', 2);
hold on;
s = plot(time, Xpred(2,:) + 3*sqrt(sigma_yy), '--k', 'LineWidth', 2);
hold on;
t = plot(time, Xpred(2,:) - 3*sqrt(sigma_yy), '--k', 'LineWidth', 2);
hold off;
title('True and Predicted x2[ft/s] vs t[s]','Interpreter','latex');
grid;
u = legend([p q r s ],'$x2_{true}$','$x2_{pred}$', '$e_{x2}$', '$\pm3\sigma$');
set(u,'Interpreter','latex');
% xlabel('t[s]','Interpreter','latex')
% ylabel('$y_{true}$ $\&$ $y_{pred}$ [m]','Interpreter','latex')

subplot(223)
p = plot(time, xtrue(3,:),'--r', 'LineWidth', 3);
hold on;
q = plot(time, Xpred(3,:),'g', 'LineWidth', 2);
hold on;
r = plot(time, xtrue(3,:) - Xpred(3,:), 'b', 'LineWidth', 2);
hold on;
s = plot(time, Xpred(3,:) + 3*sqrt(sigma_tt), '--k', 'LineWidth', 2);
hold on;
t = plot(time, Xpred(3,:) - 3*sqrt(sigma_tt), '--k', 'LineWidth', 2);
hold off;
title('True and Predicted x3[lb/ft2] vs t[s]','Interpreter','latex');
grid;
u = legend([p q r s ],'$x3_{true}$','$x3_{pred}$', '$e_{x3}$', '$\pm3\sigma$');
set(u,'Interpreter','latex');
% xlabel('t[s]','Interpreter','latex')
% ylabel('$\theta_{true}$ $\&$ $\theta_{pred}$ [rad]','Interpreter','latex')

subplot(224)

p = plot(time, sigma_xx, 'r', 'LineWidth', 2);
%%title('$e_{x}$ [m] vs t[m]','Interpreter','latex');
%%grid;
hold on
q = plot(time, sigma_yy, 'g', 'LineWidth', 2);
%%title('$e_{y}$ [m] vs t[m]','Interpreter','latex');
%%grid;
hold on
r = plot(time, sigma_tt, 'b', 'LineWidth', 2);
%%title('$e_{\theta}$ [m] vs t[m]','Interpreter','latex');
title(' $\Sigma_{x1x1}$, $\Sigma_{x2x2}$, $\Sigma_{x3x3}$ vs t[s]','Interpreter','latex');
grid;
hold off;
s = legend([p q r],'$\Sigma_{x1x1}$','$\Sigma_{x2x2}$', '$\Sigma_{x3x3}$');
set(s,'Interpreter','latex');
end