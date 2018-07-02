function p = plot_ellipse(Xpred, sigma_xx, sigma_yy, on, color)
a = 3*sqrt(sigma_xx(end));
b = 3*sqrt(sigma_yy(end));
h = Xpred(1, end);
k = Xpred(2, end);
theta = 0:0.01:2*pi;
X = [h*ones(size(theta)) + a*cos(theta);
     k*ones(size(theta)) + b*sin(theta)];

p = plot(X(1,:),X(2,:), color, 'LineWidth', 3);
hold on;
Line = [h, h+0.1*cos(Xpred(3,end));
        k, k+0.1*sin(Xpred(3,end))];

p1 = Line(:,1);
p2 = Line(:,2);
dp = p2 - p1;
quiver(p1(1),p1(2),dp(1),dp(2),0, color , 'LineWidth', 3);
hold on;
plot(h, k,'-mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',color,...
    'MarkerSize',10);
if on==1
   hold on;
else
   hold off;
end

end