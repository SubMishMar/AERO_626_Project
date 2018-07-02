function mu_bar = syst(xk, yk, thetak, vk, wk, Ts)

mu_bar = [xk + Ts*cos(thetak)*vk;
          yk + Ts*sin(thetak)*vk;
          thetak + Ts*wk];
end