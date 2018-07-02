function Gt = syst_jacobian(xk, yk, thetak, vk, wk, Ts)
Gt = [1, 0, -Ts*sin(thetak)*vk;
      0, 1,  Ts*cos(thetak)*vk;
      0, 0, 1];
end