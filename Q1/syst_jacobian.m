function Gt = syst_jacobian(xk, yk, zk, T, rho_0, k_rho)
Gt = [1, T, 0;
      -T*(yk^2)*rho_0*exp(-xk/k_rho)/(2*k_rho*zk), 1 + T*rho_0*exp(-xk/k_rho)*yk/zk, -T*rho_0*exp(-xk/k_rho)*yk^2/(2*zk^2);
      0, 0, 1];
end