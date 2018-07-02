function mu_bar = syst(xk, yk, zk, g, T, rho_0, k_rho)

mu_bar = [xk + T*yk;
          yk + T*rho_0*exp(-xk/k_rho)*yk^2/(2*zk) - T*g;
          zk];
end