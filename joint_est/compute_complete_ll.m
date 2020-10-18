function [ll] = compute_complete_ll(T, stats, V_q, A, E, G_s)
logdetV = 2*sum(log(diag(chol(V_q))));
logdetE = 2*sum(log(diag(chol(E))));
ll= -1/2*(T*logdetV + trace(V_q\(stats.Phi - stats.Psi*A' - A*stats.Psi' + A*stats.Sigma*A'))...
        + T*logdetE + trace(E\(stats.Omega - stats.Lambda*G_s' - G_s*stats.Lambda' + G_s*stats.Xi*G_s')));
end

