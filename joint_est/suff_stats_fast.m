function [S] = suff_stats_fast(S, q_ks, q_ks1, Y_, M_t_s,  N, nq, ny, w_1)
x_t_hat = zeros(3*nq+ny, N);
x_t_hat(1:nq,:)  = q_ks;
x_t_hat(nq+1:2*nq,:) = q_ks1;
x_t_hat(2*nq+1:ny+2*nq,:) = Y_;
x_t_hat(2*nq+1+ny:end,:) = q_ks;
x_t_hat_w =  w_1.* x_t_hat;
S = S + (x_t_hat_w * x_t_hat' + M_t_s);
end
    