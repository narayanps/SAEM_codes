function M = one_lag_gpu(A, P_kf, J_tm1, M_tp1, J_t)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
P_kf =  gpuArray(P_kf);
J_tm1=gpuArray(J_tm1);
M_tp1=gpuArray(M_tp1);
J_t=gpuArray(J_t);
A=gpuArray(A);
M = gputimes(P_kf, gputranspose(J_tm1)) + gputimes(gputimes(J_t, (M_tp1 - gputimes(A, P_kf))), gputranspose(J_tm1));      
M=gather(M);
end

