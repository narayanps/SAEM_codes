function [J, q_ks, P_ks] = kalman_smoother_gpu(P_kf_t, P_kp_t1, q_kf_t,q_kp_t1, q_ks_t1, P_ks_t1, A, N, Jt)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gpurdivide = @(A, B) pagefun(@mrdivide, A, B);  % A/B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
P_kf_t=gpuArray(P_kf_t);
P_kp_t1=gpuArray(P_kp_t1);
q_kf_t=gpuArray(q_kf_t);
q_ks_t1=gpuArray(q_ks_t1);
P_ks_t1=gpuArray(P_ks_t1);
q_kp_t1 = gpuArray(q_kp_t1);
A=gpuArray(A);
nq=size(q_kp_t1,1);
J = gpurdivide(gputimes(P_kf_t, gputranspose(A)),P_kp_t1);
q(:,1,:) = (q_ks_t1-q_kp_t1);
J_ = repmat(J, [1 1 1 Jt]);
J_=reshape(J_, nq,nq,N*Jt);
q_ks = q_kf_t + squeeze(gputimes(J_, q));
P_ks = P_kf_t + gputimes(gputimes(J, (P_ks_t1-P_kp_t1)), gputranspose(J));
q_ks = gather( squeeze(q_ks));
P_ks = gather(P_ks);
J=gather(J);
end

