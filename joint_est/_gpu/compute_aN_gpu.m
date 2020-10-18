function [aN] = compute_aN_gpu(r_pf_t, r_pf, q_kf, P_kf, V_r, ...
omega_t, lambda_t, w, N, CONST, log_det)
np=size(q_kf,1);

gamma = zeros(np, np,N);
logw_rb = zeros(N,1);
I=eye(np);
for i=1:1:N
    gamma(:,:,i) = chol(P_kf(1:np,1:np,i));
end
q_kf_r (:,1,:) = q_kf;
q_kf_gpu = gpuArray(q_kf_r);
gamma_gpu = gpuArray(gamma);
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gpurdivide = @(A, B) pagefun(@mrdivide, A, B);  % A/B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
lam = gputimes(gputimes(gputranspose(gamma_gpu),omega_t), gamma_gpu) + I;
tmp = (gputimes(gputranspose(gamma_gpu), lambda_t - gputimes(omega_t, q_kf_gpu)));
eta_ = gputimes(gputimes(gputranspose(q_kf_gpu), omega_t), q_kf_gpu) - gputimes(2*lambda_t', q_kf_gpu) - gputimes(gpurdivide(gputranspose(tmp),lam),tmp);
eta = gather(eta_);
lam=gather(lam);
for i=1:1:N
    log_r = CONST - 0.5*log_det ...
        -0.5*((squeeze(r_pf_t(:,end))-squeeze(r_pf(:,i)))'/ squeeze(V_r(:,:)) * (squeeze(r_pf_t(:,end))-squeeze(r_pf(:,i))));
    logw_rb(i) = log(w(i)) -eta(i)/2 - log(det(lam(:,:,i)))/2 + log_r;
end
logw_rb = logw_rb - max(logw_rb);
waN =  exp(logw_rb);
waN = waN./sum(waN);
aN = find(mnrnd(1,waN'));
