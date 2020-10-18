function [aN] = compute_aN_fast(r_pf_t, r_pf, q_kf, P_kf, V_r, ...
                                omega_t, lambda_t, num_sources, w, N, CONST, log_det)
lam = zeros(num_sources,num_sources,N);
eta = zeros(N,1);
logw_rb = zeros(N,1);
for i=1:1:N
    gamma = chol(P_kf(1:num_sources,1:num_sources,i));
    lam(:,:,i) = gamma'*omega_t*gamma + eye(num_sources);
    tmp = gamma'*(lambda_t - omega_t*q_kf(1:num_sources,i));
    eta(i) = (q_kf(1:num_sources,i)'*omega_t)*q_kf(1:num_sources,i)...
        - 2*lambda_t'*q_kf(1:num_sources,i) - ...
        tmp'/squeeze(lam(:,:,i))*tmp;
    log_r = CONST - 0.5*log_det ...
        -0.5*((squeeze(r_pf_t(:,end))-squeeze(r_pf(:,i)))'/ squeeze(V_r(:,:)) * (squeeze(r_pf_t(:,end))-squeeze(r_pf(:,i))));
    logw_rb(i) = log(w(i)) -eta(i)/2 - log(det(lam(:,:,i)))/2 + log_r;
end
logw_rb = logw_rb - max(logw_rb);
waN =  exp(logw_rb);
waN = waN./sum(waN);
aN = find(mnrnd(1,waN'));
end