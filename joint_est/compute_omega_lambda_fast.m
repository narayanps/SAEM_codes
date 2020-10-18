function [omega_t, lambda_t] = compute_omega_lambda_fast(A, V_q, G_pf, Y, XI, num_sources, p, N, T)
nq = num_sources*p;
F_ = chol(V_q(1:num_sources, 1:num_sources));
% F_ = zeros(nq,nq);
% F_(1:num_sources, 1:num_sources) = F;
% F_(num_sources+1:num_sources*p, 1:num_sources)=0;
% F_(1:num_sources*p, num_sources+1:end)=0;

omega_t = zeros(num_sources,num_sources,T-1);
lambda_t = zeros(num_sources,T-1);
tmp_L = squeeze(G_pf(:,:,T,N)); % repmat(zeros(size(G_pf(:,:,T,N))), 1, p-1)];
lambda_hat = tmp_L'/XI * Y(:,T);
omega_hat =  tmp_L' /XI * tmp_L;


for t = T-1:-1:1
    M_t1  = F_' * omega_hat * F_ + eye(num_sources);
    m_t1 = lambda_hat;
    omega_t(:,:,t)= A' * (eye(num_sources) - omega_hat * F_ /M_t1*F_')*omega_hat*A;
    lambda_t(:,t) = (A'*(eye(num_sources)-omega_hat*F_/M_t1*F_') * m_t1);
    tmp_L = squeeze(G_pf(:,:,t,N)) ;%repmat(zeros(size(G_pf(:,:,t,N))), 1, p-1)];
    omega_hat = omega_t(:,:,t) + tmp_L' /XI * tmp_L;
    lambda_hat = lambda_t(:,t) + tmp_L'/XI * Y(:,t);
    
end
end