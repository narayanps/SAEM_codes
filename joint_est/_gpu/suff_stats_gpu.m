
function [S] = suff_stats_gpu(q_ks, q_ks1, Y_,  N, nq, ny, w_1, J, T)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
gputimes_e = @(A, B) pagefun(@times, A, B);      % A*B on GPU

 q_ks = permute(q_ks, [1 3 4 2]);
 q_ks1 = permute(q_ks1, [1 3 4 2]);
 Y_ = permute(Y_, [1 4 3 2]);
 
 q_ks=reshape(squeeze(q_ks(:,:,:)), nq, N, J*T);
 q_ks1=reshape(squeeze(q_ks1(:,:,:)), nq, N, J*T);
 Y_=reshape(squeeze(Y_), ny, N, J*T);

x_t_hat = single(zeros(3*nq+ny, N,J*T));
x_t_hat(1:nq,:,:)  = q_ks;
x_t_hat(nq+1:2*nq,:,:) = q_ks1;
x_t_hat(2*nq+1:ny+2*nq,:,:) = Y_;
x_t_hat(2*nq+1+ny:end,:,:) = q_ks;
clear q_ks q_ks1 Y_
x_t_hat = gpuArray(x_t_hat);
x_t_hat_w =  gputimes_e(w_1, x_t_hat);
x_t_hat_w=gpuArray(x_t_hat_w);
%M_t_s = gpuArray(repmat(M_t_s, [1 1 1 J]));
%M_t_s = reshape(M_t_s, 3*nq+ny,3*nq+ny,J*T);
S =  (gputimes(x_t_hat_w , gputranspose(x_t_hat)));% + M_t_s);
clear x_t_hat x_t_hat_w
S=gather(S);
end
