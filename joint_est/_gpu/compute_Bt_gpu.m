function M_t = compute_Bt_gpu(P_ks, P_ks1, M, nq, ny, N,T)
gputimes = @(A, B) pagefun(@mtimes, A, B);   

    M_t = zeros(3*nq+ny, 3*nq+ny, N);
    M_t=gpuArray(M_t);
    M=gpuArray(M);
    P_ks=gpuArray(P_ks);
    P_ks1=gpuArray(P_ks1);
    M_t(1:nq,1:2*nq,:) = gputimes(1, gpuArray([P_ks, M]));
    M_t(1:nq,2*nq+ny+1:end,:) = gputimes(1, P_ks);
    M_t(nq+1:2*nq,1:2*nq,:) = gputimes(1, gpuArray([M, P_ks1]));
    M_t(nq+1:2*nq,2*nq+ny+1:end,:) = gputimes(1, M);
    M_t(2*nq+ny+1:end, 1:2*nq,:) = gputimes(1, gpuArray([P_ks, M]));
    M_t(2*nq+ny+1:end, 2*nq+ny+1:end,:) = gputimes(1, P_ks);
  
 M_t=gather(M_t);
    

end