function M_t = compute_Bt(P_ks, P_ks1, M, nq, ny, N)
    M_t = zeros(3*nq+ny, 3*nq+ny, N);
    M_t(1:nq,1:2*nq,:) = [P_ks, M];
    M_t(1:nq,2*nq+ny+1:end,:) = P_ks;
    M_t(nq+1:2*nq,1:2*nq,:) = [M, P_ks1];
    M_t(nq+1:2*nq,2*nq+ny+1:end,:) = M;
    M_t(2*nq+ny+1:end, 1:2*nq,:) = [P_ks, M];
    M_t(2*nq+ny+1:end, 2*nq+ny+1:end,:) = P_ks;
end

