function [M] = init_M_gpu(I, Kf, G, P_kf, A)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gpurdivide = @(A, B) pagefun(@mrdivide, A, B);  % A/B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
I=gpuArray(I);
Kf = gpuArray(Kf);
G=gpuArray(G);
P_kf = gpuArray(P_kf);
A=gpuArray(A);
M =  gputimes(gputimes((I - gputimes(Kf, G)), A), P_kf);
M=gather(M);
end

