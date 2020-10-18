function [q, P] = kalman_predict_gpu(q, P, A, V_q )
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'


q=gpuArray(q);
P=gpuArray(P);
A= gpuArray(A);
V_q=gpuArray(V_q);
q = gputimes(A,q);
P = gputimes(gputimes(A,P),gputranspose(A)) + V_q;
q=gather(q);
P=gather(P);
end