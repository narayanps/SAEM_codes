function [log_weights] = gpu_logpdf(Y,ypred,sigma, logDetSigma,N,ny)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gpurdivide = @(A, B) pagefun(@mrdivide, A, B);  % A/B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
y1(:,1,:) = ypred - repmat(Y, [1 N]);
a1=gpuArray(y1);
b = gpuArray(sigma);
logDetSigma_g(1,1,:) = logDetSigma;
logDetSigma_g = gpuArray(logDetSigma_g);
log_weights_1 = -ny * 0.5 * log(2*pi) - 0.5*logDetSigma_g - gputimes(0.5, gputimes(gpurdivide(gputranspose(a1),b),a1));
log_weights_1=gather(log_weights_1);
log_weights=squeeze(log_weights_1);
end

