function [logp] = logpdf(x, mu, Sigma)
%Author : Narayan Subramaniyam
%Aalto/NBE
m = size(x,1);

CONST = -m * 0.5 * log(2*pi);
logDetSigma = 2*sum(log(diag(chol(Sigma))));
%logp = CONST - 0.5*logDetSigma - 0.5*((x-mu)' * ...
%choleskeyinverse(chol(Sigma)) * (x-mu));
logp = CONST - 0.5*logDetSigma - 0.5*((x-mu)'/Sigma * (x-mu));
%logp =  - 0.5*((x-mu)'/Sigma * (x-mu));