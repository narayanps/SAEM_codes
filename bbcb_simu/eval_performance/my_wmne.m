function [j] = my_wmne(Y, E, G, rank_E, SNR, chan_ids, pca_flag, method)
num_chans = size(E,2);
num_sources = size(G,2);
lambda2 = 1/(SNR);
W = zeros(0, num_chans);
if ~isempty(chan_ids)
    W_meg = prep_whitener(E, chan_ids, rank_E, pca_flag);
    W_tmp = zeros(size(W_meg,1), num_chans);
    W_tmp(:,chan_ids) = W_meg;
    W = [W; W_tmp];
end

if any(isnan(W(:))) || any(isinf(W(:)))
    error('Invalid noise covariance matrix.')
end

% whiten Lead field matrix
G = W * G;

% scale source-covariance matrix
w = ones(num_sources,1);
Vj = speye(num_sources, num_sources);
Vj = spdiags(w, 0, Vj);
trclcl = trace(G * Vj * G');
Vj = Vj * (rank_E / trclcl);
Rc = chol(Vj, 'lower');
GW = G * Rc;

%svd
[U,S,V] = svd(GW,'econ');
s = diag(S);
ss = s ./ (s.^2 + lambda2);
Winv = Rc * V * diag(ss) * U';

%unweighted mne
if strcmp(method, 'mne')
Winv_mne = Winv * W;
j = Winv_mne * Y;

elseif strcmp(method, 'dspm')
%dspm 
Winv_dspm=zeros(size(Winv));
for i=1:1:num_sources
    dspmdiag = sum(Winv(i,:) .^2, 2);
    dspmdiag = sqrt(dspmdiag);
    Winv_dspm(i,:) = bsxfun(@rdivide, Winv(i,:), dspmdiag);
end
Winv_dspm = Winv_dspm * W;
j = Winv_dspm * Y;

elseif strcmp(method, 'sloreta')
%sloreta
Winv_sloreta=zeros(size(Winv));
for i=1:1:num_sources
sloretadiag = sqrt(sum(Winv(i,:) .* G(:,i)', 2));
Winv_sloreta(i,:) = bsxfun(@rdivide, Winv(i,:), sloretadiag);
end
Winv_sloreta = Winv_sloreta * W;
j = Winv_sloreta * Y;
end