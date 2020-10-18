function [loc_est, G_mcmv, id_est, data_cov_inv] = beamformer_localizer_n(Y, noise_cov, G, mesh, included_sp, num_sources)
%addpath('/m/nbe/scratch/braintrack/pnas_mne_results')
num_trials = size(Y,3);
nchan=size(Y,1);
T = size(Y,2);
data_cov=0;
for j=1:1:num_trials
    data_cov = data_cov + (squeeze(Y(:,:,j))*squeeze(Y(:,:,j))');
end
data_cov = data_cov./T;
data_cov=double(data_cov);
%alpha=.1*trace(data_cov)/length(data_cov);
%data_cov=data_cov+alpha*eye();

for i=1:1:size(G,2)
    eta_ai(i) =  compute_ai(G(:,i), noise_cov , data_cov);
end

loc(1) = find(eta_ai == max(eta_ai));

for j=1:1:num_sources

    for i=1:1:size(G,2)
        id = find(loc == i);
        if isempty(id)
            eta_mai(i) = compute_mai([G(:,loc) G(:,i)], noise_cov , data_cov);
        end
    end

eta_mai(loc) = 0;
loc(j) = find(eta_mai == max(eta_mai));
end

loc_est = reshape(mesh.p(included_sp(loc),:)', 3*num_sources,1);
%id_est = sort(loc, 'ascend');
%G_mcmv = G(:, id_est);
G_mcmv = G(:,loc);
id_est = loc;
data_cov_inv = inv(data_cov);
