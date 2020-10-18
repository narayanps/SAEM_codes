function [err, success_rate] = beamformer_localizer(truth, G, mesh, included_sp)
addpath('/m/nbe/scratch/braintrack/pnas_mne_results')
G=G(1:3:306,:);
data_cov = (truth.Y*truth.Y')./5000;
n_noise_sources=size(G,2);
bn = randn(n_noise_sources, 5000);
Y_brain_noise = G(:, :)*bn;
[Y_baseline, ~, XI_baseline, ~ ] = add_meas_noise(Y_brain_noise,3);
noise_cov = (Y_baseline*Y_baseline')./5000;

success_rate=0;

for i=1:1:7498
    eta_ai(i) =  compute_ai(G(:,i), noise_cov , data_cov);
end

first_loc = find(eta_ai == max(eta_ai));

for i=1:1:7498
    if i ~= first_loc
        eta_mai(i) = compute_mai([G(:,first_loc) G(:,i)], noise_cov , data_cov);
    end
end

second_loc = find(eta_mai == max(eta_mai));

loc_est = mesh.p(included_sp([first_loc; second_loc]),:);
loc_true = mesh.p(truth.r_id,:);
nq=2;
[comb_id, min_mean_err, err] = check_loc_error(reshape(loc_true', 3*nq,1), reshape(loc_est', 3*nq, 1), 2); 
if max(err) <= 10
    success_rate=1;
end
