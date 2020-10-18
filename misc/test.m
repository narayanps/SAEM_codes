addpath('/m/nbe/scratch/braintrack/pnas_mne_results')
addpath('/m/nbe/scratch/braintrack/pnas_codes/MNE')
load G
load mesh_ds;
load included_sp
mesh_ds.p = mesh_ds.p * 1000;
for i=1:1:100
    load (sprintf('truth_snr_0.5_type_1_sim_%d.mat', i+0))
    [loc_error(i,:), success_rate(i)] = beamformer_localizer(truth, G, mesh_ds, included_sp);
end


for i=59
    success_rate_saem(i)  = 0;
    load (sprintf('truth_snr_10.0_type_1_sim_%d.mat', i+400))
    load (sprintf('state_snr_10.0_type_1_sim_%d.mat', i+400))
    load (sprintf('params_snr_10.0_type_1_sim_%d.mat', i+400))
    nq=2;
    loc_true = mesh_ds.p(truth.r_id,:);
    loc_est = mean(squeeze(state.r_hist(:,:,end)),2);
    [comb_id, min_mean_err, err(i,:)] = check_loc_error(reshape(loc_true', 3*nq,1), loc_est, 2); 
    if max(err(i,1)) <= 10
        success_rate_saem(i)=1;
    end
end