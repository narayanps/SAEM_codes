function [TPR, FPR, no_of_correct_hits, reconst_error] = eval_performance_meas_noise (type, meas_snr, method)
path ='/scratch/work/puthann1/codes_saem_paper/results/bbcb_meas';
addpath(path);
path ='/scratch/work/puthann1/codes_saem_paper/misc';
addpath(path);
path = '/scratch/work/puthann1/codes_saem_paper/mvar_functions';
addpath(genpath(path));
load('/scratch/work/puthann1/codes_saem_paper/bbcb_simu/data/sa')
load('/scratch/work/puthann1/codes_saem_paper/bbcb_simu/data/miscdata')

mesh.p= sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :); % %sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
mesh.e=sa.cortex2K.tri;
mesh.nn = sa.cortex75K.normals(sa.cortex2K.in_from_cortex75K, :) ; %; %sa.cortex75K.normals(sa.cortex2K.in_from_cortex75K, :);
L=sa.cortex75K.MEG_V_bem(:, sa.cortex2K.in_from_cortex75K, :);
elind = unique(round(1:2.5:298));
elind = elind(1:108);
L = L(elind, :, :);

Lred = [];
for ivox = 1:size(L, 2)
    [U_, ~, ~] = svd(squeeze(L(:, ivox, :)), 'econ');
    Lred(:, ivox, :) = U_(:, 1:2);
end
L = Lred;
meas_snr_vals = [1 3 5 10];
meas_snr_arr = repmat(meas_snr_vals , [100 1]);
id_orig=find(meas_snr_arr==meas_snr);
bs_ss=100;
for bs=1:1:30
    tic
    TP=0;
    FN=0;
    FP=0;
    TN=0;
    id = datasample(id_orig, bs_ss);
for i=1:100
    load(sprintf('truth_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type,id(i)));
    Y=truth.Y;
    G=truth.model.G;
    E=truth.model.XI;
    truth.order=5;
    
    
    if strcmp(method, 'mne') ||strcmp( method, 'sloreta') || strcmp(method, 'dspm')
        chan_ids = 1:size(Y,1);
        pca_flag=0;
        [soln_inv] = my_wmne(Y, E, G, size(Y,1), meas_snr, chan_ids, pca_flag, method);
        est{i,1}.amp = soln_inv(truth.r_id,:);
        [est{i,1}.A, est{i,1}.V]=idMVAR(est{i,1}.amp,5,2);
        est{i,1}.r = truth.r(:,1);

    elseif strcmp(method, 'lcmv')
        cov=zeros(size(Y,1));
        Ntotal=zeros(size(Y,1));
        T=size(Y,2);
        Y_ = Y;
        Y_=Y_-repmat(mean(Y_,2),[1 T]);
        cov = Y_*Y_';
        Ntotal = T;
        cov = cov ./ (Ntotal-1);
        [~, A1, ~]=mkfilt_lcmv(L,cov);
        est{i,1}.amp= ((Y(:,:))'*A1(:,[truth.r_id]))';
        [est{i,1}.A, est{i,1}.V]=idMVAR(est{i,1}.amp,5,2);
        est{i,1}.r = truth.r(:,1);

       
    else
        load(sprintf('loc_est_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type,id(i)))
        load(sprintf('amp_est_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type,id(i)))
        load(sprintf('A_est_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type,id(i)))
        load(sprintf('V_q_est_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type,id(i)))
        load(sprintf('XI_est_snr_%d_type_%d_sim_%d.mat',meas_snr_arr(id(i)), type, id(i)))
        est{i,1}.r = mean(loc_est,2);
        est{i,1}.A = squeeze(A(1:2,:,end));
        est{i,1}.XI = squeeze(XI(:,:,end));
        est{i,1}.V = squeeze(V(1:2,1:2,end));
        est{i,1}.amp = amp_est(1:2,1:end-1);
        vle(i,:) = norm(truth.model.XI-est{i,1}.XI, 'fro')./norm(truth.model.XI, 'fro');

        
        % map smoothed locations on head model
        for j=1:1:2
            [~, est{i,1}.r_id(j)] = find_nearest_id(mesh.p,est{i,1}.r((j-1)*3+1 : 3*j,1)');
            est{i,1}.r((j-1)*3+1 : 3*j,1) = mesh.p(est{i,1}.r_id(j),:);
        end
        for j=1:1:8
            for k=1:1:2
                if find(inds_roi_inner_2K{1,j} == est{i,1}.r_id(k))
                    est{i,1}.rois(k) = j;
                elseif find(inds_roi_outer_2K{1,j} == est{i,1}.r_id(k))
                    est{i,1}.rois(k) = j;
                end
            end
        end
        
        
    end
    [TP, FP, TN, FN, cle(i,:), sle(i,:) ] = compute_error_stats_no_trial(est{i,1},truth, TP, FP, TN, FN);

    
end
id = find(cle ~= inf);
mean_cle = mean(cle(id));
if strcmp(method, 'saem')
    mean_sle= mean(mean(sle));
    mean_meas_var = mean(vle(id));
end

reconst_error.pdc(bs) = mean_cle; 
if strcmp(method, 'saem')
    reconst_error.dipole(bs) = mean_sle;
    reconst_error.meas_var(bs) = mean_meas_var;
end
TPR(bs) = TP / (TP + FN);
FPR(bs) = FP / (FP+ TN);
no_of_correct_hits(bs) = (FP+ TN);
toc
end

