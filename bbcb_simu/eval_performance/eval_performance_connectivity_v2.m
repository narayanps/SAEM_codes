function eval_performance_connectivity_v2 (type, arr, method, sim_res_path, save_path)
path =sim_res_path;
addpath(path);
addpath(genpath(fullfile(pwd, '..','..', 'misc')));
addpath(genpath(fullfile(pwd, '..', 'bbcb_simu')));
addpath(genpath(fullfile(pwd, '..', 'bbcb_simu/data')));

%load head model
load miscdata
load sa
mesh.p= sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
G = sa.cortex75K.MEG_V_bem_normal(:, sa.cortex2K.in_from_cortex75K);
elind = unique(round(1:2.5:298));
elind = elind(1:108);
G=G(elind,:);
bio_snr_vals = [1 3 5 10];
bio_snr_arr = repmat(bio_snr_vals , [50 1]);
id_orig=find(bio_snr_arr==bio_snr_vals(arr));
fs = 100;
bs_ss=50;
count=0;
B=100;
for jj=1:1:B
    tic
    TP=0;
    FN=0;
    FP=0;
    TN=0;
    
    id = datasample(id_orig, bs_ss);
    
    for i=1:length(id)
        load(sprintf('truth_snr_%.1f_type_%d_sim_%d.mat',bio_snr_arr(id(i)), type,id(i)));
        Y=truth.Y;
        truth.order=2;
        if strcmp(method, 'eloreta')
            P = mkfilt_eloreta(G,0.01);
            sources_el = Y'*P;
            est{i,1}.amp(1,:)=sources_el(:,truth.r_id(1));
            est{i,1}.amp(2,:)=sources_el(:,truth.r_id(2));
            [est{i,1}.A, est{i,1}.V]=idMVAR(est{i,1}.amp,truth.order,2);
            est{i,1}.r = truth.r(:,1);
            
        elseif strcmp(method, 'lcmv')
            data=Y';
            segleng=fs;                 % length of each segment to estimate the FFT
            epleng=2*fs;                % length of one epoch (can be ignored here: only important for event-related analysis)
            segshift=segleng/2;         % shift of each segment. segleng/2 leads to 50% overlap of segments
            maxfreqbin=(segleng/2)+1;   % maximum frequency bin to estimate.
            
            para=[];
            para.segave = 1;            % average over segments
            para.subave = 0;            % subtract average
            
            [cs, ~, nave] = data2cs_event(data, segleng, segshift, epleng, maxfreqbin, para);
            
            rcs1f = mean(real(cs),3);
            [~, A1, po] = mkfilt_lcmv(G,rcs1f);
            
            mdats = data*A1;
            est{i,1}.amp(1,:)=mdats(:,truth.r_id(1));
            est{i,1}.amp(2,:)=mdats(:,truth.r_id(2));
            [est{i,1}.A, est{i,1}.V]=idMVAR(est{i,1}.amp,truth.order,2);
            est{i,1}.r = truth.r(:,1);
            
        elseif strcmp(method, 'mcmv')
            data_cov_inv = (Y*Y')./size(Y,2);
            wT = inv(G(:, truth.r_id)' * data_cov_inv * G(:, truth.r_id)) * G(:, truth.r_id)' * data_cov_inv;
            est{i,1}.amp = wT*truth.Y;
            [est{i,1}.A, est{i,1}.V]=idMVAR(est{i,1}.amp,2,2);
            est{i,1}.r = truth.r(:,1);
            
        elseif strcmp(method, 'saem')
            load(sprintf('params_snr_%.1f_type_%d_sim_%d.mat',bio_snr_arr(id(i)), type,id(i)))
            load(sprintf('state_snr_%.1f_type_%d_sim_%d.mat',bio_snr_arr(id(i)), type, id(i)))
            r_ = squeeze(state.r_hist(:,:,end));
            qs_ = squeeze(state.qs_hist(:,:,end));
            est{i,1}.r = mean(r_,2);
            est{i,1}.amp = qs_; %qs_(:,1:end-1);
            est{i,1}.A = squeeze(params.A_hist(:,:,end));
            est{i,1}.V = squeeze(params.V_hist(:,:,end));
            est{i,1}.sigma_m = squeeze(params.sigma_m_hist(1,end));
            est{i,1}.sigma_b = squeeze(params.sigma_b_hist(1,end));
            nle(i) = abs(truth.sigma_m - est{i,1}.sigma_m) / truth.sigma_m;
        end
        
        [TP, FP, TN, FN, cle(i,:), sle(i,:), ~, ~ ] = compute_error_stats_no_trial(est{i,1},truth, TP, FP, TN, FN);
    end
    res.FPR(jj,1) = FP / (FP+ TN);
    res.TNR(jj,1) = TN / (FP+ TN);
    if type==1
        res.TPR(jj,1) = TP / (TP + FN);
        res.FNR(jj,1) = FN / (TP + FN);
    end
    res.noc(jj,1) = (FP+TN)./length(id);
    id_valid=find(cle~=inf);
    res.cle(jj,1) = mean(cle(id_valid));
    if strcmp(method, 'saem')
        
        res.sle(jj,:) = mean(mean(sle));

        res.nle(jj,1) = mean(nle(id_valid));
    end
    toc
end
fname = sprintf(strcat(save_path, '/res_%s_%d_%d'), method, bio_snr_vals(arr), type);
save(fname, 'res');
end
    