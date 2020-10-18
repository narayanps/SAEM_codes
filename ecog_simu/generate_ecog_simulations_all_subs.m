function  generate_ecog_simulations_all_subs(sub, meg_path)
%meg_path is the MEG data path in the folder where data from openfMRI is
%saved
meas_snr=3;
subs = {'sub002', 'sub003', 'sub004','sub006','sub007','sub008','sub009'...
     'sub010', 'sub011','sub012','sub013','sub014','sub015','sub017'...
     'sub018', 'sub019'};
sub_id = [2:4 6:15 17:19];
vals = [1:length(sub_id)];
arr = repmat(vals , [20 1]); %20 mc iteration for each subject
addpath(sprintf(strcat(meg_path,'/%s', '/matlab'),  subs{arr(sub)}))
addpath(genpath('/models_for_estim'))
addpath(genpath(fullfile(pwd, '..', 'misc')));
addpath(genpath(fullfile(pwd, '..', 'joint_est')));

%load head model

load (strcat(meg_path,'/%s', '/matlab/G_%s'),  subs{arr(sub)},subs{arr(sub)});
load (strcat(meg_path,'/%s', '/matlab/included_sp_%s'),  subs{arr(sub)},subs{arr(sub)});
mm=load (strcat(meg_path,'/%s', '/matlab/mesh_%s'),  subs{arr(sub)},subs{arr(sub)});
load (strcat(meg_path,'/%s', '/matlab/anat_ds_%s'),  subs{arr(sub)},subs{arr(sub)});
load erp_ecog_seg

[status, seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(sub);
sd = rng;
truth.sd=sd;

%anat parcel no. corresponding to source locations
anat_no = [24 14]; %32];
num_sources = 2;
mesh_ds = mm.mesh{1,1};
mesh_ds.p=mesh_ds.p*1000; %in mm
%select sources that are atleast 20 mm apart
dmin=0;
distMat = inf*eye(num_sources,num_sources);
fg_pnts = mesh_ds.p(included_sp(anat_ds.parcels{14,1}),:);
loc_pnts = mesh_ds.p(included_sp(anat_ds.parcels{24,1}),:);


%point corresponding to SO1 elec
cent=mean(loc_pnts);
[~, id] = find_nearest_id(cent, loc_pnts(:,:));
LF_norm=sqrt(sum(G(1:3:306,anat_ds.parcels{24,1}).* G(1:3:306,anat_ds.parcels{24,1}),1));
LF_norm = LF_norm./max(LF_norm);
LF_norm(id)
vert_(1) = anat_ds.parcels{24,1}(id);
truth.r(1,:) = loc_pnts(id,:);
%point corresponding to IO3 elec
%'length of fg'
fg_len = max(fg_pnts(:,2)) - min(fg_pnts(:,2));
%shift approx 30% from the posterior
y_fg_shift = min(fg_pnts(:,2)) + 0.3*fg_len;
tmp_ids =  find(y_fg_shift-1 < fg_pnts(:,2) & y_fg_shift+1 > fg_pnts(:,2));
fg_id = tmp_ids(find(fg_pnts(tmp_ids,1) == max(fg_pnts(tmp_ids,1))));
truth.r(2,:) = fg_pnts(fg_id(1),:);
vert_(2) = anat_ds.parcels{14,1}(fg_id(1));
LF_norm=sqrt(sum(G(1:3:306,anat_ds.parcels{14,1}).* G(1:3:306,anat_ds.parcels{14,1}),1));
LF_norm = LF_norm./max(LF_norm);
LF_norm(fg_id)

truth.r=reshape(truth.r', 3*num_sources,1);
truth.r_id=vert_;
J=size(erp_ecog_seg,3);
G=G(1:3:306,:); %omly mag
num_chan=size(G,1);
bio_snr = 1;
num_steps = size(erp_ecog_seg,1);
%generate MEG
for j=1:1:J
    Y_brain_signal(:,:,j) = G(:, vert_) * 10^(-9)*squeeze(erp_ecog_seg(:,:,j))';
end

T=size(Y_brain_signal,2);

for j = 1:J
    %Calculate the observed signal power
    obser_sig_covar_matrix = abs(Y_brain_signal(:,:,j)*Y_brain_signal(:,:,j)')/(size(Y_brain_signal,2));
    obser_signal_power(j) = trace(obser_sig_covar_matrix);
end

% add colored noise
for j=1:1:J
    n_noise_sources=2000; 
    ind_noise_rand=randperm(size(G,2));
    ind_noise=ind_noise_rand(1:n_noise_sources)';   
    pn = mkpinknoise(num_steps, n_noise_sources)';
    Y_brain_noise(:,:,j) = G(:, ind_noise)*pn;
    brain_noise_covar_matrix=(Y_brain_noise(:,:,j)*Y_brain_noise(:,:,j)')/(size(Y_brain_noise,2));  
    snr_sf = obser_signal_power(j)/bio_snr;
    brain_noise_power=trace(brain_noise_covar_matrix);
    snr_sf = snr_sf/brain_noise_power;
    brain_noise_adj = sqrt(snr_sf)*Y_brain_noise(:,:,j);
    Y_(:,:,j) = Y_brain_signal(:,:,j) + brain_noise_adj;
    
end


%now add sensor level noise
for j = 1:J
    obser_sig_covar_matrix = abs(Y_(:,:,j)*Y_(:,:,j)')/(size(Y_,2));
    obser_signal_power(j) = trace(obser_sig_covar_matrix);
end
randn_noise = randn(num_chan,T,J);
for j = 1:J
    obser_noise = randn_noise(:,:,j);                                                 %white noise
    
    %Adjust power noise to get a proper SNR level
    noise_covar_matrix=(obser_noise*obser_noise')/(size(obser_noise,2));      %temporal covariance matrix
    
    snr_sf = obser_signal_power(j)/meas_snr;
    noise_power=trace(noise_covar_matrix);
    snr_sf = snr_sf/noise_power;
    obser_noise_adj = sqrt(snr_sf)*obser_noise;
    
    Y_seg(:,:,j) = Y_(:,:,j) + obser_noise_adj;
end
%truth.XI=XI;
truth.q_true = permute(erp_ecog_seg, [2 1 3]);
truth.Y = Y_seg;

clear erp_ecog_seg
clear Y_seg

%constants
Niter=1000;
Ns=5;
Np=500;
P=15;
whiten_flag=0;
roi_flag = 1; % if 1 draw initial set of particles from set of source points within anatomical ROIS (ex 8 ROIS in BBCB simu, 
              %else from whole source space.
dist_thr = 20; %in mm, min distance between source dipoles in each particle
gpu_flag=1;

kappa = 1; b = 299; 
gamma = zeros(1,Niter);gamma(1:2) = 1; gamma(3:b) = 0.95;
gamma(b+1:end) = 0.95*(((0:Niter-b-1)+kappa)/kappa).^(-0.7);
lambdamax=10;
num_sources=Ns;
p=P;
while lambdamax >1
A0 = 0.9*eye(num_sources*p);
A0(num_sources+1:num_sources*p, 1:num_sources*(p-1)) = eye(num_sources*(p-1));
A0(num_sources+1:num_sources*p, 1+num_sources*(p-1):end) = zeros(num_sources*(p-1), num_sources);
lambda=eig(squeeze(A0));lambdamax=max(abs(lambda));
end

V_init_range = [0.1 0.9];
V_q0 = diag(V_init_range(1)+(V_init_range(2)-V_init_range(1))*rand(Ns*P,1));
  
V_q0(Ns+1:Ns*P, 1:Ns)=0;
V_q0(1:Ns*P, Ns+1:end)=0;

% init struct
init.A0 =A0;
init.V0 = V_q0;
init.q0 = zeros(Ns*P,1);
init.P0=1*eye(Ns*P);
init.A0 = A0;
init.V_q0 = V_q0;
init.sigma_m0=1;%chol(mean(trace(E(iChan, iChan))));
init.sigma_b0=1;


%prepare model struct
load (sprintf('/models_for_estim/G_estim_%s', subs{arr(sub)}));
G=G(1:3:306,:);
load (sprintf('/models_for_estim/included_sp_estim_%s', subs{arr(sub)}));
load (sprintf('/models_for_estim/anat_ds_estim_%s', subs{arr(sub)}))
model.G=double(G)*1e-9;
model.mesh=mesh_ds;
model.incl_verts=included_sp;
model.anat_parc = anat_ds;
clear mesh_ds
clear included_sp
clear anat_ds
clear G

%opt_params struct
opt_params.Niter =Niter;
opt_params.gamma = gamma;



Y_avg_seg = squeeze(mean(truth.Y,3));

[state, params, LL_]= saem_bbcb_simu(init, model, opt_params, Y_avg_seg, truth.Y, Ns, Np, P, dist_thr, 1, whiten_flag, roi_flag,truth);




path = '/m/nbe/scratch/braintrack/pnas_results_ecog_new_res_mul_sim_2_redone'; 
 
fname = sprintf(strcat(path,'/state_%d_%d_%.1f.mat'),sub_id(arr(sub)),sub, bio_snr);                                                                                                                                      
save(fname, 'state');
 
fname = sprintf(strcat(path,'/params_%d_%d_%.1f.mat'),sub_id(arr(sub)),sub, bio_snr);                                                                                                                                      
save(fname, 'params');

fname = sprintf(strcat(path,'/LL_%d_%d_%.1f.mat'),sub_id(arr(sub)),sub, bio_snr);                                                                                                                                    
save(fname, 'LL_');

fname = sprintf(strcat(path,'/truth_%d_%d_%.1f.mat'),sub_id(arr(sub)),sub, bio_snr);                                                                                                                                      
save(fname, 'truth');
% 
% 
