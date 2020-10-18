function  run_bbcb_sim(sim, type, gpu_flag, res_path)
meas_snr=5;
addpath(genpath(fullfile(pwd, '..', 'joint_est')));
addpath(genpath(fullfile(pwd, '..', 'misc')));
addpath('data')

%load head model
load miscdata
load sa

bio_snr_vals = [1 3 5 10];
bio_snr_arr = repmat(bio_snr_vals , [50 1]);



%load mesh details
mesh.p= sa.cortex75K.vc(sa.cortex2K.in_from_cortex75K, :);
mesh.e=sa.cortex2K.tri;
mesh.nn = sa.cortex75K.normals(sa.cortex2K.in_from_cortex75K, :) ; 

%get led-field matrix
G = sa.cortex75K.MEG_V_bem_normal(:, sa.cortex2K.in_from_cortex75K);

% set dipole locations from pre-defined 8 ROIS
for i=1:1:length(inds_roi_outer_2K)
    anat_parc.parcels{i,1} = inds_roi_outer_2K{1,i};
end

included_sp = 1:length(G);
elind = unique(round(1:2.5:298));
elind = elind(1:108);

%generate MEG
G = G(elind, :);
num_steps=1000;
num_trials=1;
true_num_sources=2;
p_true=2;

%make results reproducible
[status, seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(seed);
sd = rng;

% sample dipole locations from pre-defined 8 ROIS
for i=1:1:length(inds_roi_inner_2K)
rois_inner.parcels{i,1} = inds_roi_inner_2K{1,i};
end

[truth] = generate_GT_dip_sources(G, mesh, [1:size(mesh.p,1)]', num_steps,true_num_sources, meas_snr,bio_snr_arr(sim),p_true,rois_inner,type);
truth.sd=sd;

%SAEM iteration parameters
num_iter = 400;
b=199;
kappa = 1; 
gamma = zeros(1,num_iter);
gamma(1:2) = 1;
gamma(3:b) = 0.98;
gamma(b+1:end) = 0.98*(((0:num_iter-(b+1))+kappa)/kappa).^(-0.7);

%set number of sources, particles, order etc
Ns=2;
Np=500;
P=2;
scale=1;
whiten_flag=0;
roi_flag = 1; % if 1 draw initial set of particles from set of source points within anatomical ROIS (ex 8 ROIS in BBCB simu, 
              %else from whole source space.
dist_thr = 10; %in mm, min distance between source dipoles in each particle


%initial value for A and V
mean_A0 =[0 0];
Sigma_A0 = eye(Ns);
lambdamax=10;
while lambdamax >1       
for jj = 1:P
    A_p = mean_A0' + chol(Sigma_A0)*randn(Ns,1);         %diagonal entries of A matrix are U[a,b]
    A_p = diag(diag(A_p));
    A0(Ns*(jj-1)+1:Ns*(jj-1)+Ns,1:Ns) = A_p.*eye(Ns);
end 
A0 = A0';
A0(Ns+1:Ns*P, 1:Ns*(P-1)) = eye(Ns*(P-1));
A0(Ns+1:Ns*P, 1+Ns*(P-1):end) = zeros(Ns*(P-1), Ns);
lambda=eig(squeeze(A0));lambdamax=max(abs(lambda));
end  

V_init_range = [0.1 0.9];
V_q0 = diag(V_init_range(1)+(V_init_range(2)-V_init_range(1))*rand(Ns*P,1));
  
V_q0(Ns+1:Ns*P, 1:Ns)=0;
V_q0(1:Ns*P, Ns+1:end)=0;

% init struct for SAEM algorithm
init.A0 =A0;
init.V0 = V_q0;
init.q0 = zeros(Ns*P,1);
init.P0=1*eye(Ns*P);
init.sigma_m0 = 1; %eye(size(G,1));
init.sigma_b0 = 1;

% model struct for SAEM
model.G=double(G)*scale;
model.mesh=mesh;
model.incl_verts=included_sp;
model.anat_parc=anat_parc;


%opt_params struct
opt_params.Niter = num_iter;
opt_params.gamma = gamma;


Y=truth.Y(:,1:300,:); %use only 300 points for SAEM
Y_avg = squeeze(mean(Y,3)); 

[state, params, LL_complete]= saem(init, model, opt_params, Y_avg, Y, Ns, Np, P, dist_thr, gpu_flag, whiten_flag, roi_flag,truth);

path = res_path; 
 
fname = sprintf(strcat(path,'/state_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'state');
 
fname = sprintf(strcat(path,'/params_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'params');

fname = sprintf(strcat(path,'/truth_snr_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'truth');

fname = sprintf(strcat(path,'/LL_%.1f_type_%d_sim_%d.mat'),bio_snr_arr(sim),type,sim);                                                                                                                                      
save(fname, 'LL_complete');




