function run_saem_for_MEG(i, meg_path, save_path)
clc
%close all
%run joint estimation code for all 16 subjects (for face stimuli) and save
%the results.
addpath(genpath('/m/nbe/scratch/braintrack/pnas_codes/joint_est'))
type=0;
subs = {'sub002', 'sub003', 'sub004','sub006','sub007','sub008','sub009'...
     'sub010', 'sub011','sub012','sub013','sub014','sub015','sub017'...
     'sub018', 'sub019'};
sub_id = [2:4 6:15 17:19];
vals = [1:length(sub_id)];
arr = repmat(sub_id, [10 1]); %run 10 random init for each subject. So for every subject 10 results are available and the one with highest LL can be chosen

[status, seed] = system('od /dev/urandom --read-bytes=4 -tu | awk ''{print $2}''');
seed=str2double(seed);
rng(seed);
sd = rng;

%path for func
 all = 1:306;
 mag=1:3:306;
 grad = setdiff(1:306,mag);
 iChan = mag;
 [Y, E, G, t, mesh, included_sp, anat_parc ] = load_data(arr(i), type, meg_path);
 %if you want to use whitened data, run prep_whitner to obtain W and multiply Gmag and Y with W 
 %also set whiten_flag=1 below
 
 Gmag=G(mag,:);
 Gw=Gmag;
trials_full = size(Y,1);
rand_=randperm(trials_full);
num_tr=100;
for j=1:1:num_tr %size(Y,1)
    Yw(j,:,:) = squeeze(Y(rand_(j),iChan,:));
end
trials = rand_(1:num_tr);

mesh.p = mesh.p * 1000 ; % convert to mmema
Y = permute(Yw, [2 3 1]);
Y_avg = squeeze(mean(Y,3));

%SEGMENT DATA
t_beg = 0;
t_end = 0.25;
id_beg = find(abs(t - t_beg) == min(abs(t - t_beg)));
id_end = find(abs(t - t_end) == min(abs(t - t_end)));
Y_avg_seg  = Y_avg(:, id_beg : id_end);
Y_seg = Y(:, id_beg : id_end, :);

%constants
Niter=1000;
Ns=5;
Np=500;
P=10;
scale=1e-9;
whiten_flag=0;
roi_flag = 1; % if 1 draw initial set of particles from set of source points within anatomical ROIS (ex 8 ROIS in BBCB simu, 
              %else from whole source space.
dist_thr = 30; %in mm, min distance between source dipoles in each particle
gpu_flag=1;

kappa = 1; b = 299; 
gamma = zeros(1,Niter);gamma(1:2) = 1; gamma(3:b) = 0.95;
gamma(b+1:end) = 0.95*(((0:Niter-b-1)+kappa)/kappa).^(-0.7);
% 
% mean_A0 =[ -5 -5 -5 -5];
% Sigma_A0 = 1e-3*eye(Ns);
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
init.sigma_m0=chol(mean(trace(E(iChan, iChan))));
init.sigma_b0=5;


%prepare model struct
model.G=double(Gw)*scale;
model.mesh=mesh;
model.incl_verts=included_sp;
model.anat_parc = anat_parc;


%opt_params struct
opt_params.Niter =Niter;
opt_params.gamma = gamma;


[state, params, LL_]= saem(init, model, opt_params, Y_avg_seg, Y_seg, Ns, Np, P, dist_thr, gpu_flag, whiten_flag, roi_flag,[]);



fname = sprintf(strcat(save_path,'/state_%d_%d.mat'),(arr(i)),i); 

save(fname, 'state');

fname = sprintf(strcat(save_path,'/params_%d_%d.mat'),(arr(i)),i);
save(fname, 'params');

fname = sprintf(strcat(save_path,'/LL_%d_%d.mat'),(arr(i)),i); 

save(fname, 'LL_');

fname = sprintf(strcat(save_path,'/sd_%d_%d.mat'),(arr(i)),i);

save(fname, 'sd');

fname = sprintf(strcat(save_path,'/trials_%d_%d.mat'),(arr(i)),i);

save(fname, 'trials');
