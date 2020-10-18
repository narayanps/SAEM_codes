function [GT] = simulate_MEG(G, mesh, included_sp, num_steps, num_sources, meas_snr, bio_snr, p, type)

distMat = inf*eye(num_sources,num_sources);
dmin = 0;
source_verts = zeros(num_sources,1);
while dmin < 30 % just to ensure the dipoles drawn are not too close (in mm )
source_verts = randperm(length(included_sp),num_sources)';
    for i=1:1:num_sources
        d= mesh.p(included_sp(source_verts'),:) - repmat(mesh.p(included_sp(source_verts(i,1)),:), num_sources,1);
    for j=1:1:num_sources
        if i ~= j
        distMat(i,j) = norm(d(j,:));
        end
    end
    end
    dmin = min(distMat(:));
end

%MVAR dynamics
r_true(:,1:num_steps) = repmat(reshape(mesh.p(included_sp(source_verts'),:)',3*num_sources,1), 1, num_steps);
if type==1
    [q_true, A_true,~, ~] = gen_ar_biv(num_steps, p);
    GT.type=1;
else
    GT.type=0;
    [q_1, A_1,~, ~] = gen_ar_uni(num_steps, p);
    [q_2, A_2,~, ~] = gen_ar_uni(num_steps, p);
    
    q_true(1,:)=q_1;
    q_true(2,:)=q_2;
    A_true = zeros(num_sources,num_sources,p);
    for jj=1:1:p
        A_true(1,1,jj) = A_1(1,1,jj);
        A_true(2,2,jj) = A_2(1,1,jj);
    end
end


%generate brain signal
Y_brain_signal = G(:, source_verts) * q_true(1:num_sources,1:end) ;

%generate brain noise
if ~isempty(bio_snr)
    %generate brain noise
    n_noise_sources=size(G,2);
    bn = randn(n_noise_sources, num_steps);
    Y_brain_noise = G(:, :)*bn;
    bio_snr_=trace(Y_brain_signal*Y_brain_signal')./trace(Y_brain_noise*Y_brain_noise');
    Y_brain_noise = Y_brain_noise*sqrt(bio_snr_/bio_snr);
    Y_ = Y_brain_signal + Y_brain_noise;
    GT.sigma_b = sqrt(bio_snr_/bio_snr);

else
    
    Y_ = Y_brain_signal;
end



[Y, ~, XI, ~ ] = add_meas_noise(Y_,meas_snr);


%generate pseudo-baseline MEG/EEG , i.e. data comprising of brain noise and
%sensor noise only.
[Y_baseline, ~, XI_baseline, ~ ] = add_meas_noise(Y_brain_noise,meas_snr);

% save GT
GT.model.V_q = eye(num_sources); % covariance of process noise
GT.model.A=A_true;
GT.Y=Y;
GT.q = q_true;
GT.r = r_true;
GT.r_id =included_sp(source_verts);
GT.nSources=num_sources;
GT.sigma_m = chol(XI(1,1));
GT.Y_baseline = Y_baseline;
GT.sigma_baseline = chol(XI_baseline(1,1));



%generate brain signal
% Y_brain_signal = G(:, source_verts) * q_true(1:num_sources,1:end) ;
% 
% %generate brain noise
% 
% if ~isempty(bio_snr)
%     %generate pink noise
%     n_noise_sources=size(G,2);
%     %noise_inds = ceil(size(G, 2)*rand(n_noise_sources, 1));
%     %pn = mkpinknoise(num_steps, n_noise_sources)';
%     pn = randn(n_noise_sources, num_steps);
%     Y_brain_noise = G(:, :)*pn;
%     bio_snr_=trace(Y_brain_signal*Y_brain_signal')./trace(Y_brain_noise*Y_brain_noise');
%     Y_brain_noise = Y_brain_noise*sqrt(bio_snr_/bio_snr);
%     Y_ = Y_brain_signal + Y_brain_noise;
% else
%     
%     Y_ = Y_brain_signal;
% end
% 
% [Y, ~, XI, ~ ] = add_meas_noise(Y_,meas_snr);
% sqrt(bio_snr_/bio_snr)
% chol(XI(1,1))
% 
% 
% % save GT
% GT.model.G=G;
% GT.model.XI=XI + ((bio_snr_/bio_snr)*(G*G'));
% GT.model.V_q = eye(num_sources); % covariance of process noise
% GT.model.A=A_true;
% GT.model.mesh=mesh;
% GT.Y=Y;
% GT.q = q_true;
% GT.r = r_true;
% GT.r_id = included_sp(source_verts);
% GT.nSources=num_sources;
% GT.model.incl_vert = included_sp;
% GT.sigma_m = chol(XI(1,1));

