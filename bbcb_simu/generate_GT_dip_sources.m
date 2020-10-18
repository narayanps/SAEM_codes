function [GT] = generate_GT_dip_sources(G, mesh, included_sp, num_steps, num_sources, meas_snr, bio_snr, p, rois_inner, type)


rois_inner_sp = [];
source_vertices = zeros(num_sources,1);
for i=1:1:length(rois_inner.parcels)
   rois_inner_sp = [rois_inner_sp; rois_inner.parcels{i,1}];
end



distMat = inf*eye(num_sources,num_sources);
dmin = 0;
while dmin < 20 % just to ensure the dipoles drawn are not too close (in mm )
         for i=1:1:num_sources
                source_vertices(i,1) = rois_inner_sp(ceil(length(rois_inner_sp)*rand(1)));
         end
    for i=1:1:num_sources
        d= mesh.p(included_sp(source_vertices),:) - repmat(mesh.p(included_sp(source_vertices(i,1)),:), num_sources,1);
    for j=1:1:num_sources
        if i ~= j
        distMat(i,j) = norm(d(j,:));
        end
    end
    end
    dmin = min(distMat(:));
end

%MVAR dynamics
r_true(:,1:num_steps) = repmat(reshape(mesh.p(included_sp(source_vertices),:)',3*num_sources,1), 1, num_steps);
[q_true, A_true,~, ~] = gen_ar_biv(num_steps, p, type);
GT.type=type;

%generate brain signal
Y_brain_signal = G(:, source_vertices) * q_true(1:num_sources,1:end) ;

%generate brain noise
if ~isempty(bio_snr)
    n_noise_sources=500; 
    ind_noise_rand=randperm(size(G,2));
    ind_noise=ind_noise_rand(1:n_noise_sources)';   
    pn = mkpinknoise(num_steps, n_noise_sources)';
    Y_brain_noise = G(:, ind_noise)*pn;
    
    %n_noise_sources=size(G,2);
    %bn = randn(n_noise_sources, num_steps);
    %Y_brain_noise = G(:, :)*bn;
    bio_snr_=trace(Y_brain_signal*Y_brain_signal')./trace(Y_brain_noise*Y_brain_noise');
    Y_brain_noise = Y_brain_noise*sqrt(bio_snr_/bio_snr);
    Y_ = Y_brain_signal + Y_brain_noise;
    GT.sigma_b = sqrt(bio_snr_/bio_snr);

else
    
    Y_ = Y_brain_signal;
end

[Y, ~, XI, ~ ] = add_meas_noise(Y_,meas_snr);

% save GT
GT.model.V_q = eye(num_sources); % covariance of process noise
GT.model.A=A_true;
GT.Y=Y;
GT.q = q_true;
GT.r = r_true;
GT.r_id =(source_vertices);
GT.nSources=num_sources;
GT.sigma_m = chol(XI(1,1));

