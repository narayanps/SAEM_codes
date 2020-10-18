function [ y, noise, C ] = add_meas_noise_trials(noisefreeData, snr)
%author : Narayan Subramaniyam
%Aalto/NBE
% resonable to use this if you have just one type of sensors
J = size(noisefreeData,3);
[m] = size(noisefreeData,1);
T = size(noisefreeData,2);
cov_ = zeros(m,m);
for j=1:1:size(noisefreeData,3)
    cov_ = cov_ + squeeze(noisefreeData(:,:,j))*squeeze(noisefreeData(:,:,j))';
end
cov_ = cov_./(J*T);
avg_var = trace(cov_)/m;


C=(avg_var / snr) * eye(m);
for j=1:1:J
    noise(:,:,j) = chol(C) * randn(m,T);
end


for j=1:1:J
    y(:,:,j) = noisefreeData(:,:,j) + noise(:,:,j);
end

end


