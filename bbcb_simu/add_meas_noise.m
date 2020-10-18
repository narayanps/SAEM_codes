function [ y, noise, C, est_snr ] = add_meas_noise(noisefreeData, snr)
%author : Narayan Subramaniyam
%Aalto/NBE
% resonable to use this if you have just one type of sensors

[m] = size(noisefreeData,1);

sig_var = zeros(m,1);

for i = 1:1:m
    sig_var(i,1) = var(noisefreeData(i,:));
end
avg_var = mean(sig_var);


C=(avg_var / snr) * eye(m);

noise = chol(C) * randn(size(noisefreeData));
est_snr = 0;

for i = 1:1:size(noisefreeData,1)
    est_snr = est_snr + (var(noisefreeData(i,:)) / var(noise(i,:)));
end

est_snr = est_snr/m ;

y = noisefreeData + noise;


end


