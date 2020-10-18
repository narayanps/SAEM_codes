function [eta] = compute_ai(h, noise_cov, data_cov)
eta = (h' * inv(noise_cov)* h )/ (h' * inv(data_cov)* h );
end

