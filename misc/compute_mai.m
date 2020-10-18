function eta = compute_mai(H, noise_cov, data_cov)
inv_noise_cov = inv(noise_cov);
inv_data_cov = inv(data_cov);
eta = trace((H'*inv_noise_cov*H) * inv(H'*inv_data_cov*H));
end

