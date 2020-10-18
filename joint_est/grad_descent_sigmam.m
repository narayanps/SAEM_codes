function [sigma_m] = grad_descent_sigmam(S_R, G, T, sigma_b, sigma_m)

diff_obj = inf;
obj = 1e10;
old_obj = 0;
step_init=1.0;
tau=0.8;
iter_count=1;
max_count = 1000;
while abs(diff_obj/obj)>=1e-10
    if iter_count > max_count
        break;
    end
    R = sigma_m^2*eye(size(G,1)) + sigma_b^2*(G*G');
    grad_R = T*inv(R) - R\S_R/R;
    grad_sigmam=trace(grad_R' * (2*sigma_m*(eye(size(G,1)))));
    f = calc_obj(R, S_R, T);
    step=step_init;
    tmp_diff = inf;
    while tmp_diff > 0
        step = step*tau;
        ref = f - (step/2)*grad_sigmam^2;
        tmp_sigmam = sigma_m - step*grad_sigmam;
        R = tmp_sigmam^2*eye(size(G,1)) + sigma_b^2*(G*G');
        tmp_f = calc_obj(R, S_R, T);
        tmp_diff = tmp_f - ref;
    end
    sigma_m = abs(tmp_sigmam);
    old_obj = obj;
    obj = tmp_f;
    diff_obj = old_obj - obj;
    iter_count=iter_count+1;
end



