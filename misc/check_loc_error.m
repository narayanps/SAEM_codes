function [comb_id, min_mean_err, err] = check_loc_error(loc_true, loc_est, nq)    

    aa=perms(1:nq);
    for jj=1:1:length(aa) 
     for i=1:1:nq 
        dip_err(jj,i) = sqrt(sum((loc_est((3*i-3)+1:3*i,end) -  loc_true((3*aa(jj,i)-3)+1:3*aa(jj,i),end)).^2));
     end
    end
    
    mean_err = mean(dip_err,2);
    [min_mean_err, ind] = min(mean_err);
    id = find(mean_err == min(mean_err));
    comb_id = aa(id,:);
    err = dip_err(ind,:);