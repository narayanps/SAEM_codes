function [id, err] = check_loc_error_ds(loc_true, loc_est)    

    ne = length(loc_est)./3;
    nt = length(loc_true)./3;
    for j=1:1:nt
        loc_t = loc_true((3*j-3)+1:3*j);
        for k=1:1:ne
            loc_e = loc_est((3*k-3)+1:3*k);
            err(k,j) = norm(loc_t-loc_e); %sqrt(sum((loc_t-loc_e).^2));
        end
    end
    
%     [err,id] = min(err,[],1);
[err1,id1] = min(err(:,1));
err(id1,:)=inf;
[err2,id2] = min(err(:,2));

err = [];
err = [err1 err2];
id= [id1 id2];