function [TP, FP, TN, FN, cle, sle, pdc_int, success_flag] = compute_error_stats_no_trial(est,truth, TP, FP, TN, FN)
success_flag=0;
truth.A_true = reshape(truth.model.A, 2,4);
numsurr=100;
nfft = 128;
fc = 100;
P = 2;
[~,~,~,gpdc,~,~,~,~,~,~,f] = fdMVAR(truth.A_true,eye(2),nfft,fc);
PDC_true=abs(gpdc).^2; % partial directed coherence
PDC_est = [];


[comb_id, ~, ~] = check_loc_error(truth.r(:,1), est.r, 2) ;

if (comb_id(1) == 1 && comb_id(2)  == 2)
    sle(1) = norm(est.r(1:3) - truth.r(1:3,1));
    sle(2) = norm(est.r(4:6) - truth.r(4:6,1));
    if max(sle) <= 20
        success_flag = 1;
        [pdc_int, PDC_est, PDCsth, f] = get_significant_conn_no_trials_old(est.amp(:,1:end-1),est.A,est.V, numsurr, nfft, fc, 2, P);
        %[pdc_int] = get_significant_conn_no_trials(est.amp(:,1:end-1),est.A,est.V);
        if truth.type==1
            if pdc_int(2,1) > 0
                TP=TP+1;
            end
            if pdc_int(2,1) == 0
                FN=FN+1;
            end
            
            if pdc_int(1,2) > 0
                FP=FP+1;
            end
            if pdc_int(1,2) == 0
                TN=TN+1;
            end
        elseif truth.type==0
            if (pdc_int(2,1)==0 && pdc_int(1,2)==0)
                TN=TN+1;
            end
            
            
            if (pdc_int(2,1) > 0 || pdc_int(1,2) > 0)
                FP=FP+1;
            end
            
        end
    end
elseif (comb_id(1) == 2 && comb_id(2)  == 1)
    
    sle(1) = norm(est.r(1:3) - truth.r(4:6,1));
    sle(2) = norm(est.r(4:6) - truth.r(1:3,1));
    if max(sle) <= 20
        success_flag = 1;
        AA=[];
        for p=1:1:P
            tmp = flip(reshape(est.A(:,2*(p-1)+1:2*p),4,1));
            AA = [AA reshape(tmp',2,2)];
        end
        est.A = AA;
        
        est.amp = flip(squeeze(est.amp));
        tmp=est.V;
        tmp(1,1) = est.V(2,2);
        tmp(2,2) = est.V(1,1);
        est.V=tmp;
        
        [pdc_int, PDC_est] = get_significant_conn_no_trials_old(est.amp(:,1:end-1),est.A,est.V, numsurr, nfft, fc,2, P);
        %[pdc_int] = get_significant_conn_no_trials(est.amp(:,1:end-1),est.A,est.V);
        if truth.type==1
            if pdc_int(2,1) > 0
                TP=TP+1;
            end
            if pdc_int(2,1) == 0
                FN=FN+1;
            end
            if pdc_int(1,2) > 0
                FP=FP+1;
            end
            if pdc_int(1,2) == 0
                TN=TN+1;
            end
        elseif truth.type == 0
            if (pdc_int(2,1)==0 && pdc_int(1,2)==0)
                TN=TN+1;
            end
            
            
            if (pdc_int(2,1) > 0 ||pdc_int(1,2) > 0)
                FP=FP+1;
            end
            
            
        end
    end
end

if ~isempty(PDC_est)
    tmp = norm(squeeze(mean(PDC_true(:,:,:),3)) - pdc_int, 'fro')/norm(squeeze(mean(PDC_true,3)), 'fro');
    cle = mean(tmp);
else
    cle = inf;
    [~, ~, sle] = check_loc_error(truth.r(:,1), est.r, 2) ;
    pdc_int=[];
end

