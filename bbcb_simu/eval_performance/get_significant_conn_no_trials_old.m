function [PDC_, PDC, PDCsth, f] = get_significant_conn_no_trials_old(x,Am,S, numsurr, nfft, fc, Ns, p)
%addpath(genpath('/m/nbe/scratch/braintrack/pnas_codes/misc/biosig_ext'));
%generate surrogates
S1=S .* eye(Ns);
for is=1:numsurr
%     xs=surrVFT(x);
%     [Ams,Sus,Yp,Up]=idMVAR(xs,p,1);
%    [~,~,~,gpdcs,~,~,~,~,~,~,~] = fdMVAR(Ams,Sus,nfft,fc);
%    PDCs(:,:,:,is)=abs(gpdcs(:,:,:)).^2;
    % surrogates for DC and PDC
    for ii=1:Ns
        for jj=1:Ns
            if ii~=jj
                xs=surrVCFTd(x,Am,S1,ii,jj);
                [Ams,Sus]=idMVAR(xs,p,2);
                Sus = Sus .* eye (Ns);
                [~,~,~,gpdcs,~,~,~,~,~,~,~] = fdMVAR(Ams,Sus,nfft,fc);
                PDCs(ii,jj,:,is)=abs(gpdcs(ii,jj,:)).^2;
            end
        end
    end
end

PDCsth=prctile(PDCs,95,4);

[~,~,~,gpdc,~,~,~,~,~,~,f] = fdMVAR(Am,S1,nfft,fc);
PDC=abs(gpdc).^2; % partial directed coherence
%PDC=cohm.pdcspctrm;
% figure
% M=5;
% for i=1:1:M
%     for j=1:1:M
% subplot(M,M, M*(i-1)+j)
% %if i ~=j
% plot(f,squeeze(PDC(i,j, :)))
% % hold on
% % plot(f, squeeze(PDCsth(i,j, :)), '-.')
% %axis([0 fc/2 -0.05 1.05]);
% %ylim([0 1])
% %end
% xlabel('Frequency [Hz]')
% ylabel('PDC')
%     end
% end

PDC(PDC < PDCsth) = 0;
PDC_=squeeze(mean(PDC(:,:,:),3));
pdc_surr = mean(PDCsth(:,:,:),3);
PDC_(PDC_ < pdc_surr) = 0;
end