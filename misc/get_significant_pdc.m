function [pdc_delta, pdc_theta, pdc_alpha, pdc_beta, PDC,f] = get_significant_pdc(x,Am,S, numsurr, nfft, fc, Ns, p)

%generate surrogates

for is=1:numsurr
    is
    % surrogates for DC and PDC
    for ii=1:Ns
        for jj=1:Ns
            if ii~=jj
%                 clc; disp(['surrogate ' int2str(is) ' of ' int2str(numsurr)
%                     ', ' int2str(ii) '->' int2str(jj) ', CFTd surrogates']);
                for tr=1:1:size(x,3)
                    xs(:,:,tr)=surrVCFTd(squeeze(x(:,:,tr)),Am,S,ii,jj);
                  %  [Ams(:,:,tr),Sus(:,:,tr)]=idMVAR(xs,p,2);
                  %  Sus(:,:,tr) = Sus(:,:,tr) .* eye(Ns);
                end
                ys = permute (xs, [2 1 3]);
                [~, Ams, Sus, ~, ~, ~]=arfit(ys,p,p);
                [~,~,~,gpdcs,~,~,~,~,~,~,~] = fdMVAR(Ams,Sus,nfft,fc);
                PDCs(ii,jj,:,is)=abs(gpdcs(ii,jj,:)).^2;
            end
        end
    end
end

PDCsth=prctile(PDCs,99,4);

[~,~,gpdc,~,~,~,~,~,~,~,f] = fdMVAR(Am,S,nfft,fc);
PDC=abs(gpdc).^2; % partial directed coherence

% 
%figure
% count=0;
% for i=1:1:Ns
%     for j=1:1:Ns
%         count=count+1;
%         subplot(Ns,Ns,(Ns*i-Ns)+j);
%         plot(f, squeeze(PDC(i,j,:)),'g');  hold on;
%         %plot(f, squeeze(PDCsth(i,j,:)),'g:', 'LineWidth', 2.5);
%         axis([0 fc/2 -0.05 1.05]); %title(title_1(count));
%     end
% end

alpha = find(f > 8 & f < 13);
pdc_alpha = mean(PDC(:,:,alpha),3);
pdcs = mean(PDCsth(:,:,alpha),3);
pdc_alpha(pdc_alpha < pdcs) = 0;

theta = find(f > 4 & f < 8);
pdc_theta = mean(PDC(:,:,theta),3);
pdcs = mean(PDCsth(:,:,theta),3);
pdc_theta(pdc_theta < pdcs) = 0;

delta = find(f > 0.5 & f < 4);
pdc_delta = mean(PDC(:,:,delta),3);
pdcs = mean(PDCsth(:,:,delta),3);
pdc_delta(pdc_delta < pdcs) = 0;


beta = find(f > 13 & f < 30);
pdc_beta = mean(PDC(:,:,beta),3);
pdcs = mean(PDCsth(:,:,beta),3);
pdc_beta(pdc_beta < pdcs) = 0;




% pdc = mean(PDC(:,:,:),3);
% pdcs = mean(PDCsth(:,:,:),3);
% pdc(pdc < pdcs) = 0;

