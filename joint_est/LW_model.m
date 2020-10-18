function [ind_new, Vr] = LW_model(loc, delta, p)
n = size(loc,1);
N = size(loc,2);
nverts=size(p,1);
a = (3*delta - 1)/(2*delta);
h2 = 1 - a^2;
mean_ = a.*loc + repmat((1-a)*mean(loc,2), [1 N]);
Vt = (var(loc,0,2)'+ 1e-8).*eye(n);
Vr=h2*Vt;
loc = mean_ + chol(Vr)*randn(n,N);

% h2 = 1-((3*c-1)/2*c)^2;
% a = sqrt(1-h2);
% mean_ = a.*loc + repmat((1-a)*mean(loc,2), [1 N]);
% V = (var(loc')+0.1).*eye(n);
% V_r = h2*V;
% loc = mean_ + chol(V_r)*randn(n,N);


%V_r = diag(((1./eta) - 1) * (var(loc,0,2)'+1e-7));
%loc = loc + chol(V_r)*randn(n,N);
ind_new=zeros(n/3, N);
for i=1:1:(n/3)
    for j=1:1:N
 tmp = sqrt(sum((p - repmat(loc(3*(i-1) + 1:3*i, j)', ...
     nverts,1)).^2,2));
ind_new(i,j)=find(tmp == min(tmp));
    end
end
