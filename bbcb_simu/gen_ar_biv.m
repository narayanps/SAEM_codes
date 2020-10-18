function [data, Arsig, x, lambdamax] = gen_ar_biv(N, P,type)
% generates data according to bivariate AR model with unidirectional flow
% from first to second time series
%
% N: number of data-points
% P: order of AR-model
%
% Guido Nolte, 2006-2015
% Stefan Haufe, 2011-2015
%
% g.nolte@uke.de
% stefan.haufe@tu-berlin.de

% If you use this code for a publication, please ask Guido Nolte for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

M = 2; %number of channels;
sigma = 1; %scale of random AR-parameters

N0=1000; %length of ignored start 
r=10;
while r >= 3 % to ensure a good power ratio between signals
lambdamax=10;
while lambdamax > 1 
  Arsig=[];
  for k=1:P
    aloc = zeros(M);
    aloc([1 4]) = -(0.5+0.5*rand(1,2))'; %[randn; randn];
    if type==1
        aloc([2]) =-(0.5+0.5*rand(1,1))'; %randn(1, 1);
    end
    Arsig=[Arsig,aloc];
  end
  E=eye(M*P);AA=[Arsig;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
end

x=1*randn(M,N+N0);
y=x;
for i=P+1:N+N0;
    yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
    y(:,i)=Arsig*yloc+x(:,i);
end
data=y(:,N0+1:end);

x = x(:, N0+1:end);

Arsig = reshape(Arsig, M, M, P);

if norm(data(1,:)) > norm(data(2,:))
        r = norm(data(1,:)) / norm(data(2,:));
    else
        r = norm(data(2,:)) / norm(data(1,:));
end
end
% figure; psd(data(1, :)')
% figure; psd(data(2, :)')





