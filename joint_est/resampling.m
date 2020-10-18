function [indices] = resampling( weights, nParticles, alg )
% Narayan Subramaniyam, 07-2016
% narayan.subramaniyam@aalto.fi

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

% n - dimension of the state-space
% N - no. of particles
% m0 - initial mean
% P0 - initial covariance
%indices = ones(1,nParticles);

m = 0:nParticles-1;

if strcmp(alg, 'simplerandom')
    u = (m + rand(1,1)) / nParticles;
    c = cumsum(weights);
    c(length(weights)) = 1.0;
    j = 1;
    
    for i = 1:1:nParticles
        while(c(j) < u(i))
            j = j+1;
        end
        indices(i) = j;
    end
    
elseif strcmp(alg, 'stratified')
    u = (m + rand(1,nParticles)) / nParticles;
    c = cumsum(weights);
    i = 1;
    j = 1;
    while (i <= nParticles)
        if u(i) < c(j)
           indices(i) = j;
           i = i + 1;
        else
            j = j + 1;
        end
    end
    

elseif strcmp(alg, 'systematic')
    u = (m + rand(1,1)) / nParticles;
    c = cumsum(weights);
    i = 1;
    j = 1;
    while (i <= nParticles)
        if u(i) < c(j)
           indices(i) = j;
           i = i + 1;
        else
            j = j + 1;
        end
    end
elseif strcmp(alg, 'residual')
    
end
    
end
