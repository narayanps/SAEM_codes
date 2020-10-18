function [dist_vec, id] = find_nearest_id(pnts,pos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dist_vec=sqrt(sum((pnts - repmat(pos,[ size(pnts,1) 1])).^2,2));
id = find(dist_vec == min(dist_vec));
end

