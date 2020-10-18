function [ind_, r_, G_] = set_part_traj_noroi(pnts, G, dist_thr, num_sources, T,rois_inner, roi_flag)

if roi_flag==1
source_ids = [];
id = zeros(num_sources,1);
for i=1:1:length(rois_inner.parcels)
   source_ids = [source_ids; rois_inner.parcels{i,1}];
end
else
source_ids = 1:length(pnts);
end

ind_ = zeros(num_sources, T);
distMat = inf*eye(num_sources,num_sources);
dmin = 0;
nr = 3*num_sources;
%

while dmin < dist_thr % in meters
    for jj=1:1:num_sources
    id(jj) = source_ids(ceil(length(source_ids)*rand(1)));
    end

    for i=1:1:num_sources
        d= pnts(id, :) - repmat(pnts(id(i), :), num_sources,1);
        for j=1:1:num_sources
            if i ~= j
                distMat(i,j) = norm(d(j,:));
            end
        end
    end
    dmin = min(distMat(:));
end

for j=1:1:num_sources
    ind_(j,:) = repmat(id(j), [1 T]);
end

r_ = repmat(reshape(pnts(ind_(:,1),:)',nr,1),[1 T]);
G_ =repmat(G(:, ind_(:,1)),[1 1 T]);


end