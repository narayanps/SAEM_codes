function [id, r_, G_] = set_all_part_traj_noroi(pnts, G, dist_thr, num_sources, N, rois_inner, roi_flag)

if roi_flag==1
source_ids = [];
id = zeros(num_sources,1);
for i=1:1:length(rois_inner.parcels)
   source_ids = [source_ids; rois_inner.parcels{i,1}];
end
else
source_ids = 1:length(pnts);
end

dist_sources=zeros(1,N);
r_=zeros(3*num_sources,N);
G_=zeros(size(G,1),num_sources,N);
min_dist = 0;
while min_dist < dist_thr
    unique_flag=0;
    while unique_flag==0
        id= source_ids(ceil(length(source_ids)*rand(num_sources,N)));
        unique_id = unique(id', 'rows', 'stable');
        if size(unique_id,1) == N
            unique_flag=1;
        end
    end
    
    for j=1:1:N
        X=[];
        for i=1:1:num_sources
            X = [X;pnts(id(i,j),:)];
        end
        dist_sources(1,j) = max(pdist(X));
    end
    
    min_dist = min(dist_sources);
end
min_dist
for j=1:1:N
r_(:,j) = reshape(pnts(id(:,j),:)',3*num_sources,1);
G_(:,:,j) = G(:, id(:,j)');
end
end
