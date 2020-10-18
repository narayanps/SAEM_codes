function [v_id_mne, hemi, v_id] = map_locs_to_mne(avg_loc, mesh, incl_vert, src_l, src_r)
num_sources = size(avg_loc,1)./3;
v_id = zeros(num_sources,1);
m=size(mesh.p(incl_vert,:),1);
for i=1:1:num_sources
        d = sqrt(sum((mesh.p(incl_vert,:) - repmat(avg_loc((i-1)*3 + 1 : 3*i, 1)', m, 1)).^2, 2));
        v_id(i,1) = find(d==min(d));
    end
    v_id_full = incl_vert(v_id)';
    v_id_mne = zeros(size(v_id));
    hemi = zeros(size(v_id));
    
    
    for i=1:1:size(v_id,1)
        if v_id_full(i,1) <= 4098
            v_id_mne(i,1) = src_l.src_full.vertno((v_id_full(i,1)));
            hemi(i) = 0;
        else
            v_id_mne(i,1) = src_r.src_full.vertno((v_id_full(i,1)-4098));
            hemi(i)=1;
        end
    end
end

