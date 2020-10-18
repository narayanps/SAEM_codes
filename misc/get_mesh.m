function [mesh_ds] = get_mesh(lh, rh)

rh_use_tri = rh.use_tris + 1;
rh_use_v = rh.vertno + 1;

lh_use_tri = lh.use_tris + 1;
lh_use_v = lh.vertno + 1;

offset = size(lh_use_v,2);




for i=1:1:length(lh_use_tri)
faces_l(i,1) = find(lh_use_v == lh_use_tri(i,1));
faces_l(i,2) = find(lh_use_v == lh_use_tri(i,2));
faces_l(i,3) = find(lh_use_v == lh_use_tri(i,3));

faces_r(i,1) = find(rh_use_v == rh_use_tri(i,1)) + offset;
faces_r(i,2) = find(rh_use_v == rh_use_tri(i,2)) + offset;
faces_r(i,3) = find(rh_use_v == rh_use_tri(i,3)) + offset;
end




vertices_l=lh.rr(lh_use_v,:);
vertices_r=rh.rr(rh_use_v,:);
nn_l = lh.nn(lh_use_v,:);
nn_r = rh.nn(rh_use_v,:);

mesh_ds.p = [vertices_l; vertices_r];
mesh_ds.e= [faces_l; faces_r];
mesh_ds.nn = [nn_l ; nn_r];
end

