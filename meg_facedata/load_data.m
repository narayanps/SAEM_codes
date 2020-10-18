function [Y, XI, G, t, mesh, included_sp, anat_parc ] = load_data(sub, type, meg_path)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if sub < 10
    data_path = sprintf(strcat(meg_path,'/sub00%d/matlab'), sub);
    addpath(genpath(data_path))
    if type == 0
        fname_data_u = sprintf('Y_trials_unfamiliar_sub00%d', sub);
        Yu=load(fname_data_u); 
        Y=Yu.Y;
    elseif type == 1
        fname_data_f = sprintf('Y_trials_famous_sub00%d', sub);
        Yf=load(fname_data_f); 
        Y=Yf.Y;
    else
        fname_data_s = sprintf('Y_trials_scrambled_sub00%d', sub);
        Ys = load(fname_data_s);
        Y = Ys.Y;
    end
    fname_XI = sprintf('XI_sub00%d', sub);
    fname_G = sprintf('G_sub00%d', sub);
    fname_t = sprintf('t_sub00%d', sub);
    fname_mesh = sprintf('mesh_sub00%d', sub);
    fname_isp = sprintf('included_sp_sub00%d', sub);
    fname_anat = sprintf('anat_ds_sub00%d', sub);
else
    data_path = sprintf(strcat(meg_path,'/sub0%d/matlab'), sub);
    addpath(genpath(data_path))
    if type == 0
        fname_data_u = sprintf('Y_trials_unfamiliar_sub0%d', sub);
        Yu=load(fname_data_u);
        Y=Yu.Y;
    elseif type == 1
        fname_data_f = sprintf('Y_trials_famous_sub0%d', sub);
        Yf=load(fname_data_f);
        Y=Yf.Y;
    else
        fname_data_s = sprintf('Y_trials_scrambled_sub0%d', sub);
        Ys = load(fname_data_s);
        Y = Ys.Y;
    end
    fname_XI = sprintf('XI_sub0%d', sub);
    fname_G = sprintf('G_sub0%d', sub);
    fname_t = sprintf('t_sub0%d', sub);
    fname_mesh = sprintf('mesh_sub0%d', sub);
    fname_isp = sprintf('included_sp_sub0%d', sub);
    fname_anat = sprintf('anat_ds_sub0%d', sub);
    
end

data = load(fname_XI);
XI = data.XI;
data = load(fname_G) ;
G = double(data.G);
data = load(fname_t);
t = data.t;
data = load(fname_mesh);
mesh = data.mesh{1,1};
data = load(fname_isp);
included_sp = data.included_sp;
data=load(fname_anat);
anat_parc = data.anat_ds;




