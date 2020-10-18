#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:00:56 2020

@author: Narayan Subramaniyam / AALTO / NBE
"""
import mne
import os.path as op
from mne.minimum_norm import apply_inverse, read_inverse_operator, write_inverse_operator, make_inverse_operator
import scipy.io as io
from mne import convert_forward_solution
import numpy as np

###PLEASE UPDATE THIS AS PER YOUR PATHS####################################
your_path = ''#'' ##UPDATE_WITH_THE_PATH_WHERE_DATA_FETCH_DATA.PY_IS_SAVED
###########################################################################
path_save = '/models_for_estim/'
subjects_dir = op.join(your_path, 'subjects')
data_path = op.join(subjects_dir, 'MEG')
dir_path=op.join(your_path, 'ds117')
conductivity = (0.3,)  # for single layer
for subject_id in range(1,20):
    subject = "sub%03d" % subject_id
    preprocessed_data = op.join(your_path, 'MEG', subject)
    raw_data_path =  op.join(dir_path, subject, 'MEG')
    spacing = 'oct6'
    model = mne.make_bem_model(subject, ico=3,
                           conductivity=conductivity,
                           subjects_dir=subjects_dir)
    fname_trans = op.join(raw_data_path, '%s-trans.fif' % subject)
    fname_src = op.join(your_path, 'subjects', subject, 'bem', '%s-%s-src.fif'
                        % (subject, spacing))
    l_cutoff=0.25
    h_cutoff=40
    fname_ave = op.join(preprocessed_data, 'evoked', '%s_highpass_%sHz-ave.fif' %(subject, l_cutoff))
    bem = mne.make_bem_solution(model)
    info = mne.io.read_info(fname_ave)
    fname_cov = op.join(preprocessed_data, 'cov', '%s_highpass-%sHz-cov.fif'
                            % (subject, l_cutoff))
# Because we use a 1-layer BEM, we do MEG only
    fwd = mne.make_forward_solution(info, fname_trans, fname_src, bem,
                                    meg=True, eeg=False, mindist=5)
    cov = mne.read_cov(fname_cov)
    inv = make_inverse_operator(
        info, fwd, cov, loose=0.2, depth=0.8)
    fwd_n = convert_forward_solution(fwd, force_fixed=True, surf_ori=True, use_cps=True)
    
    patch_stats = True  # include high resolution source space
    src = mne.read_source_spaces(fname_src, patch_stats=patch_stats)
    

    nl_p = src[0]['vertno'].shape[0]
    nr_p = src[1]['vertno'].shape[0]
    nl_e = src[0]['use_tris'].shape[0]
    nr_e = src[1]['use_tris'].shape[0]
    mesh = [{'p': np.zeros((nl_p, 3), dtype=float), 'e': np.zeros((nl_e, 3), dtype=int)}, 
        {'p': np.zeros((nr_p, 3), dtype=float), 'e': np.zeros((nr_e, 3), dtype=int)}]
    mesh_ds = [{'p' : np.zeros((nl_p+nr_p, 3), dtype=float),'e': np.zeros((nl_e + nr_e, 3), dtype=int)}]

    offset = src[0]['vertno'].shape[0]
#add 1 in the end to be compatible with MATLAB
    for i in range(nl_e):
        mesh[0]['e'][i,0] = np.where(src[0]['vertno'] == src[0]['use_tris'][i,0])[0] +1
        mesh[0]['e'][i,1] = np.where(src[0]['vertno'] == src[0]['use_tris'][i,1])[0] +1
        mesh[0]['e'][i,2] = np.where(src[0]['vertno'] == src[0]['use_tris'][i,2])[0] +1
        mesh[1]['e'][i,0] = np.where(src[1]['vertno'] == src[1]['use_tris'][i,0])[0] + offset +1
        mesh[1]['e'][i,1] = np.where(src[1]['vertno'] == src[1]['use_tris'][i,1])[0] + offset +1
        mesh[1]['e'][i,2] = np.where(src[1]['vertno'] == src[1]['use_tris'][i,2])[0] + offset +1
    
    mesh[0]['p'] = src[0]['rr'] [src[0]['vertno'],:]
    mesh[1]['p'] = src[1]['rr'] [src[1]['vertno'],:]
    mesh_ds[0]['p'] = np.concatenate((mesh[0]['p'], mesh[1]['p']))
    mesh_ds[0]['e'] = np.concatenate((mesh[0]['e'], mesh[1]['e']))

    src_ = inv['src']
    orig, incl_pts_l, id_ =np.intersect1d(src[0]['vertno'], src_[0]['vertno'], return_indices=True)
    orig, incl_pts_r, id_ =np.intersect1d(src[1]['vertno'], src_[1]['vertno'], return_indices=True) 
    incl_pts_r = incl_pts_r + offset
    included_sp = np.concatenate((incl_pts_l, incl_pts_r))
    included_sp = included_sp + 1 #MATLAB
    io.savemat(path_save + 'src_left_estim_%s.mat'%(subject), {'src':src_[0]} )
    io.savemat(path_save + 'src_right_estim_%s.mat'%(subject), {'src':src_[1]} )

    #io.savemat(path_save + 'mesh_estim_%s.mat'%(subject), {'mesh':mesh_ds} )
    io.savemat(path_save + 'included_sp_estim_%s.mat'%(subject), {'included_sp':included_sp} )
    io.savemat(path_save + 'G_estim_%s.mat'%(subject), {'G':fwd_n['sol']['data']} )
