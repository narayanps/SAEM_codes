#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 15:15:21 2019

@author: puthann1
"""
import mne
import os.path as op
import numpy as np 
import scipy.io as io
res_path = ''##UPDATE_WITH_THE_PATH_WHERE_DATA_FETCH_DATA.PY_IS_SAVED
for sub in range(1,20):
    subj_id = 'sub%03d'%(sub)
    spacing = 'oct6'
    preprocessed_data=op.join(res_path, 'MEG', subj_id)
    fname_src = op.join(res_path, 'subjects', subj_id, 'bem', '%s-%s-src.fif'
                        % (subj_id, spacing))
    fname_inv = op.join(preprocessed_data, 'inverse', '%s-meg-eeg-%s-inv.fif'
                        % (subj_id, spacing))

    patch_stats = True  # include high resolution source space
    src = mne.read_source_spaces(fname_src, patch_stats=patch_stats)
    inv = mne.minimum_norm.read_inverse_operator(fname_inv)

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

    dir_name = op.join(preprocessed_data, 'matlab')
    path_sub = dir_name + '/'
    io.savemat(path_sub + 'mesh_%s.mat'%(subj_id), {'mesh':mesh_ds} )
    io.savemat(path_sub + 'included_sp_%s.mat'%(subj_id), {'included_sp':included_sp} )