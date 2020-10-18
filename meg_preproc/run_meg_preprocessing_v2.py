# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 14:14:22 2017

@author: Narayan
"""
#import mne
import sys
import os
import os.path as op
#from mne.parallel import parallel_func
from preproc_pipeline import MegSetPreProcessing
    

#event ids
#event_id = dict(famous_1=5, famous_2=6, famous_3=7, unknown_1=13, unknown_2=14, unknown_3=15, 
#                    scrambled_1=17, scrambled_2=18, scrambled_3=19)




events_id = {
    'face/famous/first': 5,
    'face/famous/immediate': 6,
    'face/famous/long': 7,
    'face/unfamiliar/first': 13,
    'face/unfamiliar/immediate': 14,
    'face/unfamiliar/long': 15,
    'scrambled/first': 17,
    'scrambled/immediate': 18,
    'scrambled/long': 19,
    }
tsss=None
total_runs=6
subj_id = 'sub%03d'%(int(sys.argv[1]))
print(subj_id)
res_path = '' #UPDATE_WITH_THE_PATH_WHERE_DATA_FETCH_DATA.PY_IS_SAVED
dir_path = op.join(res_path, 'ds117')
min_duration = 0.003
stim_channel = 'STI101'
l_cutoff=0.8
h_cutoff=40
spacing = 'oct6'
data = MegSetPreProcessing(subj_id, dir_path, res_path, events_id, total_runs, l_cutoff, h_cutoff, spacing)
data.extract_events(stim_channel='STI101', min_duration=0.003)
    #if subj_id == 'sub003':
    #    data.maxwell_filter(l_cutoff=None, h_cutoff=h_cutoff, method='fir',fir_design='firwin')
data.run_filter(l_cutoff=l_cutoff, h_cutoff=h_cutoff)
data.run_ica(tsss=None, n_components=0.999, method='fastica', reject=dict(grad=4000e-13, mag=4e-12), random_state=42)
    
    #if subj_id=='sub003':
    #    parallel, run_func, _ = parallel_func(data.run_ica, n_jobs=10)
    #    parallel(run_func(tsss, method='fastica',n_components=0.99, reject=dict(grad=4000e-13, mag=4e-12), random_state=42) for tsss in (10, 1)) 
data.make_epochs(l_cutoff=l_cutoff, tsss=None, tmin=-0.2, tmax=2.9, reject_tmax=0.8, 
                    reject=dict(grad=400e-12, mag=3e-12),
                    decim=5, random_state = 42)
    
    #if subj_id == 'sub003':
    #     parallel, run_func, _ = parallel_func(data.make_epochs, n_jobs=10)
    #     parallel(run_func(tsss, tmin=-0.2, tmax=2.9, reject_tmax=0.8,baseline=(None, 0),reject=dict(grad=400e-12, mag=4e-12),decim=5, random_state = 42) for tsss in (10, 1)) # Maxwell filtered data
        
data.make_evoked(tsss=None)
    
    #if subj_id == 'sub003':
    #     parallel, run_func, _ = parallel_func(data.make_evoked, n_jobs=10)
    #     parallel(run_func(tsss) for tsss in (10, 1))
data.make_cov(tsss=None)
    
    #if subj_id == 'sub003':
    #    parallel, run_func, _ = parallel_func(data.make_cov, n_jobs=10)
    #    parallel(run_func(tsss) for tsss in (10, 1))
data.make_forward(tsss=None)
    
    #if subj_id == 'sub003':
    #    parallel, run_func, _ = parallel_func(data.make_forward, n_jobs=10)
    #    parallel(run_func(tsss) for tsss in (10, 1))
    
data.make_inverse(tsss=None)
    
    #if subj_id == 'sub003':
    #    parallel, run_func, _ = parallel_func(data.make_inverse, n_jobs=10)
    #    parallel(run_func(tsss) for tsss in (10, 1))
    
data.save_matlab()
