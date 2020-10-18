# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:55:21 2017

@author: Narayan Subramaniyam
MEGSET data preprocessing pipeline -> filtering, epoching, artifact correction with ica ....
"""
import os, mne, time, glob
import numpy as np
import os.path as op
from mne.minimum_norm import apply_inverse, read_inverse_operator, write_inverse_operator, make_inverse_operator
from mne.parallel import parallel_func
from mne import pick_types
from mne.preprocessing import create_ecg_epochs, create_eog_epochs, read_ica
from autoreject import get_rejection_threshold
from sklearn.model_selection import KFold
class MegSetPreProcessing:
    
    """
    This class contains all the necessary functions for
    megset data preprocessing
    
    useage : data = MegSetPreProcessing()
    
    """

    
    def __init__(self, subj_id, dir_path, res_path, event_id, total_runs, l_cutoff, h_cutoff, spacing):
      self.subj_id = subj_id
      self.dir_path = dir_path 
      self.res_path = res_path
      self.events_id = event_id
      self.total_runs = total_runs
      self.subjects_dir = op.join(self.dir_path, subj_id)
      self.l_cutoff=l_cutoff
      self.h_cutoff=h_cutoff
      self.spacing = spacing
      
      # path where maxfiltered meeg data is stored
      self.data_path = op.join(self.subjects_dir, 'MEG')
      
      # make directory to save processed stuff
      self.preprocessed_data = op.join(self.res_path, 'MEG', subj_id)
      for directory in [self.preprocessed_data]:
              if not op.isdir(directory):
                  os.makedirs(directory)
      self.processing_info = 'Preprocessing starts - May the force be with you: %s\n' %time.strftime('%d %b %Y, %H:%M:%S') + 48*'-' + '\n'
      
      
      
    def extract_events(self, stim_channel='STI101', min_duration=0.003):
        for run in range(1,self.total_runs+1):
             file_name = op.join(self.data_path,'run_%02d_raw.fif' %(run))
             raw = mne.io.read_raw_fif(file_name, verbose='INFO')
             mask = 4096 + 256  # mask for excluding high order bits
             events = mne.find_events(raw, stim_channel=stim_channel, consecutive='increasing', mask=mask,
                                 mask_type='not_and', min_duration=min_duration)
             #save event info
             dir_name = op.join(self.preprocessed_data, 'events')
             for directory in [dir_name]:
                 if not op.isdir(directory):
                     os.makedirs(directory)
                 target = op.join(self.preprocessed_data, 'events', 'events_run%02d-eve.fif' %(run))
                 mne.event.write_events(target, events)
        
        del events
      
    def run_filter (self, l_cutoff=None, h_cutoff=40, method='fir',fir_design='firwin'):
        
        """
        This function loads raw files. Note that, this function assumes that
        you are giving the maxfiltered data (obtained from neuromag software)
        as input. Thus, make sure to have maxfiltered your data before calling this func!
        
        """
        # load maxfiltered data and run bandpass and notch filters

        
        for run in range(1,self.total_runs+1):
            file_name = op.join(self.data_path,'run_%02d_sss.fif' %(run))
            raw = mne.io.read_raw_fif(file_name, preload=True, verbose=False)
            raw.set_channel_types({'EEG061':'eog', 
                                   'EEG062':'eog',
                                   'EEG063':'ecg',
                                   'EEG064': 'misc'})
            raw.rename_channels({'EEG061': 'EOG061',
                                 'EEG062': 'EOG062',
                                 'EEG063': 'ECG063'})
            
            #raw.notch_filter(np.arange(50, 451, 50), filter_length='auto',
                                 #phase='zero', fir_window='hamming')
           
            raw.filter(l_cutoff, h_cutoff, method=method, 
                           fir_design=fir_design, filter_length='auto', 
                           fir_window='hamming', phase='zero', 
                           l_trans_bandwidth='auto', 
                           h_trans_bandwidth='auto')
        
        # High-pass EOG to get reasonable thresholds in autoreject
            picks_eog = pick_types(raw.info, meg=False, eog=True)
            raw.filter(1.0, None, picks=picks_eog, 
                          method=method, fir_design=fir_design, 
                          filter_length='auto', fir_window='hann', 
                          l_trans_bandwidth='auto', 
                           h_trans_bandwidth='auto')   
            dir_name = op.join(self.preprocessed_data, 'filtered')
          
            for directory in [dir_name]:
                if not op.isdir(directory):
                  os.makedirs(directory)
            target = op.join(self.preprocessed_data, 'filtered', 
                           'filt_sss_run%02d_%s-%sHz.fif' 
                           %( run, str(l_cutoff), str(h_cutoff)))
            raw.save(target, overwrite=True)
            #self.make_notes("Raw data filtered from %d to %d Hz with %s method." %(l_cutoff, h_cutoff, method))
        
        del raw
    
    
    def run_ica(self, tsss=None, n_components=0.999, method='fastica', reject=dict(grad=4000e-13, mag=4e-12), random_state=42):

        raws=list()
        for run in range(1,self.total_runs+1):
            if tsss:
                   file_name = op.join(self.preprocessed_data,'tSSS', 'run_%02d_filt_tsss_%d_raw.fif' 
                                %(run, tsss))

            else:
                    file_name = op.join(self.preprocessed_data,'filtered', 'filt_sss_run%02d_%s-%sHz.fif' 
                                %(run, str(self.l_cutoff), str(self.h_cutoff)))
            raws.append(mne.io.read_raw_fif(file_name))
        raw = mne.concatenate_raws(raws)
        ica = mne.preprocessing.ICA(n_components, method=method, random_state=random_state)
        picks = mne.pick_types(raw.info, meg=True, eeg=False, eog=False, ecg=False, stim=False, exclude='bads')
        ica.fit(raw, picks=picks, decim=11, reject=reject)
        print('  Fit %d components (explaining at least %0.1f%% of the variance)'
          % (ica.n_components_, 100 * n_components))  
            #save ica stuff
            
        dir_name = op.join(self.preprocessed_data, 'ica')
        for directory in [dir_name]:
            if not op.isdir(directory):
                    os.makedirs(directory)
        if tsss:
            target = op.join(self.preprocessed_data, 'ica', 'run_concat-tsss_%d-ica.fif' %(tsss))
        else:
            target = op.join(self.preprocessed_data, 'ica', 'run_concat-ica.fif' )
        ica.save(target)
        
       

        
    def make_epochs(self, l_cutoff=None, tsss=None, tmin=-0.2, tmax=2.9, reject_tmax=0.8, 
                    reject=dict(grad=400e-12, mag=4e-12),
                    decim=5, random_state = 42):
        baseline = (None, 0) 
        events_list=list()
        raws_list=list()   
        for run in range(1,self.total_runs+1):
            if tsss:
                 file_name = op.join(self.preprocessed_data, 'tSSS', 'run_%02d_filt_tsss_%d_raw.fif' 
                                %(run, tsss))
            else:
                file_name = op.join(self.preprocessed_data,'filtered', 'filt_sss_run%02d_%s-%sHz.fif' 
                                %(run, str(self.l_cutoff), str(self.h_cutoff)))
            raw = mne.io.read_raw_fif(file_name, preload=True)
            events = mne.read_events(op.join(self.preprocessed_data, 'events', 
                                             'events_run%02d-eve.fif' % (run)))
            delay = int(round(0.0345 * raw.info['sfreq']))
            events[:, 0] = events[:, 0] + delay
            events_list.append(events)
            raws_list.append(raw)
            del raw
            del events
        raw, events = mne.concatenate_raws(raws_list, events_list=events_list)
        del raws_list
        
        picks = mne.pick_types(raw.info, meg=True, eog=True, eeg=True, 
                               stim=True, exclude=())
        
        epochs = mne.Epochs(raw, events, self.events_id, tmin, tmax,proj=True,
                        picks=picks, baseline=baseline, preload=False,
                        decim=decim, reject=None, reject_tmax=reject_tmax)
        #ICA
        if tsss:
            file_name = op.join(self.preprocessed_data, 'ica', 'run_concat-tsss_%d-ica.fif' %(tsss))
        else:
            file_name = op.join(self.preprocessed_data,'ica', 'run_concat-ica.fif' )
        ica = read_ica(file_name)
        ica.exclude=[]        
        ecg_epochs = create_ecg_epochs(raw, tmin=-.3, tmax=.3, preload=False)
        eog_epochs = create_eog_epochs(raw, tmin=-.5, tmax=.5, preload=False)
        del raw
        
        
        filter_label = '_highpass-%sHz' % self.l_cutoff
        n_max_ecg = 3  # use max 3 components
        ecg_epochs.decimate(5)
        ecg_epochs.load_data()
        ecg_epochs.apply_baseline((None, None))
        ecg_inds, scores_ecg = ica.find_bads_ecg(ecg_epochs, method='ctps',
                                             threshold=0.8)
        print('    Found %d ECG indices' % (len(ecg_inds),))
        ica.exclude.extend(ecg_inds[:n_max_ecg])
        ecg_epochs.average().save(op.join(self.preprocessed_data,'ica', '%s%s-ecg-ave.fif'
                                      % (self.subj_id, filter_label)))
        np.save(op.join(self.preprocessed_data,'ica', '%s%s-ecg-scores.npy'
                    % (self.subj_id, filter_label)), scores_ecg)
        
        del ecg_epochs

        n_max_eog = 3  # use max 2 components
        eog_epochs.decimate(5)
        eog_epochs.load_data()
        eog_epochs.apply_baseline((None, None))
        eog_inds, scores_eog = ica.find_bads_eog(eog_epochs)
        print('    Found %d EOG indices' % (len(eog_inds),))
        ica.exclude.extend(eog_inds[:n_max_eog])
        eog_epochs.average().save(op.join(self.preprocessed_data,'ica', '%s%s-eog-ave.fif'
                                      % (self.subj_id, filter_label)))
        np.save(op.join(self.preprocessed_data,'ica', '%s%s-eog-scores.npy'
                    % (self.subj_id, filter_label)), scores_eog)
        
        del eog_epochs
        
        if tsss:
            ica_out_name = op.join(self.preprocessed_data,'ica', 'run_concat-tsss_%d-ica.fif' 
                                %(tsss))
        else:
            ica_out_name = op.join(self.preprocessed_data,'ica', 'run_concat_highpass-%sHz-ica.fif' 
                                %(self.l_cutoff))
        ica.save(ica_out_name)
        

        epochs.load_data()
        ica.apply(epochs)

        print('  Getting rejection thresholds')
        reject = get_rejection_threshold(epochs.copy().crop(None, reject_tmax),
                                     random_state=random_state)
        epochs.drop_bad(reject=reject)
        print('The rejection dictionary is %s' % reject)
        print('  Dropped %0.1f%% of epochs' % (epochs.drop_log_stats(),))
        
        # save epochs
        dir_name = op.join(self.preprocessed_data, 'epochs')
        for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
        if tsss:
            target = op.join(self.preprocessed_data, 'epochs', '%s-tsss_%d-epo.fif' %(self.subj_id, tsss))
        else:
            target = op.join(self.preprocessed_data, 'epochs', '%s_highpass-%sHz-epo.fif' %(self.subj_id, self.l_cutoff))
        epochs.save(target)
        
        
    def make_evoked(self, tsss=None):
        if tsss:
            file_name = op.join(self.preprocessed_data,'epochs', '%s-tsss_%d-epo.fif' 
                                %(self.subj_id, tsss))
        else:
            file_name = op.join(self.preprocessed_data,'epochs', '%s_highpass-%sHz-epo.fif' 
                                %(self.subj_id, self.l_cutoff))
        epochs = mne.read_epochs(file_name, preload=True)
        evoked_famous = epochs['face/famous'].average()
        evoked_scrambled = epochs['scrambled'].average()
        evoked_unfamiliar = epochs['face/unfamiliar'].average()

    # Simplify comment
        evoked_famous.comment = 'famous'
        evoked_scrambled.comment = 'scrambled'
        evoked_unfamiliar.comment = 'unfamiliar'

        contrast = mne.combine_evoked([evoked_famous, evoked_unfamiliar,
                                   evoked_scrambled], weights=[0.5, 0.5, -1.])
        contrast.comment = 'contrast'
        faces = mne.combine_evoked([evoked_famous, evoked_unfamiliar], 'nave')
        faces.comment = 'faces'

    # let's make trial-count-normalized ones for group statistics
        epochs_eq = epochs.copy().equalize_event_counts(['face', 'scrambled'])[0]
        evoked_faces_eq = epochs_eq['face'].average()
        evoked_scrambled_eq = epochs_eq['scrambled'].average()
        assert evoked_faces_eq.nave == evoked_scrambled_eq.nave
        evoked_faces_eq.comment = 'faces_eq'
        evoked_scrambled_eq.comment = 'scrambled_eq'
        #save evoked
        dir_name = op.join(self.preprocessed_data, 'evoked')
        for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
        if tsss:
            target = op.join(self.preprocessed_data, 'evoked', '%s_tsss_%d-ave.fif' %(self.subj_id, tsss)) 
        else:
            target = op.join(self.preprocessed_data, 'evoked', '%s_highpass_%sHz-ave.fif' %(self.subj_id, self.l_cutoff)) 
        mne.evoked.write_evokeds(target, [evoked_famous, evoked_scrambled,
                                          evoked_unfamiliar, contrast, faces,
                                         evoked_faces_eq, evoked_scrambled_eq])
    
    def make_cov(self, tsss=None):
        if tsss:
            fname_epo = op.join(self.preprocessed_data, 'epochs', '%s-tsss_%d-epo.fif'
                            % (self.subj_id,tsss))
        else:
            fname_epo = op.join(self.preprocessed_data, 'epochs', '%s_highpass-%sHz-epo.fif'
                            % (self.subj_id,self.l_cutoff))
        dir_name = op.join(self.preprocessed_data, 'cov')
        for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
        if tsss:
            target = op.join(self.preprocessed_data, 'cov', '%s-tsss_%d.fif'
                            % (self.subj_id, tsss))
        else:
            target = op.join(self.preprocessed_data, 'cov', '%s_highpass-%sHz-cov.fif'
                            % (self.subj_id, self.l_cutoff))
        print('  Computing regularized covariance')
        epochs = mne.read_epochs(fname_epo, preload=True)
        cv = KFold(3, random_state=42)  # make sure cv is deterministic
        cov = mne.compute_covariance(epochs, tmax=0, method='shrunk', cv=cv)
        cov.save(target)
        
    def make_forward(self, tsss=None):
        
        dir_name = op.join(self.preprocessed_data, 'forward')
        for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
        if tsss:
            fname_ave = op.join(self.preprocessed_data, 'evoked', '%s_tsss_%d-ave.fif' %(self.subj_id, tsss))
        else:
            fname_ave = op.join(self.preprocessed_data, 'evoked', '%s_highpass_%sHz-ave.fif' %(self.subj_id, self.l_cutoff))
            
        fname_fwd = op.join(self.preprocessed_data, 'forward', '%s-meg-eeg-%s-fwd.fif'
                        % (self.subj_id, self.spacing))
        fname_trans = op.join(self.data_path, '%s-trans.fif' % self.subj_id)
        fname_src = op.join(self.res_path, 'subjects', self.subj_id, 'bem', '%s-%s-src.fif'
                        % (self.subj_id, self.spacing))
    # Here we only use 1-layer BEM because the 3-layer is unreliable
        fname_bem = op.join(self.res_path, 'subjects', self.subj_id, 'bem',
                        '%s-5120-bem-sol.fif' % self.subj_id)

        info = mne.io.read_info(fname_ave)
    # Because we use a 1-layer BEM, we do MEG only
        fwd = mne.make_forward_solution(info, fname_trans, fname_src, fname_bem,
                                    meg=True, eeg=False, mindist=5)
        mne.write_forward_solution(fname_fwd, fwd, overwrite=True)
      
    def make_inverse(self, tsss=None):

      dir_name = op.join(self.preprocessed_data, 'inverse')
      for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
      if tsss:
          fname_ave = op.join(self.preprocessed_data, 'evoked', '%s_tsss_%d-ave.fif' %(self.subj_id, tsss))
          fname_cov = op.join(self.preprocessed_data, 'cov', '%s-tsss_%d.fif'
                            % (self.subj_id, tsss))
      else:
          fname_ave = op.join(self.preprocessed_data, 'evoked', '%s_highpass_%sHz-ave.fif' %(self.subj_id, self.l_cutoff))
          fname_cov = op.join(self.preprocessed_data, 'cov', '%s_highpass-%sHz-cov.fif'
                            % (self.subj_id, self.l_cutoff))
      fname_fwd = op.join(self.preprocessed_data, 'forward', '%s-meg-eeg-%s-fwd.fif'
                        % (self.subj_id, self.spacing))
      fname_inv = op.join(self.preprocessed_data, 'inverse', '%s-meg-eeg-%s-inv.fif'
                        % (self.subj_id, self.spacing))

      evokeds = mne.read_evokeds(
        fname_ave, condition=['scrambled', 'unfamiliar', 'famous',
                              'faces', 'contrast',
                              'faces_eq', 'scrambled_eq'])
      cov = mne.read_cov(fname_cov)
      forward = mne.read_forward_solution(fname_fwd)

    # This will be an MEG-only inverse because the 3-layer BEMs are not
    # reliable, so our forward only has MEG channels.
    
      info = evokeds[0].info
      inverse_operator = make_inverse_operator(
        info, forward, cov, loose=0.2, depth=0.8)
      write_inverse_operator(fname_inv, inverse_operator)

    # Apply inverse
      snr = 3.0
      lambda2 = 1.0 / snr ** 2
      dir_name = op.join(self.preprocessed_data, 'stc_file')
      for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)

      for evoked in evokeds:
          stc = apply_inverse(evoked, inverse_operator, lambda2, "dSPM",
                            pick_ori='vector')
          stc.save(op.join(self.preprocessed_data, 'stc_file', 'mne_dSPM_inverse-%s'
                         % (evoked.comment)))
          
    def maxwell_filter(self, l_cutoff=None, h_cutoff=40, method='fir',fir_design='firwin'):
        ctc = op.join(self.data_path, 'ct_sparse.fif')
        cal = op.join(self.data_path, 'sss_cal.dat')
        N_JOBS=10
        dir_name = op.join(self.preprocessed_data, 'tSSS')
        for directory in [dir_name]:
            if not op.isdir(directory):
               os.makedirs(directory)
        sss_fname_out = op.join(self.preprocessed_data, 'tSSS', 'run_%02d_filt_tsss_%d_raw.fif')
        raw_fname_in = op.join(self.data_path,'run_%02d_raw.fif')
        sss_fname_in = op.join(self.data_path,'run_%02d_sss.fif')
        info = mne.io.read_info(sss_fname_in % 4)
        destination = info['dev_head_t']
    # Get the origin they used
        origin = info['proc_history'][0]['max_info']['sss_info']['origin']
        
        for run in range(1, 7):
            raw_in = raw_fname_in % run
            try:
                raw = mne.io.read_raw_fif(raw_in)
            except AttributeError:
            # Some files on openfmri are corrupted and cannot be read.
                warn('Could not read file %s. '
                 'Skipping run %s from subject %s.' % (raw_in, run, self.subj_id))
                continue
            print('  Run %s' % (run,))

            raw.set_channel_types({
            'EEG061': 'eog',
            'EEG062': 'eog',
            'EEG063': 'ecg',
            'EEG064': 'misc'
            })  # EEG064 free floating el.
            raw.rename_channels({
            'EEG061': 'EOG061',
            'EEG062': 'EOG062',
            'EEG063': 'ECG063'
            })
            raw.fix_mag_coil_types()

        # Read bad channels from the MaxFilter log.
            with open(
                op.join(self.data_path,
                        'run_%02d_sss_log.txt' % run)) as fid:
                for line in fid:
                    if line.startswith('Static bad channels'):
                        chs = line.split(':')[-1].split()
                        bads = ['MEG%04d' % int(ch) for ch in chs]
                        break
            raw.info['bads'] += bads
            raw_sss = mne.io.read_raw_fif(sss_fname_in % run, verbose='error')
            chpi_picks = mne.pick_types(raw_sss.info, meg=False, chpi=True)
            assert len(chpi_picks) == 9
            head_pos, t = raw_sss[chpi_picks]
            # Add first_samp.
            t = t + raw_sss.first_samp / raw_sss.info['sfreq']
            # The head positions in the FIF file are all zero for invalid positions
            # so let's remove them, and then concatenate our times.
            mask = (head_pos != 0).any(axis=0)
            head_pos = np.concatenate((t[np.newaxis], head_pos)).T[mask]
            # In this dataset due to old MaxFilter (2.2.10), data are uniformly
            # sampled at 1 Hz, so we save some processing time in
            # maxwell_filter by downsampling.
            skip = int(round(raw_sss.info['sfreq']))
            head_pos = head_pos[::skip]
            del raw_sss
       
            for st_duration in (10, 1):
                    print('    st_duration=%d' % (st_duration,))
                    raw_sss = mne.preprocessing.maxwell_filter(
                    raw, calibration=cal, cross_talk=ctc, st_duration=st_duration,
                    origin=origin, destination=destination, head_pos=head_pos)

                    raw_out = sss_fname_out % (run, st_duration)

            # Here we only low-pass MEG (assuming MaxFilter has high-passed
            # the data already), but we still need to band-pass EEG:
                    picks_meg = mne.pick_types(raw.info, meg=True, exclude=())
                    raw_sss.filter(
                            None, 40, l_trans_bandwidth='auto', h_trans_bandwidth='auto',
                            filter_length='auto', phase='zero', fir_window='hamming',
                            fir_design='firwin', n_jobs=N_JOBS, picks=picks_meg)
                    picks_eeg = mne.pick_types(raw.info, eeg=True, exclude=())
                    raw_sss.filter(
                            l_cutoff, 40, l_trans_bandwidth='auto', h_trans_bandwidth='auto',
                            filter_length='auto', phase='zero', fir_window='hamming',
                            fir_design='firwin', n_jobs=N_JOBS, picks=picks_eeg)
            # High-pass EOG to get reasonable thresholds in autoreject
                    picks_eog = mne.pick_types(raw.info, meg=False, eog=True)
                    raw_sss.filter(
                            l_cutoff, None, picks=picks_eog, l_trans_bandwidth='auto',
                            filter_length='auto', phase='zero', fir_window='hann',
                            fir_design='firwin')
                    raw_sss.save(raw_out, overwrite=True)
                    
    def save_matlab(self):
        import scipy.io as io
        from mne import convert_forward_solution
        fname_epo = op.join(self.preprocessed_data,'epochs', '%s_highpass-%sHz-epo.fif' 
                                %(self.subj_id, self.l_cutoff))
        fname_evo = op.join(self.preprocessed_data, 'evoked', '%s_highpass_%sHz-ave.fif' %(self.subj_id, self.l_cutoff))
        fname_cov = op.join(self.preprocessed_data, 'cov', '%s_highpass-%sHz-cov.fif'
                            % (self.subj_id, self.l_cutoff))
        fname_fwd = op.join(self.preprocessed_data, 'forward', '%s-meg-eeg-%s-fwd.fif'
                        % (self.subj_id, self.spacing))
        fname_inv = op.join(self.preprocessed_data, 'inverse', '%s-meg-eeg-%s-inv.fif'
                        % (self.subj_id, self.spacing))
        fname_src = op.join(self.res_path, 'subjects', self.subj_id, 'bem', '%s-%s-src.fif'
                        % (self.subj_id, self.spacing))

        evoked = mne.read_evokeds(fname_evo)
        epochs = mne.read_epochs(fname_epo)
        noise_cov = mne.read_cov(fname_cov)
        patch_stats = True  # include high resolution source space
        src_full = mne.read_source_spaces(fname_src, patch_stats=patch_stats)
        fwd = mne.read_forward_solution(fname_fwd)
        inv = read_inverse_operator(fname_inv)
        labels = mne.read_labels_from_annot(self.subj_id, 'aparc', 'both', subjects_dir= op.join(self.res_path, 'subjects'))
        anat_parcels=[]
        
        for i in range(68):
            anat_parcels.append(labels[i].vertices)

        fwd_n = convert_forward_solution(fwd, force_fixed=True, surf_ori=True, use_cps=True)
        famous_evo, scrambled_evo, unfamiliar_evo, contrast_evo, faces_evo = evoked[:5]
        
        
        dir_name = op.join(self.preprocessed_data, 'matlab')
          
        for directory in [dir_name]:
            if not op.isdir(directory):
                os.makedirs(directory)
        path_sub = dir_name + '/'
        src = inv['src']
        io.savemat(path_sub + 'src_left_%s.mat'%(self.subj_id), {'src':src[0]} )
        io.savemat(path_sub + 'src_right_%s.mat'%(self.subj_id), {'src':src[1]} )
        io.savemat(path_sub + 't_%s.mat'%(self.subj_id), {'t':scrambled_evo.times} )
        io.savemat(path_sub + 'XI_%s.mat'%(self.subj_id), {'XI':noise_cov['data']} )
        io.savemat(path_sub + 'src_full_left_%s.mat'%(self.subj_id), {'src_full':src_full[0]} )
        io.savemat(path_sub + 'src_full_right_%s.mat'%(self.subj_id), {'src_full':src_full[1]})
        io.savemat(path_sub + 'G_%s.mat'%(self.subj_id), {'G':fwd_n['sol']['data']} )
        io.savemat(path_sub + 'anat_labels_%s.mat'%(self.subj_id), {'anat_labels':anat_parcels})
        io.savemat(path_sub + 'Y_trials_famous_%s.mat'%(self.subj_id), {'Y':epochs['face/famous'].pick_types(meg=True).apply_baseline().get_data()} )
        io.savemat(path_sub + 'Y_trials_unfamiliar_%s.mat'%(self.subj_id), {'Y':epochs['face/unfamiliar'].pick_types(meg=True).apply_baseline().get_data()} )
        io.savemat(path_sub + 'Y_trials_scrambled_%s.mat'%(self.subj_id), {'Y':epochs['scrambled'].pick_types(meg=True).apply_baseline().get_data()} )


        
        
        

           
        

        
    
    
