close all
clear all
addpath('/l/fieldtrip-20190403')
ft_defaults
addpath('/m/nbe/scratch/braintrack/ecog/SubjectNY394')
load epoch_data_clean_chan
fid = fopen('NY394_MRI_coor.txt');
elec_info = textscan(fid,'%s %f %f %f %s');
fclose(fid);

% create FieldTrip electrode structure
% elec       = [];
elec.label   = elec_info{1};
elec.elecpos = [elec_info{2} elec_info{3} elec_info{4}];
elec.unit    = 'mm';

%% load pial surface
load('NY394_MRI_rh_pial_surface.mat');

% create FieldTrip surface mesh structure
mesh1      = [];
mesh1.pos  = surface.pos;
mesh1.tri  = surface.tri;
mesh1.unit = 'mm';

figure
%% plot surface and electrodes
ft_plot_mesh(mesh1, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none')
view([90 25])
lighting gouraud
material shiny
camlight

% plot electrodes
hs = ft_plot_sens(elec, 'style', 'ko', 'label', 'on');

cfg = [];
cfg.resamplefs  = 128;
cfg.demean      = 'no';
cfg.detrend     = 'no';
resampled = ft_resampledata(cfg, epoch_data_clean_chan);


% calculate ERPs
cfg                  = [];
cfg.keeptrials       = 'yes'; % keep trials for statistics
cfg.preproc.lpfilter = 'yes';
cfg.preproc.lpfreq   = 40;    % smooth ERP with low-pass filter
cfg.preproc.hpfilter = 'yes';
cfg.preproc.hpfreq   = 0.1;     % reduce slow drifts
cfg.preproc.detrend  = 'yes';
cfg.polyremoval = 'yes';
cfg.preproc.hpfiltdir =  'onepass-zerophase';
cfg.preproc.hpfilttype =  'fir';
cfg.preproc.lpfiltdir =  'onepass-zerophase';
cfg.preproc.lpfilttype =  'fir';
cfg.preproc.reref = 'yes';
cfg.preproc.refchannel = 'all';
cfg.preproc.lpfiltord = 30;
cfg.preproc.hpfiltord = 30;
epoch_data_clean_chan=resampled;
cfg.trials = find(epoch_data_clean_chan.trialinfo == 3); % select only 'object' trials (event code 3)
ERP_object = ft_timelockanalysis(cfg, epoch_data_clean_chan);

cfg.trials = find(epoch_data_clean_chan.trialinfo == 7); % select only 'face' trials (event code 7)
ERP_face   = ft_timelockanalysis(cfg, epoch_data_clean_chan);



% baseline correction
cfg          = [];
cfg.baseline = [-.5 -.01];

ERP_object_bl = ft_timelockbaseline(cfg,ERP_object);
ERP_face_bl   = ft_timelockbaseline(cfg,ERP_face);

erp_ecog=ERP_face_bl.trial;
erp_ecog=permute(erp_ecog, [3 2 1]);
t0 = 0;
t0_ind = find(abs(ERP_face_bl.time-t0) == min(abs(ERP_face_bl.time-t0)));
tend_ind  = t0_ind + (0.5*128); %(500 ms after 0 ms, with 128 Hz as Fs)
chan = [ 76 72];
erp_ecog_seg = erp_ecog(t0_ind:tend_ind,chan,:); %0 to 500 ms

figure
plot(epoch_data_clean_chan.time{1,1}(t0_ind:tend_ind)*1000, squeeze(mean(erp_ecog_seg(:,:,:),3)), 'LineWidth', 3)


