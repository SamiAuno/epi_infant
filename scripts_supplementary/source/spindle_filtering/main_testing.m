
% Testing and debugging scipt for spindle extraction


% sFiles = {'/media/sami/ADATA SD700/EPI_INFNT_BS/brainstorm_db/EPI_Infant/data/EPI_02_1v/@rawEPI_2_Egal__Ayanle_2017-06-08_11-36-20_band_bl_resample/data_0raw_EPI_2_Egal__Ayanle_2017-06-08_11-36-20_band_bl_resample.mat'};
sFiles = {'/media/sami/ADATA SD700/EPI_INFNT_BS/brainstorm_db/EPI_Infant/data/EPI_09_1v/@rawEPI_9_LUKAS_2018-05-17_09-46-07_band_03_bl_resample_interpbad_02/data_0raw_EPI_9_LUKAS_2018-05-17_09-46-07_band_03_bl_resample_interpbad_02.mat'};
% sFiles = {'/media/sami/ADATA SD700/EPI_INFNT_BS/brainstorm_db/EPI_Infant/data/EPI_09_1v/@rawEPI_9_LUKAS_2018-05-17_09-46-07/data_0raw_EPI 9_LUKAS_2018-05-17_09-46-07.mat'};
csvName = 'Analysis_Epochs.csv';
Epochs_table = readEpochTable(csvName);

eeg = rerefBIP(sFiles{1},Epochs_table);

% process_spindles(sFiles);
%% These are with the new
%  data = [channels, time]
%  hdr  = structure that MUST HAVE ...
%       hdr.info.sfreq      = sampling frequency [Hz]
%       hdr.info.ch_names   = cell of channel names.

data = eeg.signal;
hdr.info.sfreq = eeg.srate;
hdr.info.ch_names = eeg.chName;
options.MinPeakProminence = 1e-3;
tic
spindle_prob = LSM_spindle_probabilities(data, hdr,options);
toc
%%
spindle_det = LSM_spindle_detections(spindle_prob);


% Visualize
channel = eeg.chName{12};
figure(2);
clf
LSM_spindle_visualizer(data, hdr, spindle_det, channel)



%% These are with the old spindle detection script
% 
% [eeg, OUTPUT] = extract_spindles(eeg);
% 
% %%
% 
OUTPUT_new = removeLowPowerSpindles(OUTPUT,eeg);
% 
% %%
% select_spindles(eeg,OUTPUT_new);

%%




