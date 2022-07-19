% 25/10/2019
% EPI_Infant project
%
% For spindle filtering, preprocess a file until spindle filtering, then
% delete the intermittent files. The resulting file can be used to test
% spindle filtering.

%% NOTE THE DIFFERENT BROAD-BAND FILTER THAN IN THE ACTUAL PROCESS
%% THIS IS ONLY FOR VISUALIZATION

%% ======= SETTINGS ======= %%

% Path to Source Files
sFiles = {'/media/sami/ADATA SD700/EPI_INFNT_BS/brainstorm_db/EPI_Infant/data/EPI_09_1v/@rawEPI_9_LUKAS_2018-05-17_09-46-07/data_0raw_EPI 9_LUKAS_2018-05-17_09-46-07.mat'};


% DEFINE STANDARD CHANNEL LABEL LIST
ch_standard = {'Fp1';'Fpz';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';'FC2';'FC6';'M1';'T7';'C3';'Cz';'C4';'T8';'M2';'CP5';'CP1';'CP2';'CP6';'P7';'P3';'Pz';'P4';'P8';'POz';'O1';'O2';'EOG';'AF7';'AF3';'AF4';'AF8';'F5';'F1';'F2';'F6';'FC3';'FCz';'FC4';'C5';'C1';'C2';'C6';'CP3';'CP4';'P5';'P1';'P2';'P6';'PO5';'PO3';'PO4';'PO6';'FT7';'FT8';'TP7';'TP8';'PO7';'PO8';'Oz'};

% Options
Option.samplingFrequency = 250;

path_inverseSolution = {'/media/sami/ADATA SD700/EPI_INFNT_BS/Models_for_source_analysis/201919/results_dSPM_EEG_KERNEL_191129_1541.mat'};
path_globalChannelFile = {'/media/sami/ADATA SD700/EPI_INFNT_BS/Models_for_source_analysis/201919/channelFile.txt'};


%% ======= STEP 1 ======= %%
% Check that the number of channels match


% Start a new report
bst_report('Start', sFiles);

% Go through all studies and get the channel label list for each
Nfiles = length(sFiles);
channelLists = cell(Nfiles,1);

% Get the channel list of each measurement
for si = 1:Nfiles
    [sStudy, iStudy] = bst_get('AnyFile', sFiles{si});
    ChannelFile = bst_get('ChannelFileForStudy', sFiles{si});
%     ChannelFile = sStudy.Channel(1).FileName;
    ChannelMat = in_bst_channel(ChannelFile);
    channelLists{si} = {ChannelMat.Channel.Name}';
end

% Check that each channel list contains the standard channel labels
trueNumberOfChannels = true(Nfiles,1);
extraChannelLabels = [];
for ci = 1:Nfiles
    
    missingChannelLabels = find(~ismember(ch_standard,channelLists{ci}));
    
    extraChannelLabels = [extraChannelLabels;channelLists{ci}(~ismember(lower(channelLists{ci}),lower(ch_standard)))];
    
    if ~isempty(missingChannelLabels)
        trueNumberOfChannels(ci) = false;
    end
end

% Get unique extra channels and reorganize to single string so that it can
% be fed to process_channel_setbad function.
extraChannelLabels = unique(extraChannelLabels);
extraChannelString = strjoin(extraChannelLabels,', ');

% Remove those measurements that lack the correct number of channels from
% further analysis and save the names of the excluded files
excluded_measurements = sFiles(~trueNumberOfChannels);
sFiles = sFiles(trueNumberOfChannels);


% Channels, that do not appear in the standard channel label list are set
% to Misc-type. This ensures that they are not included in later analysis.
sFiles = bst_process('CallProcess', 'process_channel_settype', sFiles, [], ...
    'sensortypes', extraChannelString, ...
    'newtype',     'Misc');

% In case you want to set channels as BAD, you can use this function.
% sFiles = bst_process('CallProcess', 'process_channel_setbad', sFiles, [], ...
%     'sensortypes', extraChannelString);


%% ======= STEP 3 ======= %%
% Remove 50 Hz notch

% sFiles_notch = bst_process('CallProcess', 'process_notch', sFiles_bl, [], ...
%     'sensortypes', 'EEG', ...
%     'freqlist',    50, ...      % Notch filter frequency, Hz
%     'cutoffW',     1, ...       % -3 dB Notch bandwidth, Hz
%     'useold',      0, ...
%     'read_all',    0);


%% ======= STEP 4 ======= %%
% Broad-band band-pass filter
% NOTE THE DIFFERENT BROAD-BAND FILTER THAN IN THE ACTUAL PROCESS
% THIS IS ONLY FOR VISUALIZATION
sFiles_band = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
    'sensortypes', 'EEG', ...
    'highpass',    1, ...
    'lowpass',     40, ...
    'tranband',    0, ...
    'attenuation', 'relax', ...  % 'relax' = 40dB, 'strict' = 60dB
    'ver',         '2019', ...  % 2019
    'mirror',      0, ...
    'read_all',    0);

%% ======= STEP 2 ======= %%
% Remove DC

sFiles_bl = bst_process('CallProcess', 'process_baseline', sFiles_band, [], ...
    'baseline',    [], ...
    'sensortypes', 'EEG', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'read_all',    0);

%% ======= STEP 5 ======= %%
% Resample data

sFiles_resample = bst_process('CallProcess', 'process_resample', sFiles_bl, [], ...
    'freq',     Option.samplingFrequency, ...
    'read_all', 0);


%% ======= STEP N ======= %%
% Interpolate bad channels
sFiles_interp = bst_process('CallProcess', 'process_interpolateBadChan', sFiles_resample, [], ...
    'invfile',     {path_inverseSolution{1}, 'ARRAY-TIMES'}, ...
    'chanfile',    {path_globalChannelFile{1}, 'ARRAY-TIMES'}, ...
    'sensortypes', 'EEG');



%% 
bst_process('CallProcess', 'process_delete',sFiles_bl,[],'target',2);
% bst_process('CallProcess', 'process_delete',sFiles_notch,[],'target',2);
bst_process('CallProcess', 'process_delete',sFiles_band,[],'target',2);
bst_process('CallProcess', 'process_delete',sFiles_resample,[],'target',2);


