%% Compute connectivity
% University of Helsinki, BABA Center
% Sami Auno 2019


clear
%% ---------- SET THESE ---------- %%

%% Save folder
% Where the connectivity matrices will be saved
path_folder_ConnectivitySave = '/path_to_savefolder/Connectivity_results';

%% Paths to Head model and Inverse model files

path_headModel      = '/modeling_data/headmodel.mat'; % Headmodel generated in Brainstorm
path_inverseModel   = '/modeling_data/dSPM_EEG_KERNEL.mat'; % Inverse model generated in Brainstorm
%% Epoch selection options
path_EpochTable = '/scripting_data/Analysis_Epochs.csv'; % File of epochs to be included in analysis

OPTIONS.epochLength = 180;              % The length of the epochs to be analysed in seconds.
OPTIONS.epochSegmentLength = 30;        % If there is more data than OPTIONS.epochLength per patient per sleep state, the data is divided into segments of length OPTIONS.epochSegmentLength that are equidistant from each other so that their total sum is OPTIONS.epochLength

%% Path to saved spindels. Used to cut spindels out from epochs.
spindelPath = '/path_to_spindle_matrix/Spindles_2subj.mat';

%% Filter Options

OPTIONS.NewFiltering = 1;               % Checks if filter bank filtered data already exists, and if it does: 0 = don't filter again but use the saved, 1 = filter again.
OPTIONS.NewSources = 1;                 % Checks if collapsed source data already exists, and if it does: 0 = don't calculate again but use the saved, 1 = calculate again.
OPTIONS.isRelax = true;
OPTIONS.TranBand = 0;
OPTIONS.numberOfBandwidths = 25;
OPTIONS.step_Bandwidth = 1.2;
OPTIONS.Save2Brainstorm = 0;            % Save steps to brainstorm. Not recommended.
OPTIONS.cutSpindlesOut = 0;             % 1 = yes, 0 = no. If yes, is performed if spindles calculated and saved in the previous processing step.


%% Paths to source files
paths_sFiles = {...    
    'path_to_cleaned_EEG_subject1.mat', ...
    'path_to_cleaned_EEG_subject2.mat', ...
    'path_to_cleaned_EEG_subject3.mat'};
%% ---------- SCRIPT STARTS HERE ---------- %%

% Change the paths to full path if relative.
paths_sFiles = file_fullpath(paths_sFiles);

N_sFiles = length(paths_sFiles);

%% Generate bandwidth bank
% save to options.
[OPTIONS.cutoffRange,OPTIONS.stopRange] =  generate_filter_bank(OPTIONS.numberOfBandwidths, OPTIONS.step_Bandwidth);

%% Read Analysis_Epochs.csv
% This contains the patient name, sleep stage (type) of the epoch, and the time stamps of the epoch that you
% want to analyse. This file must locate in the same folder as this script.
currentScriptPath = matlab.desktop.editor.getActiveFilename;
[currentFolderPath,~,~] = fileparts(currentScriptPath);

Epochs_table = readEpochTable(path_EpochTable);
% Epochs_table = readtable(fullfile(currentFolderPath,'Analysis_Epochs.csv'));
Epochs = table2cell(Epochs_table);

size_Epochs_table = size(Epochs_table);

%% Read spindels-data if available
if isfile(spindelPath) && logical(OPTIONS.cutSpindlesOut)
    tmp_spindle_struct = load(spindelPath);
    OPTIONS.spindles = tmp_spindle_struct.spindles;
else
    OPTIONS.spindles = [];
end

%% Load data
HeadModel = load(path_headModel);
kernel = load(path_inverseModel);
InverseModel = kernel.ImagingKernel;
channelIndecis = kernel.GoodChannel;
ForwardModel = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);

%% Make save folder
% Check that save folder exists
if ~isfolder(path_folder_ConnectivitySave)
    disp('path_ConnectivitySave -folder does not exist.');
    disp('Please correct the path path_ConnectivitySave to correct path.')
    disp('Terminating');
    return
end

% Make new time-stamped folder for this analysis run
path_folder_currentRun = fullfile(path_folder_ConnectivitySave,['Connectivity_run_',datestr(now,'yymmdd_HHMM')]);
mkdir(path_folder_currentRun);

%% Read source files one file at a time.
for s = 1:length(paths_sFiles)
    
    %% Get the patient name
    patient_name = getPatientName(paths_sFiles{s});
    if isempty(patient_name)
        continue
    end
    %% Check that this patient has corresponding epochs
    % If not, then pass to next patient
    if ~any(strcmp(Epochs_table.PatientName,patient_name))
        disp(['Patient ',patient_name, ' does not have epochs.']);
        continue
    else
        disp(['Processing patient ',patient_name,' (', num2str(s),'/',num2str(N_sFiles),')']);
    end
    
    
    %% Folder name
    sFolder = bst_fileparts(paths_sFiles{s});
    
    %% Read data
    [sMatrix, ~] = in_bst(paths_sFiles{s});
    
    sMatrix.sFreq = length(sMatrix.Time)/(sMatrix.Time(end)-sMatrix.Time(1));   % Sampling frequency
    sMatrix.F = sMatrix.F(channelIndecis,:);    % Take the appropriate channels
    %% Load the surface file
    SurfaceMat = in_tess_bst(kernel.SurfaceFile);

    iAtlas = SurfaceMat.iAtlas;     % NOTE! This assumes that the correct parcellation has been selected! Check
    Parcellation = SurfaceMat.Atlas(iAtlas).Scouts;
    
    N_parcels = length(Parcellation);
    
    %% Get epochs for this patient
    
    [epochs_index_mask,epochs_timeStamp,Epoch_type_array] = getEpochIndecis(Epochs_table,patient_name, sMatrix.Time,sMatrix.sFreq,OPTIONS);
    
    
    %% ---------- PROCESSING STARTS HERE ---------- %%
    
    % Make save folders for intermediate data, eg. filtered data and
    % collapsed source space data
    path_folder_intermediateData = fullfile(path_folder_ConnectivitySave,'intermediate_data');
    if ~isfolder(path_folder_intermediateData)
        disp('intermediate_data -folder does not exist.');
        disp('Making new folder for intermediate data.');
        mkdir(path_folder_intermediateData);
    end
    
    uniqStates = unique(Epoch_type_array);
    N_UniqEpochs = length(uniqStates);              % Number of unique epochs
    
    % Process one epoch at a time
    for epoch_i = 1:N_UniqEpochs
        
        disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'): Processing  epoch ', num2str(epoch_i),'/', num2str(N_UniqEpochs), '(',uniqStates{epoch_i},')']);
        
        %% Epoch time stamp
        % The time stamp is the in format HHMMSS-HHMMSS so that the first
        % time string before dash is the begining time of the first epoch
        % segment of the current state, and the last time string is the end
        % time of the last epoch segment of the current state.
        epoch_segment_indices = find(contains(Epoch_type_array,uniqStates{epoch_i}));
        epoch_timestamp = [datestr(epochs_timeStamp(epoch_segment_indices(1),1),'HHMMSS'),'-',datestr(epochs_timeStamp(epoch_segment_indices(end),2),'HHMMSS')];
        
        %% Get data for this epoch
        epoch_sMatrix = sMatrix;
        
        % Add 5 seconds to the begining and end of the measurement. There
        % are removed after filtering to avoid possible edge effects. Just
        % in case, you know.
        edgeSamples = round(5*sMatrix.sFreq);
        
        % Find the first and last index of the epoch. Used to add edge
        % samples
        first_last_index_epoch = [find(epochs_index_mask(epoch_i,:),1,'first') , find(epochs_index_mask(epoch_i,:),1,'last')];
        
        % Get the epochs
        epoch_sMatrix.F = sMatrix.F(:,epochs_index_mask(epoch_i,:));
        epoch_sMatrix.Time = sMatrix.Time(1,epochs_index_mask(epoch_i,:));
        
        
        % Add the edge samples
        epoch_sMatrix.F = [sMatrix.F(:,first_last_index_epoch(1) - 1 - edgeSamples : first_last_index_epoch(1) - 1 ) , epoch_sMatrix.F , sMatrix.F(:, first_last_index_epoch(2) + 1 : first_last_index_epoch(2) + 1 + edgeSamples ) ];
        epoch_sMatrix.Time = [sMatrix.Time(1,first_last_index_epoch(1) - 1 - edgeSamples : first_last_index_epoch(1) - 1 ) , epoch_sMatrix.Time , sMatrix.Time(1, first_last_index_epoch(2) + 1 : first_last_index_epoch(2) + 1 + edgeSamples ) ];
        
        %% Bandpass filter the file by filter bank.
        % produces [N_channels, N_samples, N_filters] sized matrix
        
        % Name of the saved filtered data file
        filename_filtered_data = ['bw_filtered_data_',patient_name,'_',epoch_timestamp,'.mat'];
        path_filtered_data = fullfile(path_folder_intermediateData,filename_filtered_data);
        
        % Calculate filter bank if:
        % 'bw_filtered_data_[Patient]_[Epochtimestamp].mat' does not exists
        % OR OPTIONS.NewFiltering is true

        % else load 'bw_filtered_data_[Patient]_[Epochtimestamp].mat'
        
        if ~isfile(path_filtered_data) || OPTIONS.NewFiltering
            
            disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Bandpass filtering.']);
            
            % Filter data
            bw_filtered_data = bandpass_filter(epoch_sMatrix, OPTIONS);
            
            % Remove the edge samples
            bw_filtered_data = bw_filtered_data(:,1 + edgeSamples : end - edgeSamples, :);
            
            % Save filtered data to file
            save(path_filtered_data,'bw_filtered_data');
            
            disp('Data bandpass filtered and saved.');
            
        else
            disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Loading pre-filtered data.']);
            % Load precalculated filtered data
            load(path_filtered_data);
            disp('Data loaded.');

        end
        
        % Remove edge samples from epoch_sMatrix just to be consistent.
        epoch_sMatrix.F     = epoch_sMatrix.F(:,1 + edgeSamples : end - edgeSamples);
        epoch_sMatrix.Time  = epoch_sMatrix.Time(1,1 + edgeSamples : end - edgeSamples);
        
        % Calculate number of samples
        N_samples = size(epoch_sMatrix.F,2);
        
        %% Calculate source time series collapsed to parcels.
        % produces [N_parcels, N_samples, N_filters] sized matrix
        
        % Name of the saved filtered data file
        filename_parcelTimeSeries = ['parcelTimeSeries_',patient_name,'_',epoch_timestamp,'.mat'];
        path_parcelTimeSeries = fullfile(path_folder_intermediateData,filename_parcelTimeSeries);
        
        % Calculate sources if:
        % 'parcelTimeSeries_[Patient]_[Epochtimestamp].mat' does not exists
        % OR OPTIONS.NewSources is true

        % else load 'parcelTimeSeries_[Patient]_[Epochtimestamp].mat'
        
        if ~isfile(path_parcelTimeSeries) || OPTIONS.NewSources
            
            disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Calculating source projection.']);
            
            % Initialise matrix
            parcelTimeSeries = zeros(N_parcels,N_samples,OPTIONS.numberOfBandwidths);

            
            % Loop through each bandwidth
            for bw_i = 1:OPTIONS.numberOfBandwidths
                parcelTimeSeries(:,:,bw_i) = process_collapse2parcels(paths_sFiles{s},bw_filtered_data(:,:,bw_i), path_headModel,path_inverseModel,OPTIONS.Save2Brainstorm);
            end
            
            % Save to file
            save(path_parcelTimeSeries,'parcelTimeSeries');
            
            disp('Sources calculated and collapsed to parcels.');
            
        else
            
            disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Loading pre-calculated source projections.']);
            
            % Load precalculated from file
            load(path_parcelTimeSeries);
            
            disp('Data loaded.');

        end
        
        
        %% Calculate connectivity metrics

        disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Calculating connectivities.']);

        % Initialise connectivity metric matrices
        connectivity.CC      = zeros(N_parcels,N_parcels,OPTIONS.numberOfBandwidths);
        connectivity.oCC     = zeros(N_parcels,N_parcels,OPTIONS.numberOfBandwidths);
        connectivity.dbwPLI  = zeros(N_parcels,N_parcels,OPTIONS.numberOfBandwidths);
        connectivity.PLV     = zeros(N_parcels,N_parcels,OPTIONS.numberOfBandwidths);
        
        % Connectivities are here calculate one bandwidth at a time
        for bw_i = 1:OPTIONS.numberOfBandwidths
            
            %% (1) Calculate the analytic signal of the parcel time series
                
            % hilbert() applies the tranforms along the first dimension, but we
            % need it along the second: transpose first and then hilbert()
            parcelTimeSeries_Hilbert = hilbert( (parcelTimeSeries(:,:,bw_i)).' );
            
            
            %% (2) Calculate Phase Locking Value
            
            tmp = parcelTimeSeries_Hilbert.'; % Transpose
            connectivity.PLV(:,:,bw_i) = abs( (tmp./abs(tmp)) * (tmp./abs(tmp))' ) / N_samples;
            
            %% (3) Calculate debiased weighted Phase Lag Index
            
            % Debiased -> true
            debias = true;
            
            % Check that fidelity operator exists in the HeadModel
            % structure. If not, make empty fidelity operator.
            if ~isfield(HeadModel,'FidelityOperator')
                HeadModel.FidelityOperator = [];
            end
            
            % Calculate debiased weighted Phase Lag Index
            connectivity.dbwPLI(:,:,bw_i) = calculate_wPLI(parcelTimeSeries_Hilbert, HeadModel.FidelityOperator , debias);
            
            %% (4) Calculate Correlation Coefficient
            
            CC_temp = corr(abs(parcelTimeSeries_Hilbert), abs(parcelTimeSeries_Hilbert));
            
            % Fidelity correction
            if isempty(HeadModel.FidelityOperator)
                connectivity.CC(:,:,bw_i) = CC_temp;
            else
                connectivity.CC(:,:,bw_i) = CC_temp .* HeadModel.FidelityOperator;
            end
            
            %% (5) Calculate orthogonal Correlation Coefficient
            
            WindowLength = round(sMatrix.sFreq);  % Use one second window length
            
            connectivity.oCC(:,:,bw_i) = get_oCC(parcelTimeSeries_Hilbert,WindowLength, HeadModel.FidelityOperator);
            
            
        end     % Loop, bandwidths
        
        % Save connectivites to path_folder_currentRun
        filename_connectivity = ['connectivity_',patient_name,'_',epoch_timestamp,'_',uniqStates{epoch_i},'.mat'];
        path_connectivity = fullfile(path_folder_currentRun,filename_connectivity);
        
        save(path_connectivity,'-struct','connectivity');
        
        disp([patient_name,' (', num2str(s),'/',num2str(N_sFiles),'), epoch (', num2str(epoch_i),'/', num2str(N_UniqEpochs),'): ', 'Connectivities calculated and saved.']);
        
        
    end     % Loop, epochs
end     % Loop, source files / patients










%% Band pass filter by filterbank

function filtered_data = bandpass_filter(sMatrix, OPTIONS)

% Initialise matrix
filtered_data = zeros(size(sMatrix.F,1),size(sMatrix.F,2),OPTIONS.numberOfBandwidths);

% Loop through each bandwidth
for i = 1:OPTIONS.numberOfBandwidths
    [filtered_data(:,:,i), ~, ~] = bst_bandpass_hfilter(sMatrix.F, sMatrix.sFreq, OPTIONS.cutoffRange(i,1), OPTIONS.cutoffRange(i,2), 0, OPTIONS.isRelax, [], OPTIONS.TranBand, []);
end

end



%% Weighted Phase Lag Index

% A wrapper function
function wpli = calculate_wPLI(A, fidelity , debias)
    % Assumes that data in A:
    % (1) is analytical
    % (2) has channels/parcels in dimension 2
    % (3) has samples in dimension 1
    
    N_channels = size(A,2);
    
    wpli = zeros(N_channels,N_channels);
    
    % Loop over channels/parcels
    for ch_x = 1:N_channels
        
        X = A(:,ch_x);
        
        % Loop again over channels
        for ch_y = ch_x:N_channels
            
            % Check if fidelity matrix is given or not, and if it is, that
            % the value between these channels is 1
            if isempty(fidelity) || fidelity(ch_x,ch_y)
                Y = A(:,ch_y);
                wpli_xy = get_wPLI(X, Y, debias);
                wpli(ch_x,ch_y) = wpli_xy;
            end
        end
        
        
    end
    
    
    
end

% wPLI calculation
function wpli = get_wPLI(X, Y, debias)

    % init

    X = real(X(:));
    Y = real(Y(:));
    L = length(X);

    % cross-spectral density
    Pxy = cpsd(X, Y, L, 1, [], []);

    % compute wpli
    Pxy = imag(Pxy);         % make everything imaginary  

    outsum = nansum(Pxy, 1); % compute the sum; this is 1 x size(2:end)
    outsumW = nansum(abs(Pxy), 1); % normalization of the WPLI

    if debias == 1
        outssq = nansum(Pxy .^ 2, 1);
        wpli = (outsum .^ 2 - outssq) ./ (outsumW .^ 2 - outssq); % do the pairwise 
    else
        wpli = outsum ./ outsumW; % estimator of E(Im(X))/E(|Im(X)|)
    end    

end % end


%% Orthogonal Correlation Coefficient


% 15/02/2019, Anton Tokariev
% Modified by Sami Auno, 11/12/2019
% Filters non-zero time series of the 3D array along the 3rd dimension
function oCC = get_oCC(A,WindowLength, FidelityOperator)

% A_orth - envelopes of mutually orthogonolized signals [ch x ch x samples]
% A - original signals                             [samples x ch]: complex 

A_orth = orthogonalize_signals(A, WindowLength, FidelityOperator);


d1 = size(A_orth, 1); % ch 
d2 = size(A_orth, 2); % ch
%d3 = size(A, 3); % length

N_Windows = fix(size(A, 1)/WindowLength);

A(WindowLength*N_Windows+1:end, :) = []; % cut extra samples

A = abs(A); % envelope from complex bandpass filtered signals

% Compute oCC for amp.frequency(k) 
oCC = zeros(d1, d2);

for ChA = 1:d1

    X = A(:, ChA); % orig signal

    for ChB = 1:d2

        Y = squeeze(A_orth(ChA, ChB, :)); % orth signals

        if isempty(FidelityOperator) || FidelityOperator(ChA, ChB)
            oCC(ChA, ChB) = corr(X, Y);    
        end

    end % ChB

end % ChA
   
end

function A_orth = orthogonalize_signals(A, WindowLength, FidelityOperator)

  % A - parcel/EEG signals [samples x ch] in complex form
  % A_orth - envelopes
  
    N_Windows = fix(size(A, 1)/WindowLength);
  
    A(WindowLength*N_Windows+1:end, :) = []; % cut extra samples
    
      L = size(A, 1);
    Nch = size(A, 2); 
      
    A_orth = zeros(Nch, Nch, L); % init 3D array 
   
    for ChA = 1:Nch % dim 1
     
         X = A(:, ChA);     % [samples x 1]
                
         for ChB = 1:Nch % dim 2
             
             Y = A(:, ChB); % [samples x 1]
             
             if isempty(FidelityOperator) || FidelityOperator(ChA, ChB)
            
              % Reshape arrays: 1 Column = 1 window  
                X_blocks = reshape(X, WindowLength, N_Windows);
                Y_blocks = reshape(Y, WindowLength, N_Windows);
                
              % Orthogonalize Y re X >> Y_orth put to A_orth(ChA, ChB, :)
              % Method based on Brookes et al. (2012) paper   
                B = sum((real(X_blocks)./repmat(sum(real(X_blocks).^2, 1), WindowLength, 1) ).*real(Y_blocks), 1);
       
                Y_orth_blocks = Y_blocks - repmat(B, WindowLength, 1) .* X_blocks;
       
                A_orth(ChA, ChB, :) = reshape(Y_orth_blocks, L, 1);
                 
             end
             
         end % ChB loop
         
    end % ChA loop
    
    A_orth = abs(A_orth); % envelopes from analytic signal
      
end

