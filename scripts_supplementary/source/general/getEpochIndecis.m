%% Get epochs for this patient

function [epochs_index_mask,epochs_timeStamp,Epoch_type_array] = getEpochIndecis(Epochs_table,patient_name, timeArray,sFreq,OPTIONS)

Epochs = table2cell(Epochs_table);
size_Epochs_table = size(Epochs_table);

patient_Epoch_indecis = find(contains(Epochs_table.PatientName,patient_name));

% Epoch type
Epoch_type_array = Epochs(patient_Epoch_indecis,2);

% List of unique states this patient has.
uniqStates = unique(Epoch_type_array);

% Number of epochs
N_epochs = length(patient_Epoch_indecis);       % Number of all epochs
N_UniqEpochs = length(uniqStates);              % Number of unique epochs
disp([patient_name,': ', num2str(N_UniqEpochs) ,' unique epochs found.']);

% Epoch sample mask. A vector of 0s and 1s. Used to mask out those
% samples that are not part of the epoch.
epochs_index_mask = zeros(N_UniqEpochs,length(timeArray));  

% Epoch durations array.
Epoch_time_array = reshape([Epochs{patient_Epoch_indecis,3:end}],N_epochs,size_Epochs_table(2)-2);

% If spindles calculated, remove them from epochs. Here, get the spindel
% mask
if ~isempty(OPTIONS.spindles)
    spindelMask = false(1,length(timeArray));
    if isfield(OPTIONS.spindles,patient_name)
        spindelIndecis = OPTIONS.spindles.(patient_name).spindelIndecis;
        spindelMask(spindelIndecis) = true;
        disp([patient_name,': Spindles...']);
    else
        disp([patient_name,': No spindles found.']);
    end
end


% For time stamping. Note, these are only the first segments of each
% epoch.
epochs_timeStamp = Epoch_time_array(:,1:2);

% Convert times to seconds
epochs_seconds              = seconds(Epoch_time_array);                 % Epoch start and end times in seconds
epochs_seconds_from_start   = epochs_seconds - timeArray(1);          % How many seconds from the start of the measurement

% Loop through unique epochs of this patient
for epoch_i = 1:N_UniqEpochs
%         N_epoch_segments = nnz(epochs_seconds_from_start(epoch_i,1:2:end));   % nnz = number of nonzero elements. Epoch doesn't need to be continuous, in which case it is in multiple shorter segments. This check how many segments there are in this epoch
    epochs_index_mask_temp = zeros(1,length(timeArray)); 

    % Loop through epochs of the same epoch segment and concanate them
    % together
    epoch_segment_indices = find(contains(Epoch_type_array,uniqStates{epoch_i}));
    N_epoch_segments = length(epoch_segment_indices);
    for epoch_segment_i = 1:N_epoch_segments

        segment_indecis = round(epochs_seconds_from_start(epoch_segment_indices(epoch_segment_i),:).* sFreq);

        epochs_index_mask_temp(1, segment_indecis(1):segment_indecis(2)) = 1; 

    end
    
    % If spindles already calculated and saved, remove from epochs_index_mask_temp
    if ~isempty(OPTIONS.spindles)
        epochs_index_mask_temp(1,spindelMask) = 0;
        disp([patient_name,': ... removed.']);
    end
    
    % This is now only the indecis of those samples that are wihtin the epoch.
    epochs_indecis_temp = find(epochs_index_mask_temp);
    
    % If the epochLength is not defined, return the whole epoch
    if isempty(OPTIONS.epochLength)
        epochs_index_mask(epoch_i, epochs_indecis_temp) = 1;
    else
    
        % Next we segment the concanated epoch into equidistant segments of length OPTIONS.epochSegmentLength
        totLengthEpoch = length(epochs_indecis_temp);

        % Take the central index of each equidistant segment
        segLength = floor((OPTIONS.epochSegmentLength * sFreq)/2)*2 + 1;    % The segment length as number of samples rounded to nearest odd integer.
        N_segments_in_epoch = OPTIONS.epochLength/OPTIONS.epochSegmentLength;
        centralIndexSegment = round(linspace(1 + (segLength - 1)/2, totLengthEpoch - (segLength - 1)/2, N_segments_in_epoch));

        % Loop through each segment of this epoch and fill it to the epochs_index_mask
        for segment_i = 1:N_segments_in_epoch
            ci = centralIndexSegment(segment_i);    % central index
            epochs_index_mask(epoch_i, epochs_indecis_temp(ci-(segLength - 1)/2 : ci+(segLength - 1)/2) ) = 1;
        end
    end
end

% Make the epoch index mask logical array so that it can be used in
% indexing.
epochs_index_mask = logical(epochs_index_mask);

% Fill this structure with info that can be used later.
patientEpochInformation.Epoch_type_array = Epoch_type_array;

end
    
