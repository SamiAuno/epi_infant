%% Re-refrence to transverse bipolar montage
% Returns a structure eeg
% eeg.srate == sampling rate
% eeg.signal == bipolar signal 
% eeg.times == time array in seconds

function [eeg,spindels] = rerefBIP(sFile,Epochs_table)

% Get data matrix
[sMatrix, ~] = in_bst(sFile);
DataMat = in_bst_data( sFile,'F');

% Get channel names
ChannelFile = bst_get('ChannelFileForStudy', sFile);
ChannelMat = in_bst_channel(ChannelFile);

% Channel names
CH = {ChannelMat.Channel.Name}';

% Make anonymous function that returns the index of given channeö name from
% the CH cell array.
c = @(CH_name) find(strcmp(CH,CH_name));


% Make transverse bipolar montage.

eeg.signal(1,:) = sMatrix.F( c('F8'),:)  -  sMatrix.F( c('F6'),:);
eeg.signal(2,:) = sMatrix.F( c('F6'),:)  -  sMatrix.F( c('F4'),:);
eeg.signal(3,:) = sMatrix.F( c('F4'),:)  -  sMatrix.F( c('F2'),:);
eeg.signal(4,:) = sMatrix.F( c('F2'),:)  -  sMatrix.F( c('Fz'),:);
eeg.signal(5,:) = sMatrix.F( c('Fz'),:)  -  sMatrix.F( c('F1'),:);
eeg.signal(6,:) = sMatrix.F( c('F1'),:)  -  sMatrix.F( c('F3'),:);
eeg.signal(7,:) = sMatrix.F( c('F3'),:)  -  sMatrix.F( c('F5'),:);
eeg.signal(8,:) = sMatrix.F( c('F5'),:)  -  sMatrix.F( c('F7'),:);

eeg.signal(9,:)  = sMatrix.F( c('FT8'),:)  -  sMatrix.F( c('FC6'),:);
eeg.signal(10,:) = sMatrix.F( c('FC6'),:)  -  sMatrix.F( c('FC4'),:);
eeg.signal(11,:) = sMatrix.F( c('FC4'),:)  -  sMatrix.F( c('FC2'),:);
eeg.signal(12,:) = sMatrix.F( c('FC2'),:)  -  sMatrix.F( c('FCz'),:);
eeg.signal(13,:) = sMatrix.F( c('FCz'),:)  -  sMatrix.F( c('FC1'),:);
eeg.signal(14,:) = sMatrix.F( c('FC1'),:)  -  sMatrix.F( c('FC3'),:);
eeg.signal(15,:) = sMatrix.F( c('FC3'),:)  -  sMatrix.F( c('FC5'),:);
eeg.signal(16,:) = sMatrix.F( c('FC5'),:)  -  sMatrix.F( c('FT7'),:);

eeg.signal(17,:) = sMatrix.F( c('T8'),:)  -  sMatrix.F( c('C6'),:);
eeg.signal(18,:) = sMatrix.F( c('C6'),:)  -  sMatrix.F( c('C4'),:);
eeg.signal(19,:) = sMatrix.F( c('C4'),:)  -  sMatrix.F( c('C2'),:);
eeg.signal(20,:) = sMatrix.F( c('C2'),:)  -  sMatrix.F( c('Cz'),:);
eeg.signal(21,:) = sMatrix.F( c('Cz'),:)  -  sMatrix.F( c('C1'),:);
eeg.signal(22,:) = sMatrix.F( c('C1'),:)  -  sMatrix.F( c('C3'),:);
eeg.signal(23,:) = sMatrix.F( c('C3'),:)  -  sMatrix.F( c('C5'),:);
eeg.signal(24,:) = sMatrix.F( c('C5'),:)  -  sMatrix.F( c('T7'),:);

eeg.signal(25,:) = sMatrix.F( c('TP8'),:)  -  sMatrix.F( c('CP6'),:);
eeg.signal(26,:) = sMatrix.F( c('CP6'),:)  -  sMatrix.F( c('CP4'),:);
eeg.signal(27,:) = sMatrix.F( c('CP4'),:)  -  sMatrix.F( c('CP2'),:);
eeg.signal(28,:) = sMatrix.F( c('CP2'),:)  -  sMatrix.F( c('CP1'),:);
eeg.signal(29,:) = sMatrix.F( c('CP1'),:)  -  sMatrix.F( c('CP3'),:);
eeg.signal(30,:) = sMatrix.F( c('CP3'),:)  -  sMatrix.F( c('CP5'),:);
eeg.signal(31,:) = sMatrix.F( c('CP5'),:)  -  sMatrix.F( c('TP7'),:);

eeg.signal(32,:) = sMatrix.F( c('P8'),:)  -  sMatrix.F( c('P6'),:);
eeg.signal(33,:) = sMatrix.F( c('P6'),:)  -  sMatrix.F( c('P4'),:);
eeg.signal(34,:) = sMatrix.F( c('P4'),:)  -  sMatrix.F( c('P2'),:);
eeg.signal(35,:) = sMatrix.F( c('P2'),:)  -  sMatrix.F( c('Pz'),:);
eeg.signal(36,:) = sMatrix.F( c('Pz'),:)  -  sMatrix.F( c('P1'),:);
eeg.signal(37,:) = sMatrix.F( c('P1'),:)  -  sMatrix.F( c('P3'),:);
eeg.signal(38,:) = sMatrix.F( c('P3'),:)  -  sMatrix.F( c('P5'),:);
eeg.signal(39,:) = sMatrix.F( c('P5'),:)  -  sMatrix.F( c('P7'),:);

eeg.signal(40,:) = sMatrix.F( c('PO6'),:)  -  sMatrix.F( c('PO4'),:);
eeg.signal(41,:) = sMatrix.F( c('PO4'),:)  -  sMatrix.F( c('POz'),:);
eeg.signal(42,:) = sMatrix.F( c('POz'),:)  -  sMatrix.F( c('PO3'),:);
eeg.signal(43,:) = sMatrix.F( c('PO3'),:)  -  sMatrix.F( c('PO5'),:);

eeg.signal(44,:) = sMatrix.F( c('PO8'),:)  -  sMatrix.F( c('O2'),:);
eeg.signal(45,:) = sMatrix.F( c('O2'),:)  -  sMatrix.F( c('Oz'),:);
eeg.signal(46,:) = sMatrix.F( c('Oz'),:)  -  sMatrix.F( c('O1'),:);
eeg.signal(47,:) = sMatrix.F( c('O1'),:)  -  sMatrix.F( c('PO7'),:);

% Naming
eeg.chName{1,1} = 'F8-F6';
eeg.chName{2,1} = 'F6-F4';
eeg.chName{3,1} = 'F4-F2';
eeg.chName{4,1} = 'F2-Fz';
eeg.chName{5,1} = 'Fz-F1';
eeg.chName{6,1} = 'F1-F3';
eeg.chName{7,1} = 'F3-F5';
eeg.chName{8,1} = 'F5-F7';

eeg.chName{9,1}  = 'FT8-FC6';
eeg.chName{10,1} = 'FC6-FC4';
eeg.chName{11,1} = 'FC4-FC2';
eeg.chName{12,1} = 'FC2-FCz';
eeg.chName{13,1} = 'FCz-FC1';
eeg.chName{14,1} = 'FC1-FC3';
eeg.chName{15,1} = 'FC3-FC5';
eeg.chName{16,1} = 'FC5-FT7';

eeg.chName{17,1} = 'T8-C6';
eeg.chName{18,1} = 'C6-C4';
eeg.chName{19,1} = 'C4-C2';
eeg.chName{20,1} = 'C2-Cz';
eeg.chName{21,1} = 'Cz-C1';
eeg.chName{22,1} = 'C1-C3';
eeg.chName{23,1} = 'C3-C5';
eeg.chName{24,1} = 'C5-T7';

eeg.chName{25,1} = 'TP8-CP6';
eeg.chName{26,1} = 'CP6-CP4';
eeg.chName{27,1} = 'CP4-CP2';
eeg.chName{28,1} = 'CP2-CP1';
eeg.chName{29,1} = 'CP1-CP3';
eeg.chName{30,1} = 'CP3-CP5';
eeg.chName{31,1} = 'CP5-TP7';

eeg.chName{32,1} = 'P8-P6';
eeg.chName{33,1} = 'P6-P4';
eeg.chName{34,1} = 'P4-P2';
eeg.chName{35,1} = 'P2-Pz';
eeg.chName{36,1} = 'Pz-P1';
eeg.chName{37,1} = 'P1-P3';
eeg.chName{38,1} = 'P3-P5';
eeg.chName{39,1} = 'P5-P7';

eeg.chName{40,1} = 'PO6-PO4';
eeg.chName{41,1} = 'PO4-POz';
eeg.chName{42,1} = 'POz-PO3';
eeg.chName{43,1} = 'PO3-PO5';

eeg.chName{44,1} = 'PO8-O2';
eeg.chName{45,1} = 'O2-Oz';
eeg.chName{46,1} = 'Oz-O1';
eeg.chName{47,1} = 'O1-PO7';

% Calculate sampling rate.
eeg.srate = 1/(sMatrix.Time(2)-sMatrix.Time(1));

% Remove BAD segments from the data.
events = DataMat.F.events;
badEventLabels = find(contains({events.label},'BAD') | contains({events.label},'transient_bandpass'));
if badEventLabels
    
    bad_mask = true(1,size(sMatrix.Time,2));
    
    for i = badEventLabels
        event_times = events(i).times;
        event_times = round((event_times - sMatrix.Time(1)) * eeg.srate + 1);      % Change times to samples/indicis
        
        % loop through each
        for s = 1:size(event_times,2)
            startEvent = event_times(1,s);
            endEvent = event_times(2,s);
            
            if endEvent > length(bad_mask)
                bad_mask(1,startEvent:end) = false;
            else
                bad_mask(1,startEvent:endEvent) = false;
            end
            
        end
        
    end
    
    eeg.signal = eeg.signal(:,bad_mask);
    eeg.times = sMatrix.Time(1,bad_mask);
    
else
    eeg.times = sMatrix.Time;
end


%% Get Epochs
% Get the patient name
patient_name = getPatientName(sFile);

if ~isempty(patient_name)
    % Check that the patient has a corresponding epoch
    if any(strcmp(Epochs_table.PatientName,patient_name))
        
        OPTIONS.epochLength = [];              % Set to empty to get the whole epoch
        OPTIONS.epochSegmentLength = [];
        OPTIONS.spindles = [];                  % This is used later, not relevant now.
        
        [epochs_index_mask,~,~] = getEpochIndecis(Epochs_table,patient_name, eeg.times,eeg.srate,OPTIONS);
        
        
        % Collapse the epochs_index_mask to just one mask
        % If there were N number of sleep stages (epochs), this collapses
        % them all into the same mask.
        epochs_index_mask = logical(sum(epochs_index_mask,1));
        
        eeg.epochIndecis = find(epochs_index_mask);
        eeg.signal = eeg.signal(:,epochs_index_mask);
        eeg.times = eeg.times(1,epochs_index_mask);
    else
        warning(['Patient ',patient_name,' does not have sleep epochs.']);
    end
end



%% Add spindle events to the eeg structure
spindleEventLabel = find(contains({events.label},'spindle') | contains({events.label},'spindles') | contains({events.label},'spindel') | contains({events.label},'spindels'));
if spindleEventLabel
    eeg.spindleEvents = events(spindleEventLabel).times;
    if badEventLabels
        % Loop through each spindle event
        for spi = 1:size(eeg.spindleEvents,2)
            eeg.spindleEvents(1,spi) = find(eeg.spindleEvents(1,spi) <= eeg.times,1,'first');
            eeg.spindleEvents(2,spi) = find(eeg.spindleEvents(2,spi) <= eeg.times,1,'first');
        end
    else
        eeg.spindleEvents = round((eeg.spindleEvents - eeg.times(1)) * eeg.srate + 1);      % Change times to samples/indicis
    end
else
    eeg.spindleEvents = [];
end














end