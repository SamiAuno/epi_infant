%% Process_spindles
% process_spindles(sFiles)
% This script utilises Vivianas script to first find spindels in each
% channel, then it groups spindels by some rules, and saves them to a file.

% sFiles are Brainstorm filepaths to preprocessed data files

% The data needs to have DC removed, Broad-band band-pass
% filtered and downsampled before this process




function process_spindles(sFiles)

if ~iscell(sFiles)
    sFiles = {sFiles};
end

N_sFiles = length(sFiles);

% Process each file separately
for ns = 1:N_sFiles
    csvName = 'Analysis_Epochs.csv';
    Epochs_table = readEpochTable(csvName);
    
    patient_name = getPatientName(sFiles{ns});
    
    if ~any(strcmp(Epochs_table.PatientName,patient_name))
        disp(['Patient ',patient_name,' does not have sleep epochs. Skip patient.']);
        continue
    end
    
    % Re-reference to transverse bipolar montage
    eeg = rerefBIP(sFiles{ns},Epochs_table);

    % Extract spindels
    [eeg, OUTPUT] = extract_spindles(eeg);
    
    % Group and threshold spindles into groups
    OUTPUT = removeLowPowerSpindles(OUTPUT,eeg);
    
    OUTPUT.removedLength = sum([OUTPUT.groups.length_seconds])/60;
    OUTPUT.remaininLength = (length(eeg.epochIndecis)/250)/60 - OUTPUT.removedLength;
    
    spindles.(patient_name) = OUTPUT;
    
    
    
end

fileSaveName = 'Spindles.mat';

save(fileSaveName,'spindles');

end











