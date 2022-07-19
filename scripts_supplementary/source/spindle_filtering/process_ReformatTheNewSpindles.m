
function spindles = process_ReformatTheNewSpindles()

% Get the path to the raw files
sFiles = getSFiles();

N_sFiles = length(sFiles);

csvName = 'Analysis_Epochs.csv';
Epochs_table = readEpochTable(csvName);

global path_headModel;
global path_inverseModel;
[path_headModel, path_inverseModel] = getPathsToHeadAndInverseModels();

% spindleRunFolder = '/media/sami/ADATA SD700/EPI_INFNT_BS/parcel_spindles/Spindle_run_210305_1814';
spindleRunFolder = 'D:\EPI_INFNT_BS\parcel_spindles\Spindle_run_210305_1814';

spindleFiles = dir([spindleRunFolder,filesep,'*.mat']);

% Process each file separately
for ns = 1:N_sFiles
    patient_name = getPatientNameFromPath(sFiles{ns});
    
    if ~any(strcmp(Epochs_table.PatientName,patient_name))
        disp(['Patient ',patient_name,' does not have sleep epochs. Skip patient.']);
        continue
    end
    
    % load the pre-calculated spindles
    spindleFileIndex = find(contains({spindleFiles.name}',patient_name),1);
    if isempty(spindleFileIndex)
        disp(['Patient ',patient_name,' does not have pre-calculated spindles. Skip patient.']);
        continue
    end
    load(fullfile(spindleRunFolder,spindleFiles(spindleFileIndex).name),'spindle_det');
    % Fetch the eeg
%     eeg = rerefBIP(sFiles{ns},Epochs_table);
    
    % Get the parcel level signal
    [eeg.signal,eeg.times] = getParcelTimeSeries(sFiles{ns});
    eeg.spindles = eeg.signal;
    eeg.srate = 1/(eeg.times(2) - eeg.times(1));
    eeg.epochIndecis = getEpochIndecisFromTable(Epochs_table,patient_name, eeg);
    % Extract spindels
%     [eeg, OUTPUT] = extract_spindles(eeg);
    % Rearrange the spindles to this format
    OUTPUT = rearrangeSpindles(spindle_det);
    
    % Group and threshold spindles into groups
    OUTPUT = removeLowPowerSpindles(OUTPUT,eeg);
    
    OUTPUT.removedLength = sum([OUTPUT.groups.length_seconds])/60;
    OUTPUT.remaininLength = (length(eeg.epochIndecis)/250)/60 - OUTPUT.removedLength;
    
    spindles.(patient_name) = OUTPUT;
    
    
    
end

fileSaveName = 'EPI_Infant_Spindles_NewSet.mat';

save(fileSaveName,'spindles');

end

function epochIndecis = getEpochIndecisFromTable(Epochs_table,patient_name, eeg)
OPTIONS.epochLength = [];              % Set to empty to get the whole epoch
OPTIONS.epochSegmentLength = [];
OPTIONS.spindles = [];                  % This is used later, not relevant now.
[epochs_index_mask,~,~] = getEpochIndecis(Epochs_table,patient_name, eeg.times,eeg.srate,OPTIONS);

% Collapse the epochs_index_mask to just one mask
% If there were N number of sleep stages (epochs), this collapses
% them all into the same mask.
epochs_index_mask = logical(sum(epochs_index_mask,1));

epochIndecis = find(epochs_index_mask);
end


function OUTPUT = rearrangeSpindles(spindle_det)

maxNumberOfSpindles = max(cellfun(@length,({spindle_det.startSample})));
numberOfParcel = length(spindle_det);

OUTPUT.start = nan(numberOfParcel,maxNumberOfSpindles);
OUTPUT.end = nan(numberOfParcel,maxNumberOfSpindles);


for parc_i = 1:numberOfParcel
    nSpindlesInParcel = length(spindle_det(parc_i).startSample);
    if nSpindlesInParcel > 0
        OUTPUT.start(parc_i,1:nSpindlesInParcel) = spindle_det(parc_i).startSample;
        OUTPUT.end(parc_i,1:nSpindlesInParcel) = spindle_det(parc_i).endSample;
    end
end



end

function patient_name = getPatientNameFromPath(paths_sFiles)

% In Windows, the file separator is '\' which happens to be the same symbol
% % used in regexp commands... Thus, if pc, change to linux convetion
if ispc
    paths_sFiles = strrep(paths_sFiles,'\','/');
end

expression = ['(?<name>\w*)','/@.'];
str = paths_sFiles;
string_name = regexp(str,expression,'names');
if isempty(string_name)
    disp(paths_sFiles);
    disp('Patient name could not be parsed from the file path!')
    patient_name = [];
else
    patient_name = string_name.name;
end
    
end

function [parcelTimeSeries,times] = getParcelTimeSeries(sFiles)
global path_headModel;
global path_inverseModel;
% Read data
kernel = load(path_inverseModel);
channelIndecis = kernel.GoodChannel;
[sMatrix, ~] = in_bst(sFiles);
sMatrix.sFreq = length(sMatrix.Time)/(sMatrix.Time(end)-sMatrix.Time(1));   % Sampling frequency
sMatrix.F = sMatrix.F(channelIndecis,:);    % Take the appropriate channels

times = sMatrix.Time;

[filtered_data, ~, ~] = bst_bandpass_hfilter(sMatrix.F, sMatrix.sFreq, 10, 16, 0, true, [], 0, []);
% Get the parcel level signal
parcelTimeSeries = process_collapse2parcels(sFiles,filtered_data, path_headModel,path_inverseModel,0);
end


function paths_sFiles = getSFiles()

paths_sFiles = {...
     'EPI_01_1v/@rawEPI_1_2017-06-29_12-29-07_band_06_bl_resample_interpbad/data_0raw_EPI_1_2017-06-29_12-29-07_band_06_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_02_1v/@rawEPI_2_Egal__Ayanle_2017-06-08_11-36-20_band_07_bl_resample_interpbad/data_0raw_EPI_2_Egal__Ayanle_2017-06-08_11-36-20_band_07_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_03_1v/@rawEPI_3_Hawraa_Madalawi_2017-08-24_12-20-49_band_05_bl_resample_interpbad/data_0raw_EPI_3_Hawraa_Madalawi_2017-08-24_12-20-49_band_05_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_04_1v/@rawEPI_4_Huhta_Dami_2017-11-02_12-48-37_band_04_bl_resample_interpbad/data_0raw_EPI_4_Huhta_Dami_2017-11-02_12-48-37_band_04_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_05_1v/@rawEPI_5_Puchkina__Evgenia_2017-09-14_12-17-58_band_04_bl_resample_interpbad/data_0raw_EPI_5_Puchkina__Evgenia_2017-09-14_12-17-58_band_04_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_06_1v/@rawEPI_6_1v_EELI_2017-10-19_12-52-37_band_02_bl_resample_interpbad_02/data_0raw_EPI_6_1v_EELI_2017-10-19_12-52-37_band_02_bl_resample_interpbad_02_2021-02-18.mat';...
     'EPI_07_1v/@rawEPI_7_Vaino_Rintasaari_2017-08-17_12-51-01_band_04_bl_resample_interpbad/data_0raw_EPI_7_Vaino_Rintasaari_2017-08-17_12-51-01_band_04_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_08_1v/@rawEPI_8_ZAVORINA_ANASTASIA_2018-03-08_10-19-44_band_02_bl_resample_interpbad/data_0raw_EPI_8_ZAVORINA_ANASTASIA_2018-03-08_10-19-44_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_09_1v/@rawEPI_9_LUKAS_2018-05-17_09-46-07_band_03_bl_resample_interpbad/data_0raw_EPI_9_LUKAS_2018-05-17_09-46-07_band_03_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_10_1v/@rawEPI_10_Karhu__Nelli_2017-09-21_10-15-36_band_02_bl_resample_interpbad/data_0raw_EPI_10_Karhu__Nelli_2017-09-21_10-15-36_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_11_1v/@rawEPI_11_Lilja_Alexandersson_2017-12-07_12-38-01_band_02_bl_resample_interpbad/data_0raw_EPI_11_Lilja_Alexandersson_2017-12-07_12-38-01_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_12_1v/@rawEPI_12_Geissler__Fjella_2017-09-21_12-12-42_band_02_bl_resample_interpbad/data_0raw_EPI_12_Geissler__Fjella_2017-09-21_12-12-42_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_16_1v/@rawEPI_16_Rezkallah__Khaled_2018-03-15_11-28-19_band_02_bl_resample_interpbad/data_0raw_EPI_16_Rezkallah__Khaled_2018-03-15_11-28-19_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_17_1v/@rawEPI_17_ANNI_2018-06-07_14-55-27_band_02_bl_resample_interpbad/data_0raw_EPI_17_ANNI_2018-06-07_14-55-27_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_18_1v/@rawEPI_18_1V_RONJA_2018-08-16_14-34-43_band_02_bl_resample_interpbad/data_0raw_EPI_18_1V_RONJA_2018-08-16_14-34-43_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_19_1v/@rawEPI19_JASON_2018-06-07_12-11-04_band_02_bl_resample_interpbad/data_0raw_EPI19_JASON_2018-06-07_12-11-04_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_21_1v/@rawEPI21_JUHO_2018-06-13_15-56-33_band_02_bl_resample_interpbad/data_0raw_EPI21_JUHO_2018-06-13_15-56-33_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_22_1v/@rawepi22_justus_2018-10-11_12-43-55_band_02_bl_resample_interpbad/data_0raw_epi22_justus_2018-10-11_12-43-55_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_23_1v/@rawEPI23_4_PIHLA_2018-07-26_14-33-14_band_02_bl_resample_interpbad/data_0raw_EPI23_4_PIHLA_2018-07-26_14-33-14_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_25_1v/@rawEPI25_VILJO_2018-07-26_12-08-31_band_02_bl_resample_interpbad/data_0raw_EPI25_VILJO_2018-07-26_12-08-31_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_26_1v/@rawEPI26_ASSI_2018-06-14_10-11-36_band_02_bl_resample_interpbad/data_0raw_EPI26_ASSI_2018-06-14_10-11-36_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_28_1v/@rawisla_EPI28_1V_2018-10-23_16-13-05_band_02_bl_resample_interpbad/data_0raw_isla_EPI28_1V_2018-10-23_16-13-05_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_29_1v/@rawEPI_29_1y_band_02_bl_resample_interpbad/data_0raw_EPI_29_1y_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_30_1v/@rawEPI_30_IISA_2018-05-17_13-21-46_band_02_bl_resample_interpbad/data_0raw_EPI_30_IISA_2018-05-17_13-21-46_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_31_1v/@rawEPI31_4_STELLA_2019-02-21_11-13-11_band_02_bl_resample_interpbad/data_0raw_EPI31_4_STELLA_2019-02-21_11-13-11_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_32_1v/@rawEPI_32_1v_2018-12-20_12-42-22_band_02_bl_resample_interpbad_02/data_0raw_EPI_32_1v_2018-12-20_12-42-22_band_02_bl_resample_interpbad_02_2021-02-18.mat'; ...
     'EPI_33_1v/@rawEPI33_4_MASSINISSA_2019-03-07_11-52-05_band_02_bl_resample_interpbad/data_0raw_EPI33_4_MASSINISSA_2019-03-07_11-52-05_band_02_bl_resample_interpbad_2021-02-10.mat';...
     'EPI_34_1v/@rawEPI34_4_HERTTA_2019-02-28_10-52-09_band_bl_resample_03_interpbad_02/data_0raw_EPI34_4_HERTTA_2019-02-28_10-52-09_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_35_1v/@rawEPI35_3_KURTEN_LOVA_2018-12-13_12-29-02_band_bl_resample_05_interpbad_02/data_0raw_EPI35_3_KURTEN_LOVA_2018-12-13_12-29-02_band_bl_resample_05_interpbad_02_2021-02-10.mat';...
     'EPI_36_1v/@rawEPI36_5_stina_2019-08-13_16-36-05_band_bl_resample_03_interpbad_02/data_0raw_EPI36_5_stina_2019-08-13_16-36-05_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_37_1v/@rawEPI37_4_1V_TOMI_2019-01-31_12-43-19_band_bl_resample_03_interpbad_02/data_0raw_EPI37_4_1V_TOMI_2019-01-31_12-43-19_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_38_1v/@rawEPI38_MATEO_2019-01-10_12-16-23_band_bl_resample_03_interpbad_02/data_0raw_EPI38_MATEO_2019-01-10_12-16-23_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_40_1v/@rawEPI40_4_WILHELMIINA_2019-03-27_16-22-57_band_bl_resample_03_interpbad_02/data_0raw_EPI40_4_WILHELMIINA_2019-03-27_16-22-57_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_42_1v/@rawEPI42_ELIN_2018-10-18_12-22-45_band_bl_resample_05_interpbad_02/data_0raw_EPI42_ELIN_2018-10-18_12-22-45_band_bl_resample_05_interpbad_02_2021-02-10.mat';...
     'EPI_43_1v/@rawEPI_43_1v_uni_band_bl_resample_03_interpbad_02/data_0raw_EPI_43_1v_uni_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_44_1v/@rawEPI44_4_VAINO_2019-05-23_12-01-38_band_bl_resample_03_interpbad_02/data_0raw_EPI44_4_VAINO_2019-05-23_12-01-38_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_46_1v/@rawEPI46_1v_Segment_0_band_03_bl_resample_interpbad/data_0raw_EPI46_1v_Segment_0_band_03_bl_resample_interpbad_2021-02-18.mat';...
     'EPI_47_1v/@rawEPI47_VALTTERI_2019-03-21_09-07-21_band_bl_resample_03_interpbad_02/data_0raw_EPI47_VALTTERI_2019-03-21_09-07-21_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_50_1v/@rawEPIvauva50_Liljeqvist_Lukas_2019-08-20_16-39-19_band_bl_resample_03_interpbad_02/data_0raw_EPIvauva50_Liljeqvist_Lukas_2019-08-20_16-39-19_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_52_1v/@rawEPI52_URHO_2019-06-19_16-34-03_band_bl_resample_03_interpbad_02/data_0raw_EPI52_URHO_2019-06-19_16-34-03_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_54_1v/@rawEPI54_4_Seela_2019-09-05_10-10-20_band_bl_resample_03_interpbad_02/data_0raw_EPI54_4_Seela_2019-09-05_10-10-20_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_55_1v/@rawEPI_55_1y.edf_band_bl_resample_interpbad_02/data_0raw_EPI_55_1y.edf_band_bl_resample_interpbad_02_2021-02-10.mat';...
     'EPI_56_1v/@rawepi_56_3_Adrian_2019-08-08_13-38-51_band_bl_resample_03_interpbad_02/data_0raw_epi_56_3_Adrian_2019-08-08_13-38-51_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_57_1v/@rawepi57_3_samira_2019-08-01_12-16-34_band_bl_resample_03_interpbad_02/data_0raw_epi57_3_samira_2019-08-01_12-16-34_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_59_1v/@rawEPI_59_1v_Alex_2019-10-17_12-04-52_band_bl_resample_03_interpbad_02/data_0raw_EPI_59_1v_Alex_2019-10-17_12-04-52_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_60_1v/@rawEPI_60_1y_band_bl_resample_interpbad_02/data_0raw_EPI_60_1y_band_bl_resample_interpbad_02_2021-02-10.mat';...
     'EPI_61_1v/@rawEPI_61_1V_AIDAN_2019-10-10_12-47-09_band_bl_resample_03_interpbad_02/data_0raw_EPI_61_1V_AIDAN_2019-10-10_12-47-09_band_bl_resample_03_interpbad_02_2021-02-10.mat';...
     'EPI_62_1v/@rawEPI_62_1y.edf_band_bl_resample_interpbad_02/data_0raw_EPI_62_1y.edf_band_bl_resample_interpbad_02_2021-02-10.mat';...
     'EPI_63_1v/@rawEPI_63_1y.edf_band_bl_resample_interpbad_02/data_0raw_EPI_63_1y.edf_band_bl_resample_interpbad_02_2021-02-10.mat'};
 
 if ispc
     paths_sFiles = strrep(paths_sFiles,'/','\');
 end

end