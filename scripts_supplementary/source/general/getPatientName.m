%% Get the patient name
% Use regular expression to find the subject/patient name
% This assumes that the source files are in the brainstorm folder
% structure and that the path is full path!
% The subject name is the string between two fileseparators ('/' in
% linux, '\' in windows) and after 'data'.

function patient_name = getPatientName(paths_sFiles)

% In Windows, the file separator is '\' which happens to be the same symbol
% used in regexp commands... Thus, if pc, change to linux convetion
if ispc
    paths_sFiles = strrep(paths_sFiles,'\','/');
end

expression = ['data','/','(?<name>\w*)','/'];
str = paths_sFiles;
string_name = regexp(str,expression,'names');
if isempty(string_name)
    disp('Patient name could not be parsed from the file path!')
    patient_name = [];
else
    patient_name = string_name.name;
end
    
end