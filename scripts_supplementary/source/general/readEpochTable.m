
function Epochs_table = readEpochTable(csvName)
% currentScriptPath = matlab.desktop.editor.getActiveFilename;
currentScriptPath = mfilename('fullpath');
[currentFolderPath,~,~] = fileparts(currentScriptPath);

Epochs_table = readtable(fullfile(currentFolderPath,csvName));
end