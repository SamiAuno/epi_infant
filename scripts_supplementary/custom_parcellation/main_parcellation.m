%% main file to make custom parcellation for EPI Infant projects

%% Paths and options
% Path to smoothened and downsampled brainstorm cortex file
path2smoothCortex = '/home/sami/brainstorm_db/CNT_testaus/anat/@default_subject/tess_cortex_pial_10000V_smooth.mat';

path2scalp = '/home/sami/brainstorm_db/CNT_testaus/anat/@default_subject/tess_head_bem_1922V.mat';

scalp = load(path2scalp);
smoothCortex = load(path2smoothCortex);

% Number of parcels per hemisphere
N_hemiparcels = 33;

% Atlas name
nameNewAtlas = 'EPI-Infant';


%% Make parcellation
Parcellation = Cx_atlas(path2smoothCortex, N_hemiparcels);


%% Hemisphere and other info
% Run this part multiple times, everytime checking the parcellation
% indecis! Do it until you are happy.

Hemisphere = cell(2*N_hemiparcels,1);

Hemisphere(1:N_hemiparcels) = {'L'};
Hemisphere(N_hemiparcels+1:end) = {'R'};
Areas = cell(2*N_hemiparcels,1);

% These indecis are for Parcellation_20191101
% frontalInd = [1,2,4,9,10,12,13,15,18,19,20,24,26]; % these indeces are only for one hemi
% centralInd = [5,7,17,21,29,31];
% temporalInd = [8,25,28];
% occipitalInd = [3,11,14,16,22,30,32];
% removeInd = [6,23,27];

% These indecis are for Parcellation_20191125
frontalInd = [2,16,18,21,22,27,31,33]; % these indeces are only for one hemi
centralInd = [4,11,12,13,20,25,28,30,32];
temporalInd = [6,10,14,19,29];
occipitalInd = [1,5,8,15,17,23,24];
removeInd = [3,7,9,26];


Areas([frontalInd,frontalInd+N_hemiparcels],1) = {'F'};
Areas([centralInd,centralInd+N_hemiparcels],1) = {'C'};
Areas([temporalInd,temporalInd+N_hemiparcels],1) = {'T'};
Areas([occipitalInd,occipitalInd+N_hemiparcels],1) = {'O'};
Areas([removeInd,removeInd+N_hemiparcels],1) = {'X'};

Parcellation.Hemisphere = Hemisphere;
Parcellation.Areas = Areas;

% Visualise
Cx_atlas_split(Parcellation,scalp,smoothCortex);



%% Make new atlas

% Make new color palet
frontalCol = [randi([80,100],[length(frontalInd),1]),randi([0,20],[length(frontalInd),1]),randi([0,20],[length(frontalInd),1])]/100;
centralCol = [randi([0,20],[length(centralInd),1]),randi([0,20],[length(centralInd),1]),randi([80,100],[length(centralInd),1])]/100;
temporalCol = [randi([0,20],[length(temporalInd),1]),randi([80,100],[length(temporalInd),1]),randi([0,20],[length(temporalInd),1])]/100;
occipitalCol = [randi([80,100],[length(occipitalInd),1]),randi([80,100],[length(occipitalInd),1]),randi([0,20],[length(occipitalInd),1])]/100;
removedCol = [0,0,0];

% Parcellation.Scouts([frontalInd,frontalInd+N_hemiparcels],1).


newAtlas.Name = nameNewAtlas;
newAtlas.Scouts = Parcellation.Scouts;

% Make color
for co = 1:length(frontalInd)
    newAtlas.Scouts(frontalInd(co)).Color = frontalCol(co,:);
    newAtlas.Scouts(frontalInd(co)+N_hemiparcels).Color = frontalCol(co,:);    
end
for co = 1:length(centralInd)
    newAtlas.Scouts(centralInd(co)).Color = centralCol(co,:);
    newAtlas.Scouts(centralInd(co)+N_hemiparcels).Color = centralCol(co,:);    
end
for co = 1:length(temporalInd)
    newAtlas.Scouts(temporalInd(co)).Color = temporalCol(co,:);
    newAtlas.Scouts(temporalInd(co)+N_hemiparcels).Color = temporalCol(co,:);    
end
for co = 1:length(occipitalInd)
    newAtlas.Scouts(occipitalInd(co)).Color = occipitalCol(co,:);
    newAtlas.Scouts(occipitalInd(co)+N_hemiparcels).Color = occipitalCol(co,:);    
end
for co = 1:length(removeInd)
    newAtlas.Scouts(removeInd(co)).Color = removedCol;
    newAtlas.Scouts(removeInd(co)+N_hemiparcels).Color = removedCol;    
end

%% Make Label, Function, Region and Handles
for i = 1:length(newAtlas.Scouts)
    
    if i <= N_hemiparcels
        hemis = 'L';
    else
        hemis = 'R';
    end
    
    
    % Make label
    if find(ismember([frontalInd,frontalInd+N_hemiparcels],i))
        
        newAtlas.Scouts(i).Label = [sprintf( '%03d', i ),'_',hemis,'_frontal',];
        
    elseif find(ismember([centralInd,centralInd+N_hemiparcels],i))
        
        newAtlas.Scouts(i).Label = [sprintf( '%03d', i ),'_',hemis,'_central',];
        
    elseif find(ismember([temporalInd,temporalInd+N_hemiparcels],i))
            
        newAtlas.Scouts(i).Label = [sprintf( '%03d', i ),'_',hemis,'_temporal',];
        
    elseif find(ismember([occipitalInd,occipitalInd+N_hemiparcels],i))
                
        newAtlas.Scouts(i).Label = [sprintf( '%03d', i ),'_',hemis,'_occipital',];
        
    elseif find(ismember([removeInd,removeInd+N_hemiparcels],i))
        
        newAtlas.Scouts(i).Label = '000';
        
    end
    
    newAtlas.Scouts(i).Function = 'Mean';
    newAtlas.Scouts(i).Region = [hemis,'U'];
    newAtlas.Scouts(i).Handles = [];
end

% Total number of vertices
newAtlas.TessNbVertices = length([newAtlas.Scouts.Vertices]);


% Remove scouts that are marked to be removed.
newAtlas.Scouts([removeInd,removeInd+N_hemiparcels]) = [];

%% Save as Brainstorm compatible scouts

save([newAtlas.Name,'_',num2str(2*N_hemiparcels),'_scout.mat'],'-struct','newAtlas');






