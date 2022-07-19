function Cx_atlas_split(MyAtlas,scalp,flat_cx)



% Load data
 %load('C:\Users\Anton\Documents\MATLAB\BabyConn Toolbox V3 src\src\MyAtlas.mat');
%   load('C:\Users\Anton\Documents\MATLAB\BabyConn Toolbox V4 src\src\MyAtlas_ROI.mat');
%   
%   load('C:\Users\Anton\Documents\MATLAB\BabyConn Toolbox V3 src\src\flat_cx.mat');
%   
%   load('C:\Users\Anton\Documents\MATLAB\BabyConn Toolbox V3 src\src\scalp.mat');
  
  
  
% ======================================================================= %  
  %Hemi = {'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'L'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'; 'R'};
  
  %Areas = {'T'; 'C'; 'T'; 'F'; 'O'; 'T'; 'C'; 'O'; 'C'; 'C'; 'X'; 'F'; 'F'; 'O'; 'F'; 'T'; 'F'; 'F'; 'O'; 'C'; 'X'; 'X'; 'F'; 'T'; 'O'; 'O'; 'F'; 'T'; 'C'; 'C'; 'T'; 'F'; 'T'; 'C'; 'T'; 'F'; 'O'; 'T'; 'C'; 'O'; 'C'; 'C'; 'X'; 'F'; 'F'; 'O'; 'F'; 'T'; 'F'; 'F'; 'O'; 'C'; 'X'; 'X'; 'F'; 'T'; 'O'; 'O'; 'F'; 'T'; 'C'; 'C'; 'T'; 'F'};

  %L_F2 = find(cell2mat(Areas) == 'F')'; 
  
  Area_F = find(cell2mat(MyAtlas.Areas) == 'F')'; % frontal
  Area_C = find(cell2mat(MyAtlas.Areas) == 'C')'; % central
  Area_T = find(cell2mat(MyAtlas.Areas) == 'T')'; % temporal
  Area_O = find(cell2mat(MyAtlas.Areas) == 'O')'; % ocipital
  
  Area_X = find(cell2mat(MyAtlas.Areas) == 'X')'; % bad
  
  figure;
  hold on
  
  set(gcf, 'Color', 'w');
  
  
  
  patch('Vertices', scalp.Vertices, 'Faces', scalp.Faces, 'FaceColor', [0.96 0.92 0.92], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
  
  
  
 % Frontal 
   for k = 1:length(Area_F)
       scatter3(flat_cx.Vertices(MyAtlas.Parcels{Area_F(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{Area_F(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{Area_F(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
   end

 % Central 
   for k = 1:length(Area_C)
       scatter3(flat_cx.Vertices(MyAtlas.Parcels{Area_C(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{Area_C(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{Area_C(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
   end  
   
 % Temporal 
   for k = 1:length(Area_T)
       scatter3(flat_cx.Vertices(MyAtlas.Parcels{Area_T(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{Area_T(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{Area_T(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
   end
 
 % Ocipital 
   for k = 1:length(Area_O)
       scatter3(flat_cx.Vertices(MyAtlas.Parcels{Area_O(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{Area_O(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{Area_O(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y');
   end 
   
 % Bad 
   for k = 1:length(Area_X)
       scatter3(flat_cx.Vertices(MyAtlas.Parcels{Area_X(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{Area_X(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{Area_X(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
   end
  
  
  
  
  
  axis off
  
  
   
    
%      
%    
%   
%   
%   Np = size(MyAtlas.Centroids, 1);
%   
%  % Frontal 
%    L_F = [4 12 13 15 17 18 23 27 32];
%    
%    L_F = [L_F L_F + 32];
%    
%  % Central 
%    L_C = [2 7 9 10 20 29 30];
%    
%    L_C = [L_C L_C + 32];
%    
%  % Temporal 
%    L_T = [1 3 6 16 24 28 31];
%    
%    L_T = [L_T L_T + 32];
%    
%  % Ocipital 
%    L_O = [5 8 14 19 25 26];
%    
%    L_O = [L_O L_O + 32];
%    
%  % Not included 
%    L_X = [11 21 22];
%    
%    L_X = [L_X L_X + 32];
%    
%    
%    %for k = 1:length(L_F)
%    %    scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'r');
%    %end
%   
%   figure;
%   hold on
%   
%   set(gcf, 'Color', 'w');
%   
%   
%   
%   patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 1.0, 'EdgeColor', 'none');
%   
%   
%   patch('Vertices', scalp.Vertices, 'Faces', scalp.Faces, 'FaceColor', [0.96 0.92 0.92], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
%   
%   
%   
%   scatter3(1.1*MyAtlas.Centroids(1:32, 1), 1.1*MyAtlas.Centroids(1:32, 2), 1.1*MyAtlas.Centroids(1:32, 3), 50, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g'); 
%   
%   
%   for j = 1:Np/2
%       text(1.15*MyAtlas.Centroids(j, 1), 1.15*MyAtlas.Centroids(j, 2), 1.15*MyAtlas.Centroids(j, 3), num2str(j), 'FontSize', 8, 'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
%   end 
%   
%  % L_Frontal 
%    for k = 1:length(L_F)
%        scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_F(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
%    end
% 
%  % L_Central 
%    for k = 1:length(L_C)
%        scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_C(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_C(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_C(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%    end  
%    
%  % L_Temporal 
%    for k = 1:length(L_T)
%        scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_T(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_T(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_T(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g');
%    end
%  
%  % L_Ocipital 
%    for k = 1:length(L_O)
%        scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_O(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_O(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_O(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'y');
%    end 
%    
%   % L_Bad 
%    for k = 1:length(L_X)
%        scatter3(flat_cx.Vertices(MyAtlas.Parcels{L_X(k), 1}, 1)  , flat_cx.Vertices(MyAtlas.Parcels{L_X(k), 1}, 2), flat_cx.Vertices(MyAtlas.Parcels{L_X(k), 1}, 3), 30, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
%    end  
%    
% 
%   
   
end 



 








