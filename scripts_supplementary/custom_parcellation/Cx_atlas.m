% Make custom parcellation atlas
% fullpath = path to brainstorm cortex file
% eg. '/home/sami/brainstorm_db/CNT_testaus/anat/@default_subject/tess_cortex_pial_10000V_smooth.mat'
% N = parcels per hemisphere

function Parcellation = Cx_atlas(path, N)


load(path);


Atlas = zeros(size(Vertices, 1), 1);



  figure;
  hold on
  set(gcf, 'Color', 'w');

% axis off

%patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
%patch('Vertices', conv_cx.Vertices, 'Faces', conv_cx.Faces, 'FaceColor', [1.0 0.0 0.0], 'FaceAlpha', 0.6, 'EdgeColor', 'none');

 

  scatter3(Vertices(:, 1), Vertices(:, 2), Vertices(:, 3), '.', 'k');

  plot3([-0.06 0.06], [0 0], [0 0], 'r');
  plot3([0 0], [-0.06 0.06], [0 0], 'r');
  plot3([0 0], [0 0], [0 0.06], 'r');

% Take only 1 hemisphere
  ind_Hemi_L = find(Vertices(:, 2) > 0);
  
  Hemi_L = Vertices(Vertices(:, 2) > 0, :); % y-axis
  
  
  ind_Hemi_R = find(Vertices(:, 2) < 0);
  
  Hemi_R = Vertices(Vertices(:, 2) < 0, :); % y-axis
  
  
  
  
  scatter3(Hemi_L(:, 1), Hemi_L(:, 2), Hemi_L(:, 3), 'o', 'r');
  
  hold off

 
  
   
  
  
% Cluster one Hemisphere  
  [ind_L, C_L, ~, D_L] = kmeans(Hemi_L, N, 'MaxIter', 5000); 
  
  Atlas(ind_Hemi_L) = ind_L;
  
  [~, ind_seed_L] = min(D_L);
  
  seed_L = zeros(N, 1);
  
  
  
  
% Find seeds for parcels (L)  
  for j = 1:N
      
      buf = zeros(size(ind_L, 1), 1);
      
      buf(ind_seed_L(j)) = 1;
      
      seed_buf = zeros(size(Vertices, 1), 1);
      
      seed_buf(ind_Hemi_L) = buf;
      
      seed_L(j) = find(seed_buf == 1);
      
  end
  
  
  figure;
  hold on
  set(gcf, 'Color', 'w');
  
%   xlim([-0.06 0.06]);
%   ylim([-0.06 0.06]);
%   zlim([0 0.08]);
  
  % scatter3(Hemi_L(:, 1), Hemi_L(:, 2), Hemi_L(:, 3), '.', 'k');
  % scatter3(C_L(:, 1), C_L(:, 2), C_L(:, 3), 'o', 'r');
  
  
  
  clr = rand(N, 3); % colors
    
  for k = 1:N
      
    % Parcel  
      scatter3(-Hemi_L(ind_L == k, 1),Hemi_L(ind_L == k, 3), -Hemi_L(ind_L == k, 2) , 'o', 'MarkerFaceColor', clr(k, :));
    % Centroid  
      scatter3(-C_L(k, 1), C_L(k, 3),-C_L(k, 2), 'o', 'MarkerFaceColor', clr(k, :));
      
  end
  
  hold off
  
  
 
  
  
  
  % Mirror Centroids to another Hemisphere
  
  C_R = C_L;
  C_R(:, 2) = C_R(:, 2) .* (-1);
  
  

  figure;
  hold on
  set(gcf, 'Color', 'w'); 
  
  patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
  
  scatter3(C_L(:, 1), C_L(:, 2), C_L(:, 3), 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
  
  scatter3(C_R(:, 1), C_R(:, 2), C_R(:, 3), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
 
  C = [C_L; C_R];
  
  for j = 1:size(C, 1)
      text(1.1*C(j, 1), 1.1*C(j, 2), 1.1*C(j, 3), num2str(j), 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
  end
  
  
  hold off
  
  
  
  
  
  
  
% Cluster another Hemisphere
  [ind_R, ~, ~, D_R] = kmeans(Hemi_R, N, 'MaxIter', 1, 'Start', C_R); 
  %ind_R = kmeans(Hemi_R, N, 'MaxIter', 1, 'Start', C_R);
  
    
  [~, ind_seed_R] = min(D_R);
  
  seed_R = zeros(N, 1);
  
  
  
  
% Find seeds for parcels (R)  
  for j = 1:N
      
      buf = zeros(size(ind_R, 1), 1);
      
      buf(ind_seed_R(j)) = 1;
      
      seed_buf = zeros(size(Vertices, 1), 1);
      
      seed_buf(ind_Hemi_R) = buf;
      
      seed_R(j) = find(seed_buf == 1);
      
  end
  
  

% Seeds together
% ======================================================================= %
  seeds = [seed_L; seed_R];
% ======================================================================= %  
  
   
  figure;
  hold on
  set(gcf, 'Color', 'w');
  
  xlim([-0.08 0.08]);
  ylim([-0.08 0.08]);
  zlim([0 0.12]);
  
  patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.4, 'EdgeColor', 'none');
  
  clr = rand(N, 3); % colors
    
  for k = 1:N
      
    % Parcel  
      scatter3(Hemi_L(ind_L == k, 1), Hemi_L(ind_L == k, 2), Hemi_L(ind_L == k, 3), 100, 'o', 'MarkerFaceColor', clr(k, :));
      scatter3(Hemi_R(ind_R == k, 1), Hemi_R(ind_R == k, 2), Hemi_R(ind_R == k, 3), 100, 'o', 'MarkerFaceColor', clr(k, :));
      
    % Centroid  
      %%%scatter3(C_L(k, 1), C_L(k, 2), C_L(k, 3), 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none');
      %%%scatter3(C_R(k, 1), C_R(k, 2), C_R(k, 3), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
      
  end
 
  axis off
  hold off
  
  
   
  
  Atlas(ind_Hemi_R) = ind_R + N;
  

% % Pack parcels  
  parcels{2*N, 1} = [];
% % ======================================================================= %  
  for k = 1:2*N
      parcels{k, 1} = find(Atlas == k)';
  end
% % ======================================================================= %
% 
% 


  for k = 1:2*N
      
      Scouts(k).Vertices = parcels{k, 1};
      Scouts(k).Seed = seeds(k, 1);
      
  end
  
Parcellation.Centroids = C;
Parcellation.Parcels = parcels;  
Parcellation.Scouts = Scouts;

%     
%  
% 
% % Feed to Brainstorm
% 
%   load('C:\Users\Anton\Documents\MATLAB\Brain_model_test_CSD\scout_scouts_Atlas_rand_symm_64_8000.mat');
%   
%   for k = 1:2*N
%       
%       Scouts(k).Vertices = parcels{k, 1};
%       Scouts(k).Seed = seeds(k, 1);
%       
%   end
%   
%   save('C:\Users\Anton\Documents\MATLAB\Brain_model_test_CSD\scout_scouts_Atlas_rand_symm_64_8000.mat', 'Scouts', 'TessNbVertices', 'Name');
    
      
   
end 



 








