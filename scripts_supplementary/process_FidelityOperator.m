% University of Helsinki, BABA Center
% Anton Tokariev, Sami Auno, 2019
%
% EPI_Infant project
% Calculate adjacency correlation matrix
% Requirements:
% (1) You need to calculate the head model first in Brainstorm. Use openMEEG
% (2) Also calculate sources for one measurement to get inverse model.
% (3) The Collapse Operator has to be calculated before and saved into the
% headmodel (as is done in step4

% Dummy signal generation settings
options.fs = 250;                               % Hz
options.f = 10.4167;                            % So that lag is 6.00
options.L = 15000;                              % samples, 60 s = 15k/250 Hz
options.iter = 200;                             % N of iterations
options.lag = round(options.fs/(4*options.f));  % 90 degrees

% Paths to files. These are brainstorm generated files that are named as:
% headmodel_[yymmdd]_[hhmm].mat
% results_[your inverse method eg. dSPM]_EEG_KERNEL_[yymmdd]_[hhmm].mat

path_headModel = '/modeling_data/headmodel.mat';
path_inverseModel = '/modeling_data/dSPM_EEG_KERNEL.mat';


%% Load data
HeadModel = load(path_headModel);
kernel = load(path_inverseModel);
InverseModel = kernel.ImagingKernel;
channelIndecis = kernel.GoodChannel;

%% Load the surface file
SurfaceMat = in_tess_bst(kernel.SurfaceFile);

iAtlas = SurfaceMat.iAtlas;     % NOTE! This assumes that the correct parcellation has been selected! Check
Parcellation = SurfaceMat.Atlas(iAtlas).Scouts;

%% Calculate Random adjacency matric

wpli_rand = adjacency_random(InverseModel,HeadModel,Parcellation,options);

% Save adjacency correlation matrix to file
% to the same folder as the head model
% save(fullfile(bst_fileparts(path_headModel),'wpli_random.mat'),'wpli_rand' );


%% Calculate coupled adjacency matrices
% Note, this takes significant amount of time. Depending on your computer,
% might take weeks.

wpli = adjacency_coupled(InverseModel,HeadModel,Parcellation,options);

% Save adjacency correlation matrix to file
% to the same folder as the head model
save(fullfile(bst_fileparts(path_headModel),'wpli.mat'),'wpli' );


%% Calculate the Fidelity Operator, i.e. threshold the adjacency matrices

[FidelityOperator,wpli_bad] = calculateFidelityOperator(wpli,wpli_rand);

% save to HeadModel
HeadModel.FidelityOperator = FidelityOperator;
save(path_headModel,'-struct','HeadModel');


















%% Functions

function wpli_rand = adjacency_random(varargin)

InverseModel = varargin{1};
HeadModel = varargin{2};
Parcellation = varargin{3};
options = varargin{4};

%% Make the head model constrained
ForwardModel = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);

% Number of parcels
N_parcels = length(Parcellation);

%% Calculate random adjacency matrices

% Make vertex_parcel_index vector.
% It is a vector of N_sources_all long
% And contains the parcel index of each source
N_sources_all = size(ForwardModel, 2);
vertex_parcel_index = zeros(N_sources_all,1);
for j = 1:N_parcels
    vertex_parcel_index(Parcellation(j).Vertices,1) = j;
end
% Some sources don't belong to any parcels.
% This indecis_of_nonzero_vertices vector
% has the index of those sources that belong to a parcel.
% This is used later to loop over only those sources that belong to a
% parcel.
indecis_of_nonzero_vertices = find(vertex_parcel_index);
  
% NO SYNCHRONY  
wpli_rand = zeros(N_parcels, N_parcels, options.iter);

tic
for iter_j = 1:options.iter


    disp(['Iteration # ', num2str(iter_j) , '/', num2str(options.iter)]);

    % (1)
    % Signals
    parcel_sig_sim = real(get_signals(N_parcels, options.L , options.f, options.fs));      % [parcels x samples] REAL!
    
    
    % (2)
    % Set PARCEL signals to SRC
    source_sig_sim = zeros(N_sources_all, options.L);
    for j = indecis_of_nonzero_vertices
        source_sig_sim(j,:) = parcel_sig_sim(vertex_parcel_index(j),:);
    end
    
    
    % (3)  
    % Compute EEG (Forward Modelling)
    EEG = ForwardModel * source_sig_sim;

    % (4)  
    % Compute SRC from EEG (Inverse Modelling)
    source_mod = InverseModel * EEG; 

    % (5)
    % 'weight' src >> multiply by Collapse Operator
    source_mod_weighted = HeadModel.CollapseOperator.*source_mod;

    % (6)
    % Reconstruct parcel signals from weighted src signals
    parcel_src_mod_wei = zeros(N_parcels, options.L);

    for j = 1:N_parcels
        parcel_src_mod_wei(j, :) = mean(source_mod_weighted(Parcellation(j).Vertices, :));
    end

    for parcelA = 1:N_parcels
        for parcelB = parcelA + 1:N_parcels

            wpli_rand(parcelA, parcelB, iter_j) = get_wPLI(parcel_src_mod_wei(parcelA, :), parcel_src_mod_wei(parcelB, :), 1);

        end
    end

end % end iter

toc
end

function wpli = adjacency_coupled(varargin)

InverseModel = varargin{1};
HeadModel = varargin{2};
Parcellation = varargin{3};
options = varargin{4};

%% Make the head model constrained
ForwardModel = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);

% Number of parcels
N_parcels = length(Parcellation);

%% Calculate adjacency matrices for perfectly coupled system

N_sources_all = size(ForwardModel, 2);

vertex_index = zeros(N_sources_all,1);

for j = 1:N_parcels
    vertex_index(Parcellation(j).Vertices,1) = j;
end
 
wpli = zeros(N_parcels, N_parcels, options.iter);

N_calls = N_parcels*(N_parcels+1)/2;

for iter_j = 1:options.iter

    calls_i = 0;
    
    signal_ref_all_calls = real(get_signals(N_calls, options.L + options.lag, options.f, options.fs));
    signal_all_all_calls = real(get_signals(N_calls*N_parcels, options.L , options.f, options.fs));
    

    disp(['Iteration # ', num2str(iter_j) , '/', num2str(options.iter)]);
    tic
    for parcelA = 1:N_parcels
        for parcelB = parcelA + 1:N_parcels
            
            calls_i = calls_i + 1;
            
            % (1)
            % Signals
            signal_ref = signal_ref_all_calls(calls_i,:);
            signal_lag = signal_ref(1, options.lag+1:end);
            signal_ref(end-options.lag+1:end) = [];
            
            signal_all = signal_all_all_calls( ((calls_i - 1)*N_parcels + 1) : (calls_i*N_parcels) ,:);
            
            % integrate synchroneous agents 
            signal_all(parcelA, :) = signal_ref;
            signal_all(parcelB, :) = signal_lag;
            
            % (2)
            % Set PARCEL signals to SRC
            source_sig_sim = zeros(N_sources_all, options.L);
            for j = find(vertex_index)
                source_sig_sim(j,:) = signal_all(vertex_index(j),:);
            end

            % (3)  
            % Compute EEG (Forward Modelling)
            EEG = ForwardModel * source_sig_sim;

            % (4)  
            % Compute SRC from EEG (Inverse Modelling)
            source_mod = InverseModel * EEG; 

            % (5)
            % 'weight' src >> multiply by Collapse Operator
            source_mod_weighted = HeadModel.CollapseOperator.*source_mod;

            % (6)
            % Reconstruct parcel signals from weighted src signals
            parcel_src_mod_wei = zeros(N_parcels, options.L);

            for j = 1:N_parcels
                parcel_src_mod_wei(j, :) = mean(source_mod_weighted(Parcellation(j).Vertices, :));
            end


            wpli(parcelA, parcelB, iter_j) = get_wPLI(parcel_src_mod_wei(parcelA, :), parcel_src_mod_wei(parcelB, :), 1);

        end
    end
    t_time = toc;
    disp(['Done in ', num2str(t_time/60), ' minutes.']);
    disp(' ');
end % end iter

end

function [FidelityOperator,wpli_bad] = calculateFidelityOperator(wpli,wpli_rand)
%% Threshold the adjacency matrix and calculate the Fidelity Operator



wpli_rand_all = [];

for k = 1:size(wpli_rand, 3)

    wpli_rand_all = [wpli_rand_all nonzeros(triu(squeeze(wpli_rand(:, :, k)),1))'];

end

wpli_rand_all = abs(wpli_rand_all);


thr = prctile(wpli_rand_all, 95);

wpli = median(abs(wpli), 3); % abs




wpli_good = double(wpli > thr);

wpli_bad = double(wpli < thr & wpli > 0); 

FidelityOperator = wpli_good;

%     topoplot_3D(wpli_bad, 'Bad edges (top view)');

end


%% TODO
function h = topoplot_3D(value, titl)

% Load data
pathFlat = '/media/sami/ADATA SD700/EPI_INFANT_BRAINSTORM/brainstorm_db/EPI_Infant/anat/@default_subject';
flat_cx = load(pathFlat);

path_headModel = '/media/sami/ADATA SD700/EPI_INFANT_BRAINSTORM/brainstorm_db/EPI_Infant/data/EPI11/@rawEPI_11_1v_Lilja_Alexandersson_2017-12-07_12-38-01/headmodel_191128_1343.mat';
HeadModel = load(path_headModel);


Np = size(MyAtlas.Centroids, 1);

h = figure;
hold on

set(gcf, 'Color', 'w');



patch('Vertices', flat_cx.Vertices, 'Faces', flat_cx.Faces, 'FaceColor', [0.8 0.8 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');


scatter3(MyAtlas.Centroids(:, 1), MyAtlas.Centroids(:, 2), MyAtlas.Centroids(:, 3), 30, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');



for j = 1:Np
text(1.1*MyAtlas.Centroids(j, 1), 1.1*MyAtlas.Centroids(j, 2), 1.1*MyAtlas.Centroids(j, 3), num2str(j), 'FontSize', 8, 'Color', 'b', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end


for chA = 1:Np
for chB = (chA + 1):Np

if value(chA, chB) > 0

line([MyAtlas.Centroids(chA, 1) MyAtlas.Centroids(chB, 1)], [MyAtlas.Centroids(chA, 2) MyAtlas.Centroids(chB, 2)], [MyAtlas.Centroids(chA, 3) MyAtlas.Centroids(chB, 3)], 'Color', 'r', 'Linewidth', 0.5*value(chA, chB));

end

end
end

title(titl, 'FontSize', 16, 'FontWeight', 'bold');

axis off
hold off

%%% saveas(gcf, [titl '.png']);
  
   

end




function X = get_signals(nChannels, N, f, fs)

    % constants
    m = 5; % wavelet parameter

    % rand samples
    data = randn(nChannels, N);

    % mixing
    % A = rand(nChannels, nChannels);
    % data = A * data;

    % filtering
    X = complex(zeros(nChannels, N), zeros(nChannels, N));

    for nChannel = 1:nChannels
        X(nChannel, :) = fft(data(nChannel, :));
    end

    % wavelet init
    [a, b] = support_wavelet_init(f, m, fs);
    W = support_wavelet_fft(a, b, N);

    % FOR each Channel
    for nChannel = 1:nChannels
        % wavelet transform
        X(nChannel, :) = ifft(X(nChannel, :) .* W);
    end

end % end





function [p11, p12] = support_wavelet_init(f0, m, fs)

    % init
    p1 = 1 / fs;
    p2 = 1 / (f0 / m * (2 * pi));
    p3 = round((p2 * 10) / p1);
    p3 = p3 + rem(p3, 2);

    % carrier
    p4 = ((0:(p3 - 1)) * p1) - (((p3 - 1) * p1) / 2);

    % formulae
    p5 = (-1) * ((p4 .^ 2) / ((p2 .^ 2) * 2));
    p6 = p4 * (2 * pi * f0);

    % shape
    p7 = complex(p5, p6);
    p8 = exp(p7) * (1 / sqrt(p2 * sqrt(pi)));
    p9 = p8 / (sum(abs(p8)) / 2);

    % split into halves
    p10 = round(p3 / 2);
    p11 = p9(1:p10);
    p12 = p9((p10 + 1):end);

end % end





function p15 = support_wavelet_fft(p11, p12, N)

    % halves into complex form
    p13 = N - (length(p11) + length(p12));
    p14 = complex(zeros(1, p13), zeros(1, p13));
    p15 = [p12, p14, p11];

    % fft
    p15 = fft(p15);

end % end





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
