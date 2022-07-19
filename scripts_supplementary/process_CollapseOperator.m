% University of Helsinki, BABA Center
% Anton Tokariev, Sami Auno, 2019
%
% Calculate weights for source downsampling to parcels
% Requirements:
% You need to calculate the head model first in Brainstorm. Use openMEEG
% Also calculate sources for one measurement to get inverse model.

% Dummy signal generation settings
fs = 250;   % Hz
f = 10.4167;
L = 15000; % samples
iter = 100; % N of iterations

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


% Remove to unused channels (BIPs and ECG and stuff) from the forward
% operator
HeadModel.Gain = HeadModel.Gain(channelIndecis',:);

%% Make the head model constrained
ForwardModel = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);

%% Load the surface file
SurfaceMat = in_tess_bst(kernel.SurfaceFile);

iAtlas = SurfaceMat.iAtlas;     % NOTE! This assumes that the correct parcellation has been selected! Check
Parcellation = SurfaceMat.Atlas(iAtlas).Scouts;
%% Calculate weights iter times

% init
N_parcels = length(Parcellation);
N_sources = size(ForwardModel, 2);
CollapseOperators = zeros(N_sources, iter); % Each iteration of collapse operators (array of weights) are saved here

% Start iterating
for N_iterations = 1:iter
    
    % (1)  
    % Generate PARCEL signals
    parcel_sig_sim = get_signals(N_parcels, L, f, fs); % [parcels x samples]
    parcel_sig_sim = real(parcel_sig_sim);
    
    
    % (2)
    % Set PARCEL signals to SRC
    source_sig_sim = zeros(N_sources, L);
    for j = 1:N_parcels
        source_sig_sim(Parcellation(j).Vertices,:) = repmat(parcel_sig_sim(j,:),[length(Parcellation(j).Vertices),1]);
    end

    % (3)  
    % Compute EEG (Forward Modelling)
    EEG = ForwardModel * source_sig_sim;

    % (4)  
    % Compute SRC from EEG (Inverse Modelling)
    source_mod = InverseModel * EEG;
    
    % (5) >> Collapse Operator
    % Compute correlation: SRC modeled  vs. initial signal >> src weights
    for j = 1:N_parcels
        
        parcel_src_mod = source_mod(Parcellation(j).Vertices, :);
        
        buf = zeros(1, size(parcel_src_mod, 1));

        for k = 1:size(parcel_src_mod, 1)

          buf(k) = corr(parcel_src_mod(k, :)', parcel_sig_sim(j, :)');

        end

        CollapseOperators(Parcellation(j).Vertices, N_iterations) = buf;

    end

  
 
  
  
end % N iterations


CollapseOperator = mean(CollapseOperators,2);
std_CollapseOperator = mean(std(CollapseOperators,0,2));



%% Save the CollapseOperator to HeadModel


HeadModel.CollapseOperator = CollapseOperator;

save(path_headModel,'-struct','HeadModel');



















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
