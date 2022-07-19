%% Select spindelse epoch centering around one spindel group highlighting it
% 
%% NOT WORKING YET

% Divides the data into epochs
% Prints the epoch centering around one spindel group highlighting it
% 
%% NOT WORKING YET


function select_spindles(eeg,OUTPUT)

epoch_length_seconds = 10;          % in seconds
% epoch_overlap = 0.5;                  % ratio, 0.0 - 1.0

N_chans = size(eeg.spindles,1);

scf = 10^(round(log10(std(eeg.spindles(:))))+1); % scaling factor for plotting.


for g = length(OUTPUT.groups)
    
    figure(1)
    clf
    hold on
    
    x_time = eeg.times;

end


for s = epochs_index
    
    figure(1)
    clf
    hold on
    
    x_time = eeg.times( epoch_time_indicis{s} );
    
    % If there are spindle events, plot them
    if ~isempty(eeg.spindleEvents)
        for spEi = 1:size(eeg.spindleEvents,2)
            event_start = eeg.times(eeg.spindleEvents(1,spEi));
            event_stop = eeg.times(eeg.spindleEvents(2,spEi));
            if  x_time(1) <= event_start && event_stop <= x_time(end)
                
                rectangle('Position',[event_start,-100,event_stop-event_start,200],'FaceColor',[70,130,180]/255,'EdgeColor',[70,130,180]/255)
%                 plot([event_start event_start], [-100 100], '-b');
%                 plot([event_stop event_stop], [-100 100], '-b');
            
            end
        end
    end
    
    pci = 1;    % this is used for scaling in plotting. number of plotted channels.
    for ci = 1:N_chans


        orig_signal = eeg.spindles(ci,epoch_time_indicis{s});
        signal = orig_signal - median(orig_signal) - ((pci - 1)*scf);
        % First plot channel
        plot(x_time , signal,'-k');
        
        % Check if the channel has any spindles
        temp_spindles = find(ismember(OUTPUT.start(ci,:),epoch_time_indicis{s}));
        if temp_spindles

            % Then plot spindle if the channel has any in this epoch

            for sm = temp_spindles
                temp_spindle_indicis = OUTPUT.start(ci,sm) : OUTPUT.end(ci,sm);
                spindle_time = eeg.times( temp_spindle_indicis );

%                 OUTPUT.signal
                orig_spindle_signal = eeg.spindles(ci,temp_spindle_indicis);
                spindle_signal = orig_spindle_signal - median(orig_signal) - ((pci - 1)*scf);
                plot(spindle_time , spindle_signal , '-r');

            end
            
        end
        pci = pci + 1;
    end
    

    
    xlim([x_time(1) x_time(end)]);
    ylim([-pci*scf  scf]);
    
end










% Epoch index list.
epoch_length = round(epoch_length_seconds * eeg.srate);
i = 1;
temp_start_ind = 1;
epoch_time_indicis = [];    % this is a cell array. Each cell contains a vector of indices that correspond to the time vector
while 1
    
    temp_end_ind = temp_start_ind + epoch_length;
    
    if temp_end_ind >= size(eeg.signal,2)
        epoch_time_indicis{i,1} = temp_start_ind : size(eeg.signal,2);
        break
    end
    
    epoch_time_indicis{i,1} = temp_start_ind : temp_end_ind;
    
    temp_start_ind = 1 + temp_start_ind + round(epoch_length*(1 - epoch_overlap));
    i = i + 1;
    
end

N_epochs = length(epoch_time_indicis);


% OUTPUT = groupSpindles(OUTPUT,eeg);


% Get the index of those epochs (index to epoch_time_indicis) that has
% atleast one spindle starting in that epoch.
epochs_index = [];
for i=1:N_epochs
    
    % Check if the epoch contains any starting index
    if sum(ismember(epoch_time_indicis{i,1},OUTPUT.start(:))) > 0
        epochs_index = [epochs_index,i];
    end
    
end


% Plot each epoch and highlight the spindles

N_chans = size(eeg.spindles,1);

scf = 10^(round(log10(std(eeg.spindles(:))))+1); % scaling factor for plotting.

for s = epochs_index
    
    figure(1)
    clf
    hold on
    
    x_time = eeg.times( epoch_time_indicis{s} );
    
    % If there are spindle events, plot them
    if ~isempty(eeg.spindleEvents)
        for spEi = 1:size(eeg.spindleEvents,2)
            event_start = eeg.times(eeg.spindleEvents(1,spEi));
            event_stop = eeg.times(eeg.spindleEvents(2,spEi));
            if  x_time(1) <= event_start && event_stop <= x_time(end)
                
                rectangle('Position',[event_start,-100,event_stop-event_start,200],'FaceColor',[70,130,180]/255,'EdgeColor',[70,130,180]/255)
%                 plot([event_start event_start], [-100 100], '-b');
%                 plot([event_stop event_stop], [-100 100], '-b');
            
            end
        end
    end
    
    pci = 1;    % this is used for scaling in plotting. number of plotted channels.
    for ci = 1:N_chans


        orig_signal = eeg.spindles(ci,epoch_time_indicis{s});
        signal = orig_signal - median(orig_signal) - ((pci - 1)*scf);
        % First plot channel
        plot(x_time , signal,'-k');
        
        % Check if the channel has any spindles
        temp_spindles = find(ismember(OUTPUT.start(ci,:),epoch_time_indicis{s}));
        if temp_spindles

            % Then plot spindle if the channel has any in this epoch

            for sm = temp_spindles
                temp_spindle_indicis = OUTPUT.start(ci,sm) : OUTPUT.end(ci,sm);
                spindle_time = eeg.times( temp_spindle_indicis );

%                 OUTPUT.signal
                orig_spindle_signal = eeg.spindles(ci,temp_spindle_indicis);
                spindle_signal = orig_spindle_signal - median(orig_signal) - ((pci - 1)*scf);
                plot(spindle_time , spindle_signal , '-r');

            end
            
        end
        pci = pci + 1;
    end
    

    
    xlim([x_time(1) x_time(end)]);
    ylim([-pci*scf  scf]);
    
end




N_full_epochs = floor(size(eeg.signal,2) / epoch_length);
% length_remainder_epoch = mod(size(eeg.signal,2),epoch_length);
% if length_remainder_epoch > 0 % If remainder is nonzero
%     epoch_lengths = [epoch_length * ones(1,N_full_epochs) , length_remainder_epoch];
% else
%     epoch_lengths = epoch_length * ones(1,N_full_epochs);
% end
% temp_indx = 1 : size(eeg.signal,2); % just a running array. 
% 
% epoch_indicis = mat2cell(temp_indx,1,epoch_lengths);







% 
% sensors = 1:size(EEG.signal, 1); % ID of channels to be filtered
% 
% % plot spindles epoch in the different channels
% 
% for kk = 1:size(sensors, 2) 
%     nepoch = size(OUTPUT.start(~isnan(OUTPUT.start(kk,:))),2);
%     if nepoch==0
%         plot(EEG.times,(EEG.mag(kk,:))-(kk*100),'Color','c') 
%     else 
%         for ww = 1:nepoch
%             x=EEG.times(OUTPUT.start(kk,ww):OUTPUT.end(kk,ww));
%             plot(x,(EEG.spindles(OUTPUT.start(kk,ww):OUTPUT.end(kk,ww))-(kk*100)),'Color','r')
% %             xlim([1 size(EEG.times,2)])
%             hold on
%         end
%         plot(EEG.times,(EEG.mag(kk,:))-(kk*100),'Color','c')
%     end
%    hold on
% end
% hold off
% %
% % saveas(gcf,[PN sprintf('spindles_extraction%s',FileName)],'tiff')
% 
% close all

end