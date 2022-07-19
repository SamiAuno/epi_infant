%% Group spindles

function OUTPUT = groupSpindles(OUTPUT,eeg)



spin_start = OUTPUT.start(:);
spin_end = OUTPUT.end(:);

spin_lngth = spin_end - spin_start;

spin_start = spin_start(~isnan(spin_lngth));
spin_end = spin_end(~isnan(spin_lngth));
spin_lngth = spin_lngth(~isnan(spin_lngth));

[spin_start,sind] = sort(spin_start);
spin_end = spin_end(sind);
spin_lngth = spin_lngth(sind);







%% Remove spindles that extend over time discontinueties
% First find those spindles that extend over discontinueties in time (those
% that are around cut epochs. Then remove them.

% timeDiffs = diff(eeg.times);
% 
% disconts = [timeDiffs>(min(diff(eeg.times))*1.5)];
% disconts_ind = find(disconts);
% 
% % loop through discontinueties
% for d = 1:length(disconts_ind)
%     
%     maskDisc = spin_start < disconts_ind(d) & spin_end > disconts_ind(d);
%     
%     spin_start = spin_start(~maskDisc);
%     spin_end = spin_end(~maskDisc);
%     
%     
% end

gi = 1;

spindleGroup = [];

while 1
    
	gr = find_group_indicis(spin_start,spin_end,1);
    
    spindleGroup(gi).start = min(spin_start(gr));
    spindleGroup(gi).end = max(spin_end(gr));
    
    spindleGroup(gi).length_seconds = eeg.times(spindleGroup(gi).end) - eeg.times(spindleGroup(gi).start);
    spindleGroup(gi).N_spindles = length(gr);
    
    % Remove those spindles that are already grouped
    spin_start = spin_start((gr(end)+1):end,1);
    spin_end = spin_end((gr(end)+1):end,1);
    
    gi = gi +1;
    
    if isempty(spin_start)
        break
    end
    
end

% % Threshold spindle groups by the number of spindles.
% % The threshold is mean number + 1 std
% N_sp = [spindleGroup.N_spindles];
% thr = mean(N_sp) + std(N_sp);
% spindleGroup = spindleGroup(N_sp>thr);




% Remove spindles from OUTPUT that are not in groups

mask = false(size(OUTPUT.start,1),size(OUTPUT.start,2));

for i=1:size(spindleGroup,2)
    Starts = spindleGroup(i).start <= OUTPUT.start;
    Ends = spindleGroup(i).end >= OUTPUT.end;
    
    mask(Starts & Ends) = true;

end


OUTPUT = trimNpadOUTPUT(OUTPUT,mask);

OUTPUT.groups = spindleGroup;



% % fill a vector with 0 where there are no spindle groups, and with the
% % number of spindles where there are groups
% nSpinArr = zeros(1,length(eeg.times));
% 
% for i = 1:length(spindleGroup)
%     nSpinArr(1,spindleGroup(i).start:spindleGroup(i).end) = spindleGroup(i).N_spindles;
% end
% 
% % Plot
% figure(1)
% clf
% plot(eeg.times , nSpinArr,'k')
% 
% % go through each event. calculate how many spindles are within events
% TP = 0;
% 
% av_sp = [spindleGroup.median];
% N_sp = [spindleGroup.N_spindles];
% thr = mean(N_sp) + std(N_sp);
% av_sp = av_sp(N_sp>thr);
% 
% 
% % for evi = 1:size(eeg.spindleEvents,2)
% %     TP = TP + length(find((av_sp >= eeg.spindleEvents(1,evi)) & (av_sp <= eeg.spindleEvents(2,evi))));
% % end
% 
% % FP = length(av_sp) - TP;
% 
% % precision = TP/(TP + FP)
% 
% 
% % fill a vector with 0 where there are no spindle groups, and with the
% % number of spindles where there are groups
% % nSpinArr = zeros(1,length(eeg.times));
% % 
% % inds = 1:length(spindleGroup);
% % inds = inds(N_sp>thr);
% % for i = inds
% %     nSpinArr(1,spindleGroup(i).start:spindleGroup(i).end) = spindleGroup(i).N_spindles;
% % end







end

function gr = find_group_indicis(spin_start,spin_end,group_indx)

gr = [];

for i = group_indx
    
    % Take the first spindle in the group
    first_ind = [spin_start(i),spin_end(i)];
    
    % Get index of spindles that start before the first spindle ends
    new_gr = find( first_ind(2) > spin_start)';
    
    % check which of these are new (don't belong to the original group)
    gr = unique([gr,new_gr]);
    
end

if gr(~ismember(gr,group_indx))
    gr = find_group_indicis(spin_start,spin_end,gr);
end

gr = unique([gr,group_indx]);


end



function OUTPUT_new = trimNpadOUTPUT(OUTPUT,mask)

max_nSpindles = max(sum(mask,2));

% Get the OUTPUT field names 
fields = fieldnames(OUTPUT);

% Loop through channels
for k = 1:size(mask,1)
    % loop through fields
    for fl = 1:length(fields)
        temp = OUTPUT.(fields{fl})(k,mask(k,:)); % new row in the field
        
        % pad with nans
        padN = max_nSpindles-size(temp,2);
        if iscell(temp)
            temp = [temp, num2cell(nan(1,padN))];
        else
            temp = [temp, nan(1,padN)];
        end
        
        OUTPUT_new.(fields{fl})(k,:) = temp;

    end
    
end


end



