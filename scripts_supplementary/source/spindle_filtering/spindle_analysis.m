

starts = OUTPUT.start;
ends = OUTPUT.end;

N_chans = size(starts,1);

L = length(eeg.signal);

I = zeros(N_chans,L,3);

for c = 1:N_chans
    
    for s = 1:sum(~isnan(OUTPUT.start(c,:)))
        
        indx = starts(c,s):ends(c,s);
        
        I(c,indx,1) = 1;
        
    end
    
    
end

figure(1)
clf
hold on
imagesc(I)


%% Calculate precision
av_sp = (OUTPUT.end + OUTPUT.start)/2;
av_sp = sort(av_sp(:));
av_sp = av_sp(~isnan(av_sp));

% go through each event. calculate how many spindles are within events
TP = 0;
for evi = 1:size(eeg.spindleEvents,2)
    TP = TP + length(find((av_sp >= eeg.spindleEvents(1,evi)) & (av_sp <= eeg.spindleEvents(2,evi))));
end

FP = length(av_sp) - TP;

precision = TP/(TP + FP)

%% Spindle density
av_sp = (OUTPUT.end + OUTPUT.start)/2;
av_sp = sort(av_sp(:));
av_sp = av_sp(~isnan(av_sp));

% divide into 1 min epochs, calculate number of spindles within each epoch

epoch_s = 20;
epoch_olp = 1;
epoch_inds = round(epoch_s * eeg.srate);

N_epochs = length(eeg.signal)/epoch_inds;

dens = [];
tims = [];

st = 1;
while 1
    en = st + epoch_inds;
    
    N_spindls = sum(av_sp > st & av_sp < en);
    
    dens = [dens;N_spindls];
    
    
    if en >= length(eeg.signal)
        tims = [tims; eeg.times(round((st+length(eeg.signal))/2)) ];
        break
    end
    
    tims = [tims; eeg.times(round((st+en)/2)) ];
    
    st = en + 1;
    st = 1 + st + round(epoch_inds*(1 - epoch_olp));
    
    if st >= length(eeg.signal)
        break
    end
    
end

dens = dens/(epoch_s / 60);

% Plot
figure(1)
clf
plot(tims/60,dens,'k')

xlabel('Time (min)')
ylabel('Spindle density (N/min)')




%% Group spindles
spin_start = OUTPUT.start(:);
spin_end = OUTPUT.end(:);

spin_lngth = spin_end - spin_start;

spin_start = spin_start(~isnan(spin_lngth));
spin_end = spin_end(~isnan(spin_lngth));
spin_lngth = spin_lngth(~isnan(spin_lngth));

[spin_start,sind] = sort(spin_start);
spin_end = spin_end(sind);
spin_lngth = spin_lngth(sind);

gi = 1;

spindleGroup = [];

while 1
    
	gr = find_group_indicis(spin_start,spin_end,1);
    
    spindleGroup(gi).start = min(spin_start(gr));
    spindleGroup(gi).end = max(spin_end(gr));
    spindleGroup(gi).median = median(((spin_start(gr) + spin_end(gr)))/2,'all');
    spindleGroup(gi).mean = mean(((spin_start(gr) + spin_end(gr)))/2,'all');
    spindleGroup(gi).std = std(spin_end(gr) - spin_start(gr),[],'all');
    
    spindleGroup(gi).time = eeg.times(round(mean(((spin_start(gr) + spin_end(gr)))/2,'all')));
    spindleGroup(gi).N_spindles = length(gr);
    
    % Remove those spindles that are already grouped
    spin_start = spin_start((gr(end)+1):end,1);
    spin_end = spin_end((gr(end)+1):end,1);
    
    gi = gi +1;
    
    if isempty(spin_start)
        break
    end
    
end

% fill a vector with 0 where there are no spindle groups, and with the
% number of spindles where there are groups
nSpinArr = zeros(1,length(eeg.times));

for i = 1:length(spindleGroup)
    nSpinArr(1,spindleGroup(i).start:spindleGroup(i).end) = spindleGroup(i).N_spindles;
end

% Plot
figure(1)
clf
plot(eeg.times , nSpinArr,'k')

% go through each event. calculate how many spindles are within events
TP = 0;

av_sp = [spindleGroup.median];
N_sp = [spindleGroup.N_spindles];
thr = mean(N_sp) + std(N_sp);
av_sp = av_sp(N_sp>thr);


for evi = 1:size(eeg.spindleEvents,2)
    TP = TP + length(find((av_sp >= eeg.spindleEvents(1,evi)) & (av_sp <= eeg.spindleEvents(2,evi))));
end

FP = length(av_sp) - TP;

precision = TP/(TP + FP)


% fill a vector with 0 where there are no spindle groups, and with the
% number of spindles where there are groups
nSpinArr = zeros(1,length(eeg.times));

inds = 1:length(spindleGroup);
inds = inds(N_sp>thr);
for i = inds
    nSpinArr(1,spindleGroup(i).start:spindleGroup(i).end) = spindleGroup(i).N_spindles;
end

% Plot
figure(2)
clf
plot(eeg.times , nSpinArr,'k')



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

































