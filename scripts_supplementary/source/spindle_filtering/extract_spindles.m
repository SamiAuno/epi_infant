
function [EEG, OUTPUT] = extract_spindles(EEG)

sensors = 1:size(EEG.signal, 1); % ID of channels to be filtered

% Spindles filter

Wp = [11 15]/(EEG.srate/2); % to change in relation to the definition of spindles
Ws = [10 16]/(EEG.srate/2);
Rp = 3;
Rs = 60;
ftype = 'bandpass';

[n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);
    
[z,p,k] = cheby2(n,Rs,Ws,ftype);
sos = zp2sos(z,p,k);

h = waitbar(0,['EEG Performing spindle Filter...']);

for kk = 1:size(sensors, 2)
    EEG.spindles(sensors(kk), :) = sosfilt(sos,EEG.signal(sensors(kk), :));
end

% % Plot EEG
% eegplot(EEG.spindles, 'srate', EEG.srate, 'winlength', 30,'eloc_file', EEG.derivations)

for kk = 1:size(sensors, 2)
    h=waitbar(kk/size(sensors, 2));
    EEG.x(sensors(kk), :) = (hilbert(EEG.spindles(sensors(kk), :)));
end
close (h)

EEG.mag=abs(EEG.x);

% set spindle tresholds
for kk = 1:size(sensors, 2)
    M(sensors(kk),:) = mean(EEG.mag(sensors(kk),:));
    SD(sensors(kk),:) = std(EEG.mag(sensors(kk),:));

%     SD3=3*(SD);
    SD3=3*(SD);
    M3=3*M;
    M6=6*M;
    MSD=M+SD;
    M3SD=M+SD3;
end

for kk = 1:size(sensors, 2)
    for km = 1:size(EEG.mag,2)
        if (EEG.mag(kk,km) >= M(kk,:))  && (EEG.mag(kk,km) < M3(kk,:))   % to change in relation to the method adopted 
            EEG.on(kk,km)=1;
        
        elseif (EEG.mag(kk,km) >= M3(kk,:))
            EEG.on(kk,km)=3;
        
        else EEG.on(kk,km)=0;
        
        end 
   
    end
end 
    
EEG.on2 = EEG.on; 
EEG.on2(EEG.on<2)=NaN;

s = diff(EEG.on, 1, 2); % define spindle epoch - if M<MAG<M3 --> -1; if M3<MAG --> -2

% extract epochs
for kk = 1:size(sensors, 2)
    en = find(s(kk,:) == -1);
    EN2 = find(s(kk,:) == -2);
    st = find(s(kk,:) == 1); 
    ST2 = find(s(kk,:) == 2); 
    
    % Just to make sure that size(ST2) = size(EN2)
    if length(ST2) > length(EN2)
        ST2 = ST2(1:length(EN2));
    elseif length(ST2) < length(EN2)
        EN2 = EN2(1:length(ST2));
    end
    
    dd=(EN2-ST2');
    dd(dd<0)=nan;
    [minDD, DD]=min(dd,[],2);
    for y=1:size(minDD,1)
        if isnan(minDD(y))
            DD(y)=NaN;
        end  
    end

    ind1=isnan(DD); DD(ind1)=[];
    en2 = unique(EN2(DD));
    try
        st2 = unique(ST2(DD));
    catch
        disp('here')
    end

    d=(st2-st');
    d(d<0)=nan;
    [minST,closestInST] = nanmin(d);
    for y=1:size(minST,2)
        if isnan(minST(y))
            closestInST(y)=NaN;
        end  
    end

    c=(en2-en');
    c(c>0)=nan;
    [minEN,closestInEN] = max(c);
    for y=1:size(minEN,2)
        if isnan(minEN(y))
            closestInEN(y)=NaN;
        end  
    end

    ind1=isnan(closestInEN); closestInEN(ind1)=[]; closestInST(ind1)=[];
    ind2=isnan(closestInST); closestInST(ind2)=[]; closestInEN(ind2)=[];
    closestST=st(closestInST); 
    closestEN=en(closestInEN); 

    for f=1:size(closestST,2)
        FF=(closestEN(f)-closestST(f));
        if FF<EEG.srate*0.3
            closestST(f)=nan;
            closestEN(f)=nan;
            %        elseif FF>EEG.srate*100
            %            closestST(f)=nan;
            %            closestEN(f)=nan;
        end
    end
    ind1=isnan(closestEN); closestEN(ind1)=[];
    ind2=isnan(closestST); closestST(ind2)=[];

    for h=1:(size(closestST,2)-1)
        i=h+1;
        HH=(closestST(i)-closestEN(h));
        if HH<EEG.srate
            closestST(i)=nan;
            closestEN(h)=nan;
        end
    end
    ind3=isnan(closestEN); closestEN(ind3)=[];
    ind4=isnan(closestST); closestST(ind4)=[];

    for f=1:size(closestST,2)
        FF=(closestEN(f)-closestST(f));
        if FF<EEG.srate*0.5
            closestST(f)=nan;
            closestEN(f)=nan;
        elseif FF>EEG.srate*10
            closestST(f)=nan;
            closestEN(f)=nan;
        end
    end
    ind5=isnan(closestEN); closestEN(ind5)=[];
    ind6=isnan(closestST); closestST(ind6)=[];

    clear ind1 ind2 ind3 ind4 ind5 ind6

    numEpochs = size(closestST,2);

    for e = 1:(numEpochs)
        
        OUTPUT.start{kk,e}=(closestST(e));
        OUTPUT.end{kk,e}=(closestEN(e));
        A=closestST(e)-EEG.srate;
        
        if A<1 
            A(e)=closestST(e);
        else
            A(e)=A;
        end
        
        B=closestEN(e)+EEG.srate;
        
        if B>size(EEG.signal,2) 
            B(e)=closestEN(e);
        elseif B==1 
            B(e)=[];
        else
            B(e)=B;
        end

        OUTPUT.signal{kk,e}= EEG.signal(kk,(A(e):B(e)));
        OUTPUT.spindles{kk,e}= EEG.spindles(kk,(A(e):B(e)));
        OUTPUT.mag{kk,e}= EEG.mag(kk,(A(e):B(e)));
    end

end
%

fx=@(x)any(isempty(x));
inds=cellfun(fx,OUTPUT.start);
OUTPUT.start(inds)={nan};
inde=cellfun(fx,OUTPUT.end);
OUTPUT.end(inde)={nan};
OUTPUT.start=cell2mat(OUTPUT.start);
OUTPUT.end=cell2mat(OUTPUT.end);
fx=@(x)any(isempty(x));
indsi=cellfun(fx,OUTPUT.signal);
OUTPUT.signal(indsi)={nan};
indsp=cellfun(fx,OUTPUT.spindles);
OUTPUT.spindles(indsp)={nan};
indma=cellfun(fx,OUTPUT.mag);
OUTPUT.mag(indma)={nan};
clear indma indsp indsi inde fx fxOUTPUT 
end