function [isTrue] = checkSusOsc(time, sol)
% Function to check for sustained oscillations in a given solution
% Inputs:
%   time - Time series
%   sol - Solution data corresponding to the time series
% Output:
%   isTrue - Boolean indicating whether sustained oscillations are present

% Avoid checking the stabilization state
quaterpoint = ceil(length(time) / 4);
time = time(quaterpoint:end);
sol = sol(quaterpoint:end);

% Initialize variables
peaks=0;
troughs=0;
index_peaks=zeros(2,1);
index_troughs=zeros(2,1);
value_peaks=zeros(2,1);
value_troughs=zeros(2,1);

% Find peaks and troughs
for z=2:length(sol)-1
    if sol(z-1)<sol(z)&&sol(z)>sol(z+1)
        peaks=peaks+1;
        index_peaks(peaks)=time(z);
        value_peaks(peaks)=sol(z);
    end
    if sol(z-1)>sol(z)&&sol(z)<sol(z+1)
        troughs=troughs+1;
        index_troughs(troughs)=time(z);
        value_troughs(troughs)=sol(z);
    end
end

% Check conditions for sustained oscillations
if peaks>=3
    half_peaks=index_peaks(index_peaks<round(length(time)/2));
    mid=size(half_peaks,1);
    if mid == 0
        mid=1;
    end
    for k=1:length(value_peaks)
        if value_peaks(k)<0.1 %Should never have peak w/ amplitude less than 0.1 (indicates noise in the system).
            isTrue=false;
            return
        end
    end
    a=find(index_peaks>=round(length(time)*0.75), 1);
    b=find(index_troughs>=round(length(time)*0.75),1);
    if isempty(a) || isempty(b) %Should still have oscillations 3/4 of the way through data.
        isTrue=false;
        return
    end
    if peaks<troughs %Case 1: Fewer peaks than troughs
        peakToTrough=zeros(peaks,1);
        for z=1:peaks
            peakToTrough(z)=value_peaks(z)/((value_troughs(z)+value_troughs(z+1))/2);
        end
        ratioEnd = peakToTrough(peaks);%peak-to-trough ratio at end of run.
        if mid<=length(peakToTrough)
            ratioMid = peakToTrough(mid); %peak-to-trough ratio at middle of run.
        else
            ratioMid = peakToTrough(peaks);
        end
    elseif peaks>troughs %Case 2: More peaks than troughs.
        peakToTrough = zeros(troughs,1);
        for z=1:troughs
            peakToTrough(z)=((value_peaks(z)+value_peaks(z+1))/2)/value_troughs(z);
        end
        ratioEnd = peakToTrough(troughs);
        if mid<=length(peakToTrough)
            ratioMid = peakToTrough(mid);
        else
            ratioMid = peakToTrough(troughs);
        end
    elseif peaks==troughs %Case 3: Same number of peaks as troughs.
        peakToTrough = zeros(peaks,1);
        for z=1:troughs
            peakToTrough(z)=value_peaks(z)/(value_troughs(z));
        end
        ratioEnd = peakToTrough(troughs);
        if mid<=length(peakToTrough)
            ratioMid = peakToTrough(mid);
        else
            ratioMid = peakToTrough(troughs);
        end
    end
    amplitudeMid=(value_peaks(mid)-value_troughs(mid))/2;%Calculate amplitude at middle point.
    if ratioEnd >= 1.5 && ratioMid >= 1.5 && ratioMid/ratioEnd <= 1.02 && ratioMid/ratioEnd>=0.98 && amplitudeMid>=0.15 %Apply all conditions.
        isTrue=true;
        return;
    else
        isTrue=false;
        return;
    end
else
    isTrue=false;
    return;
end
end