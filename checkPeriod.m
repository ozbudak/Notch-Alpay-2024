function [period] = checkPeriod(time, sol)
% Function to calculate the period of oscillations in a given solution
% Inputs:
%   time - Time series data
%   sol - Solution data corresponding to the time series
% Output:
%   period - Calculated period of the oscillations

% Initialize variables
index_peaks=zeros(2,1);
peaks=0;
sum=0;

% Find peaks
for z=2:length(sol)-1
    if sol(z-1)<sol(z)&&sol(z)>sol(z+1)
        peaks=peaks+1;
        index_peaks(peaks)=time(z);
    end
end

% Calculate period
for m=3:peaks-1
    sum=sum+index_peaks(m+1)-index_peaks(m);
end
period=(sum)/(peaks-3);

% Handle NaN case
if isnan(period)
    period=0;
end
end