function time=syncBreak(cell1,cell2)
% Function to determine the time at which synchrony between two cells breaks
% Inputs:
%   cell1 - Time series data for cell 1
%   cell2 - Time series data for cell 2
% Output:
%   time - The time index at which the synchrony breaks (NaN if no break found)

syncThreshold = 0.8; % Threshold correlation below which synchrony is considered broken

% Find peaks and troughs in cell1
[peaks, peakLocs] = findpeaks(cell1);
[troughs, troughLocs] = findpeaks(-cell1);

% Combine and sort the locations of peaks and troughs
criticalPoints = sort([peakLocs, troughLocs]);

time = NaN;
% Iterate through each pair of consecutive critical points
for i = 1:length(criticalPoints)-1
    startIdx = criticalPoints(i);
    endIdx = criticalPoints(i+1);

    segment1 = cell1(startIdx:endIdx);
    segment2 = cell2(startIdx:endIdx);

    syncScore = corr(segment1', segment2');

    % Check if the synchrony score is below the threshold
    if syncScore < syncThreshold
        time = startIdx;  % Return the starting index of the break
        return;
    end
end
end
