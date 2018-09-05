function [timeClean,potClean,idxTime2keep] = spiceCleanInitRamp(time,pot,potInit)
% SPICECLEANINITRAMP cleans initial ramping in time from the provided
% potential data 'pot'. The time of the initial ramping is defined as the
% time the potential 'pot' needs to reach 'potInit'.
%
% Input:
%   time     time axis (1-by-nt)
%   pot      potential data (np-by-nt)
%   potInit  initial value of potential (scalar)
%
% Output:
%   timeClean     time axis cleaned from ramping (1-by-ntClean)
%   potClean      potential data cleaned from ramping (np-by-ntClean)
%   idxTime2keep  indices for time axis that are kept (1-by-ntClean)
%
% See also LTspice2Matlab, runLTspice, spiceOrderOutput
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% extract number of points in space (np) and time (nt)
nt = size(pot,2);

% find time indices that shall be removed
idx2remove = find(pot < potInit);
[~,idxTime2remove] = ind2sub(size(pot),idx2remove);
idxTime2remove = unique(idxTime2remove);

% redefines time and potential according to the indices that are kept
idxTime2keep = setdiff(1:nt,idxTime2remove);
timeClean = time(:,idxTime2keep);
potClean = pot(:,idxTime2keep);

end