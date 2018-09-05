function refinedAxis = refineAxis(axisData,factor)
% REFINEAXIS refines a given axis by a specified factor.
%
% Input:
%   axis_data     vector that specifies the axis to be refined
%   factor        factor by which the axis shall be refined
%                 (optional, default: 2)
% 
% Output:
% refinedAxis   refined axis by factor
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 2, factor = 2; end

if factor == 1
    refinedAxis = axisData;
elseif factor == 2
    intervalCenters = 0.5*(axisData(1:end-1)+axisData(2:end));
    refinedAxis = sort(union(axisData,intervalCenters));
else
    % collects all points that are between the original points
    intermediatePoints = [];
    for i = 1:(factor-1)
        intermediatePoints = [intermediatePoints,axisData(1:end-1)+i/factor*(axisData(2:end)-axisData(1:end-1))];
    end
    % unifies the original axis with the additional points and sorts
    refinedAxis = sort(union(axisData,intermediatePoints));
end

end