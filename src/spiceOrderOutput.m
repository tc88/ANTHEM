function [u,idx] = spiceOrderOutput(data,prefix,suffix,verbose)
% SPICEORDEROUTPUT orders SPICE output data according to the node numbers
% such that it is consistent with the canonical numbering in FIT.
% Using prefix/suffix, a subset of all nodes can be selected.
%
% Input:
%   data    SPICE data as returned by LTspice2Matlab
%   prefix  prefix used before node number in node name
%           (optional,default:'')
%   suffix  suffix used after node number in node name
%           (optional,default:'')
%   verbose triggers console outputs and plots
%           (optional, default: [1 0])
%
% Output:
%   u       ordered result vector/matrix (np-by-1)
%   idx     vector of indices mapping SPICE order to node numbering order
%
% See also LTspice2Matlab, runLTspice, spiceCleanInitRamp
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 4, verbose = [1,0]; end
if nargin < 3, suffix = ''; end
if nargin < 2, prefix = ''; end

% input check
if ~isrow(data.variable_name_list), data.variable_name_list = data.variable_name_list'; end

tstart = tic;
if verbose(1), fprintf('reordering circuit results to be in accordance with canonical FIT numbering scheme ...\n'); end

% extract nodes according to prefix and suffix
idxListCell = num2cell(1:length(data.variable_name_list));
fun2call = @(nameListCell,idxListCell) findNodeNr(nameListCell,idxListCell,prefix,suffix);
[nodeList,idx] = cellfun(fun2call,data.variable_name_list,idxListCell);

% remove NaN entries
nodeList = nodeList(~isnan(nodeList))';
idx = idx(~isnan(idx))';

% fill result vector/matrix
u = zeros(length(nodeList),data.num_data_pnts);
u(nodeList,:) = data.variable_mat(idx,:);

if verbose(1), fprintf('finished reordering circuit results after %d seconds\n',toc(tstart)); end

function [currentNodeNr,outputIdx] = findNodeNr(currentEntry,currentIdx,prefix,suffix)
    currentNodeNr = NaN;
    outputIdx = NaN;
    % check for variable type: voltage
    if strcmp(currentEntry(1),'V')
        nodeName = currentEntry(3:end-1);
        % check whether nodeName is long enough to contain prefix and suffix
        if length(nodeName) < length(prefix) + length(suffix), return; end
        % find pre- and suffix in nodeName and check whether empty
        if isempty(prefix), idxPrefix = 1;
        else, idxPrefix = strfind(nodeName,prefix);
        end
        if isempty(suffix), idxSuffix = 2;
        else, idxSuffix = strfind(nodeName,suffix);
        end
        if isempty(idxPrefix), idxPrefix = 0; end
        if isempty(idxSuffix), idxSuffix = 0; end
        % if prefix and suffix evaluated successfully, add current nodeName to dictionary
        if idxPrefix && idxSuffix && idxSuffix > idxPrefix
            currentNodeNr = str2double(nodeName(length(prefix)+1:end-length(suffix)));
            if ~isnan(currentNodeNr), outputIdx = currentIdx; end
        end
    end
end

end