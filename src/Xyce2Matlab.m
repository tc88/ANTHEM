function data = Xyce2Matlab(filename,verbose)
% XYCE2MATLAB extracts Xyce results from ASCII result file and returns them
% as a Matlab struct.
%
% Input:
%   filename    can be given with or without file extension.
%               Default file extension is '.prn'
%   verbose     triggers console outputs and plots
%               (optional, default: [1 0])
%
% Output:
%   data        output data in the format as provided by LTspice2Matlab by
%               Paul Wagner
%
% See also importSPICEresults, LTspice2Matlab, LTspiceASCII2Matlab
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if verbose(1)
    tstart = tic;
    fprintf('extracting Xyce results ...\n');
end

% if no file extension is given, use '.raw'
[~,~,ext] = fileparts(filename);
if isempty(ext), filename = [filename,'.prn']; end

% read content from ASCII result file
fid = fopen(filename);
content = textscan(fid,'%s','delimiter','\n');
fclose(fid);

% find header and data matrix
XyceHeader = strsplit(content{1}{1});
dataSplit = cellfun(@strsplit,content{1},'UniformOutput',false);
mat = cell2mat(cellfun(@str2double,dataSplit(2:end-1),'UniformOutput',false))';

% last field of XyceHeader might be empty
if isempty(XyceHeader{end}), XyceHeader = XyceHeader(1:end-1); end

% output data
data.variable_name_list = XyceHeader';
data.variable_mat = mat;
data.num_data_pnts = size(mat,2);

if verbose(1)
    fprintf('finished extracting Xyce results after %d seconds.\n',toc(tstart));
end

end