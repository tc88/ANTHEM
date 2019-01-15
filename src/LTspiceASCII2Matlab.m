function data = LTspiceASCII2Matlab(file,verbose)
% LTSPICEASCII2MATLAB extracts LTspice results from an ASCII result file
% and returns them as a Matlab struct.
%
% Input:
%   file        can be given with or without file extension.
%               Default file extension is '.raw'
%   verbose     triggers console outputs and plots
%               (optional, default: [1 0])
%
% Output:
%   data        output data in the format as provided by LTspice2Matlab by
%               Paul Wagner. Additionally, if a log file is available,
%               information from the log file is added
%
% See also importSPICEresults, LTspice2Matlab, Xyce2Matlab
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if verbose(1)
    tstart = tic;
    fprintf('extracting LTspice ASCII results ...\n');
end

% if no file extension is given, use '.raw'
[filepath,filename,ext] = fileparts(file);
if isempty(ext), file = [file,'.raw']; end

% read content from ASCII result file
fid = fopen(file);
content = textscan(fid,'%s','delimiter','\n');
fclose(fid);
clear fid

% fill variable names and result matrix
variableNames = {};
mat = [];
colIdx = 1;
rowIdx = 1;
for i = 1:length(content{1})
    line = content{1}{i};
    fields = strsplit(line);
    % find variable names
    if numel(fields) == 3 && strcmp(fields{3},'voltage')
        var = fields{2};
        variableNames{end+1} = var;
    end
    % new data point, save its value as column index and go back to first row
    if numel(fields) == 2 && all(~isnan(str2double(fields)))
        colIdx = str2double(fields{1})+1;
        rowIdx = 1;
        mat(rowIdx,colIdx) = str2double(fields{2});
        rowIdx = rowIdx+1;
    % save single data entries to matrix
    elseif numel(fields) == 1 && ~isnan(str2double(fields))
        mat(rowIdx,colIdx) = str2double(fields{1});
        rowIdx = rowIdx+1;
    end
end
clear content

% extract information from log file, if available
logdata = struct();
logfile = [filepath filename '.log'];
if exist(logfile,'file')
    fid = fopen(logfile);
    content = textscan(fid,'%s','delimiter','\n');
    fclose(fid);
    for i = 1:length(content{1})
        line = content{1}{i};
        fields = strsplit(line);

        switch fields{1}
            case 'tnom',        logdata.tnom        = str2double(fields{3});
            case 'temp',        logdata.temp        = str2double(fields{3});
            case 'method',      logdata.method      =            fields{3};
            case 'totiter',     logdata.totiter     = str2double(fields{3});
            case 'traniter',    logdata.traniter    = str2double(fields{3});
            case 'tranpoints',  logdata.tranpoints  = str2double(fields{3});
            case 'accept',      logdata.accept      = str2double(fields{3});
            case 'rejected',    logdata.rejected    = str2double(fields{3});
            case 'fillins',     logdata.fillins     = str2double(fields{3});
            case 'solver',      logdata.solver      =            fields{3};
            case 'matrix'
                if strcmp(fields{2},'size'), logdata.matrixSize = str2double(fields{4}); end
            case 'Total'
                if strcmp([fields{2} fields{3}],'elapsedtime:'), logdata.time = str2double(fields{4}); end
        end
    end
end

% output data
data.variable_name_list = variableNames';
data.variable_mat = mat;
data.num_data_pnts = size(mat,2);
data.logdata = logdata;

if verbose(1)
    fprintf('finished extracting LTspice ASCII results after %d seconds.\n',toc(tstart));
end

end