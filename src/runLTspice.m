function [status,cmdout] = runLTspice(file,options,LTspiceEXE,wineEXE)
% RUNLTSPICE runs the LTspice software to simulate a given netlist.
% The result is saved in the LTspice binary output file (extension .raw).
%
% Input:
%   file            file containing netlist. May be given with/without path
%                   and/or with/without file extension.
%   options         options that are passed as command line option to the
%                   LTspice binary. Most common option is '-b' for blind 
%                   execution, meaning that LTspice is executed in the
%                   background. (optional, default: '')
%   LTspiceEXE      path to LTspice executable
%                   (optional, default: ~/.wine/drive_c/Program Files (x86)/LTC/LTspiceIV/scad3.exe)
%   wineEXE         path to wine executable (only required on UNIX)
%                   (optional, default: '/usr/bin/wine')
%
% Output:
%   status          command output as returned by LTspice
%   cmdout          command output as returned by LTspice
%
% See also LTspice2Matlab, spiceCleanInitRamp, spiceOrderOutput
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 4, wineEXE = '/usr/bin/wine'; end
if nargin < 3, LTspiceEXE = '~/.wine/drive_c/Program Files (x86)/LTC/LTspiceIV/scad3.exe'; end
if nargin < 2, options = ''; end

% check whether path to wine and LTspice are valid
if ~exist(wineEXE,'file'),    error('did not find wine executable in given path (%s), please adjust ''wineEXE''',wineEXE); end
if ~exist(LTspiceEXE,'file'), error('did not find LTspice executable in given path (%s), please adjust ''LTspiceEXE''',LTspiceEXE); end

tstart = tic;
fprintf('running circuit simulation using LTspice ...\n');

[~,name,~] = fileparts(file);

% escapes unescaped whitespaces and parantheses in LTspiceEXE
idx2escape = find(isspace(LTspiceEXE));
idx2escape = union(idx2escape,strfind(LTspiceEXE,'('));
idx2escape = union(idx2escape,strfind(LTspiceEXE,')'));
for i = flip(idx2escape)
    if ~strcmp(LTspiceEXE(i-1),'\')
        LTspiceEXE = [LTspiceEXE(1:i-1),'\',LTspiceEXE(i:end)];
    end
end

% fixes issue that files containing 'B' cannot be run in blind execution (see parameter 'options')
if ~isempty(strfind(name,'b')) || ~isempty(strfind(name,'B'))
    tempSuffix = '.temp';
    nameOriginal = name;
    name = strrep(nameOriginal,'b','');
    name = [strrep(name,'B',''),tempSuffix];
    movefile([nameOriginal,'.cir'],[name,'.cir']);
end

% on Linux systems, a wine installation is required, see README.md
if ispc, [status,cmdout] = system(sprintf('%s %s.cir -run %s'           ,LTspiceEXE,name,options));
else,    [status,cmdout] = system(sprintf('%s %s %s.cir -run %s',wineEXE,LTspiceEXE,name,options));
end

% restore original filename if problem for blind execution detected (see above)
if exist('tempSuffix','var')
    files = dir([name,'.*']);
    for i = 1:length(files)
        file = files(i).name;
        fileOriginal = strrep(file,name,nameOriginal);
        movefile(file,fileOriginal);
    end
end

if status ~= 0, error('LTspice execution failed!'); end
fprintf('finished circuit simulation using LTspice after %d seconds.\n',toc(tstart));

end