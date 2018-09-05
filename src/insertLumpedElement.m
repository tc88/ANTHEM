function insertLumpedElement(filename,ipn,value,suffix,idxSuffix)
% INSERTLUMPEDELEMENT inserts a lumped element at the beginning of a
% given netlist.
%
% Input:
%   filename    name of the netlist
%   ipn         indices of end points of lumped element (1-by-2)
%   value       value of the lumped element to be inserted
%   suffix      string appended to the identifier of the lumped element
%               (optional, default: 'lumped')
%   idxSuffix   string appended to the node number (optional, default: '')
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 5, idxSuffix = '';       end
if nargin < 4, suffix    = 'lumped'; end

% line number where the entry should be added to the netlist
insertLine = 2;

% read current content of netlist
fileID = fopen([filename '.cir'],'r');
netlistText = textscan( fileID, '%s', 'Delimiter','\n');
fclose(fileID);
netlistText = netlistText{1};

% write old part before the new line
fid = fopen([filename,'.cir'],'w');
for i = 1:(insertLine-1)
    fprintf(fid,'%s\n',char(netlistText(i)));
end

% write new line
fprintf(fid,'%s%d%s%d%s\t%d%s\t%d%s\t%d\n','R',ipn(1),'to',ipn(2),suffix,ipn(1),idxSuffix,ipn(2),idxSuffix,value);

% write old part after the new line
for i = insertLine:length(netlistText)
    fprintf(fid,'%s\n',char(netlistText(i)));
end

% close netlist file
fclose(fid);