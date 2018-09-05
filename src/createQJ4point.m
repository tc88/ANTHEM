function QJiString = createQJ4point(i,msh,materials)
% CREATEQJ4POINT returns the string representing the value of QJ for one
% point.
%
% Input:
%   i           index of considered point
%   msh         struct as defined by src/msh.txt
%               required fields: np,Mx,My,Mz,ipnXmin,ipnXmax,ipnYmin,
%                                ipnYmax,ipnZmin,ipnZmax
%   materials   struct as defined by src/materials.txt
%               required fields: Msigma
%
% Output:
%   QJiString   string representing QJ for point i
%
% See also createQJ4edge, createResistor, createSPICEsignal, computeQJ
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% create string for all 6 edges touching point Pi
if ismember(i,msh.ipnXmin), QJxMinus = '0';
else, QJxMinus = createQJ4edge(i-msh.Mx         ,msh,materials);
end
if ismember(i,msh.ipnXmax), QJxPlus  = '0';
else, QJxPlus  = createQJ4edge(i                ,msh,materials);
end
if ismember(i,msh.ipnYmin), QJyMinus = '0';
else, QJyMinus = createQJ4edge(msh.np+i-msh.My  ,msh,materials);
end
if ismember(i,msh.ipnYmax), QJyPlus  = '0';
else, QJyPlus  = createQJ4edge(msh.np+i         ,msh,materials);
end
if ismember(i,msh.ipnZmin), QJzMinus = '0';
else, QJzMinus = createQJ4edge(2*msh.np+i-msh.Mz,msh,materials);
end
if ismember(i,msh.ipnZmax), QJzPlus  = '0';
else, QJzPlus  = createQJ4edge(2*msh.np+i       ,msh,materials);
end
QJiString = sprintf('%s+%s+%s+%s+%s+%s',QJxMinus,QJxPlus,QJyMinus,QJyPlus,QJzMinus,QJzPlus);

end