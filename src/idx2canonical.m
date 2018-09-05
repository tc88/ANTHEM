function canonicalIdx = idx2canonical(msh,i,j,k)
% IDX2CANONICAL translates given indices (i,j,k) into the canonical
% indexing scheme. One always has to define a single point, a line or a
% cube to be translated.
%
% Input:
%   msh                 mesh struct as created by cartMesh
%                       required fields: Mx,My,Mz
%   i                   can be a scalar or a vector for the x-index
%   j                   can be a scalar or a vector for the y-index
%   k                   can be a scalar or a vector for the z-index
%
% Output:
%   idx     a vector of all indices specified by the given indices
%           i,j and k in the canonical indexing scheme
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if ~isrow(i), i = i'; end
if ~isrow(j), j = j'; end
if ~isrow(k), k = k'; end

% calculate canonical index
canonicalIdx = 1+(i-1)*msh.Mx+(j-1)*msh.My+(k-1)*msh.Mz;

end