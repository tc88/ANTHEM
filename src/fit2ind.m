function L = fit2ind(msh,Mnu)
% FIT2IND computes the lumped inductances associated with the primal edges
% of the mesh as required by the electric circuit representation of FIT.
% The inductance values are obtained from the magnetic reluctance
% matrix Mnu.
%
% Input:
%   msh     struct as defined by src/msh.txt
%           required fields: C
%   Mnu     reluctivity matrix (3np-by-3np)
%
% Output:
%   L       inductances associated with primal edges (3np-by-1)
% 
% See also cartMesh, createMeps, createMlambda, createMh, createMrhoc,
%          createMnu, and load_material
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% vector of reluctances
MnukkVec = spdiags(Mnu,0);

% dual curl matrix
Cd = msh.C';

% compute auxiliar mapping matrix
CdC = Cd.*msh.C.';

% compute vector MnuSum containing as entries MnummSum
MnuSum = CdC*MnukkVec;

% compute inductances on primary and dual edges
L = full(nullInv(MnuSum));

end