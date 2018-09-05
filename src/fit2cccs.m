function gI = fit2cccs(msh,Mnu,MmuSumVec)
% FIT2CCCS computes the current gain of the CCCSs associated with the
% primal edges of the mesh as required by the electric circuit
% representation of FIT. The current gains are obtained from the material
% permeability matrix and the curl matrices.
%
% Input:
%   msh        mesh structure as created by cartMesh
%              required fields: np,C
%   Mnu        reluctivity matrix of the medium (3np-by-3np)
%   MmuSumVec  inductances associated with primal edges (3np-by-1)
%
% Output:
%   gI         matrix of CCCS gains (3np-by-3np)
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

% dual curl matrix
Cd = msh.C';

% computing permeability matrix associated with primary edges
MmuSum = spdiags(MmuSumVec,0,3*msh.np,3*msh.np);

% gain for the CCCSs (diagonal is set to zero since no self-contribution)
gI = spdiags(sparse(3*msh.np,1),0,Cd*Mnu*msh.C*MmuSum);

end