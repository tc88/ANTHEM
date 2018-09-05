function QJ = computeQJ(msh,phi,Msigma)
% COMPUTEQJ computes the thermal losses from the electrical
% solution phi and the conductance matrix Msigma. This is done by
% computing the thermal losses first on the control volumes Vhat and then
% by relocating these losses to the dual volumes Vdual.
%
% Input:
%   msh       struct as defined by src/msh.txt
%             required fields: np,Mx,My,Mz,Sd
%   phi       electric potential (np-by-1)
%   Msigma    conductance matrix (3np-by-3np)
%
% Output:
%   QJ        thermal powers located on dual volumes (np-by-1)
%
% See also createQJ4edge, createQJ4point
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% compute integrated electrical field
ebow = msh.Sd'*phi;

% compute thermal losses
Qhat = diag(Msigma).*(ebow.^2);
QJ = cv2dv(msh,Qhat);

end