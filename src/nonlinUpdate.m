function Msigma = nonlinUpdate(msh,materials,T,verbose)
% NONLINUPDATE updates nonlinear electric conductance matrix.
%
% Input:
%   msh         struct as defined by src/msh.txt
%               required fields: np,Mx,My,Mz,DAd,DS,ipeGhost
%   materials   struct as defined by src/materials.txt
%               required fields: nonlin.T0,nonlin.Dsigma0,nonlin.Dalpha
%   T           temperature on primary nodes (np-by-1)
%   verbose     triggers console outputs and plots
%               (optional, default: [1 0])
%
% Output:
%   Msigma      updated nonlinear electric conductance matrix
%
% See also isnonlin
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 2, verbose = [1 0]; end

if verbose(1), fprintf('updating nonlinear material matrix ...\n'); end

% extract from materials struct
T0 = materials.nonlin.T0;
Dsigma0 = materials.nonlin.Dsigma0;
Dalpha = materials.nonlin.Dalpha;

% reshape diagonal matrices to vectors
drho0 = 1./diag(Dsigma0);
dalpha = diag(Dalpha);

% compute temperatures on edges
TedgeAverage = pn2pe(msh,T);

% update resistivity vector according to temperature
drho = drho0.*(1+dalpha.*(TedgeAverage - T0));

% construct electric conductance matrices
Dsigma = spdiags(1./drho,0,3*msh.np,3*msh.np);
Msigma = msh.DAd * Dsigma * nullInv(msh.DS);

end