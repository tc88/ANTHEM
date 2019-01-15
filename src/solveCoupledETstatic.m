function [phiSol,Tsol] = solveCoupledETstatic(msh,materials,idx,pots,T,verbose)
% SOLVECOUPLEDETSTATIC solves the stationary (nonlinear) electrothermal
% problem
%
% Input:
%   msh         mesh struct as created by cartMesh
%   materials   struct as defined by src/materials.txt
%               required fields: Msigma,Mlambda
%   idx         struct as defined by src/idx.txt
%               required fields: elect.dir,elect.dof,therm.dir,therm.dof
%   pots        electric potential. Entries for degrees of freedom
%               (dofs) need to be NaN while all other entries are
%               interpreted as fixed potentials (Dirichlet conditions)
%               (np-by-1)
%   T           struct for temperature data
%      .dir     vector of Dirichlet temperatures. DoFs need to be NaN while
%               all other entries are interpreted as fixed temperatures.
%               (optional, np-by-1)
%   verbose     triggers console output and plots
%               (optional, default: [1 0])
%
% Output:
%   phiSol     electric potential solution (np-by-1)
%   Tsol       temperature solution (np-by-1)
%
% See also computeQJ, solveCoupledET
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 9,  verbose = [1 0]; end

if ~isfield(idx.therm,'dir') || isempty(idx.therm.dir)
    if verbose(1), fprintf('no thermal Dirichlet nodes detected\n'); end
    T.dir = NaN*ones(msh.np,1);
end

if ~isempty(intersect(idx.therm.dof,idx.therm.dir))
    error('thermal Dirichlet nodes intersect with degrees of freedom!');
end

if ~isempty(intersect(idx.elect.dof,idx.elect.dir))
    error('electric Dirichlet nodes intersect with degrees of freedom!');
end

% extract variables from inputs
np = msh.np;
Sd = msh.Sd;
Mlambda = materials.Mlambda;
Msigma = materials.Msigma;

% initializations
phi.sol = zeros(np,1);
phi.sol(idx.elect.dir) = pots(idx.elect.dir);
phi.dir = pots(idx.elect.dir);
T.sol = zeros(np,1);
T.sol(idx.therm.dir) = T.dir(idx.therm.dir);

% thermal stiffness matrix
Kth = Sd*Mlambda*Sd';

% electrical mass and stiffness matrices
Kel = Sd*Msigma*Sd';
K11el = Kel(idx.elect.dof,idx.elect.dof);
K12el = Kel(idx.elect.dof,idx.elect.dir);
rhsEl = -K12el*phi.dir;

% solve electric system
if verbose(1), fprintf('solving electrical system ...\n'); end
phi.sol(idx.elect.dof) = K11el\rhsEl;

% compute thermal losses as a source for heat equation
Qel = computeQJ(msh,phi.sol,Msigma);

% create thermal system matrix
Ath = Kth;
A11th = Ath(idx.therm.dof,idx.therm.dof);
A12th = Ath(idx.therm.dof,idx.therm.dir);
if isempty(A12th) && isempty(idx.therm.dir)
    A12th = zeros(size(A12th,1),1);
    Tdir = 0;
else
    Tdir = T.dir(idx.therm.dir);
end

rhsTh = -A12th*Tdir + Qel(idx.therm.dof);

% solve thermal system
if verbose(1), fprintf('solving thermal system ...\n'); end
T.sol(idx.therm.dof) = A11th\rhsTh;

% return values
phiSol = phi.sol;
Tsol = T.sol;

end