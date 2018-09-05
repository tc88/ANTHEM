function bool = isnonlin(materials)
% ISNONLIN checks whether nonlinear electric conductivities are involved.
%
% Input:
%   materials   struct as defined by src/materials.txt
%
% Output:
%   bool    true if nonlinear electric conductivity is involved
%
% See also nonlinUpdate
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% check is based on existence of non-zero entries in alpha
bool = isfield(materials,'nonlin') && isfield(materials.nonlin,'alpha') && any(materials.nonlin.alpha);

end