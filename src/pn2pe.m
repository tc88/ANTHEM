function pe = pn2pe(msh,pn)
% PP2PE averages the scalar data pp located at primary nodes to the data
% pe located at primary edges.
%
% Input:
%   msh       struct as defined by src/msh.txt
%             required fields: np,Mx,My,Mz,ipeGhost
%   pn        scalar data at primary nodes (np-by-1)
%
% Output:
%   pe        scalar data at primary edges (3*np-by-1)
%
% See also cv2dv, ipe2ipn
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% get indices of non-ghost edges
ipeNonGhost = setdiff(1:3*msh.np,msh.ipeGhost);

% extract from msh struct
np = msh.np;

% dual points are always in the center of primary volumes!
pe = zeros(3*np,1);
for i = ipeNonGhost
    % average end point quantities of x-edges
    if i <= np
        pe(i) = (pn(i) + pn(i+msh.Mx))/2;
    % average end point quantities of y-edges
    elseif i > np && i <= 2*np
        pe(i) = (pn(i-np) + pn(i-np+msh.My))/2;
    % average end point quantities of z-edges
    elseif i > 2*np && i <= 3*np
        pe(i) = (pn(i-2*np) + pn(i-2*np+msh.Mz))/2;
    end
end

end