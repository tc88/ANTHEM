function ipn = ipe2ipn(msh,ipe)
% IPE2IPN finds canonical node indices ipn of the nodes that are connected
% by the edge with canonical index ipe.
%
% Input:
%   msh   struct as defined by src/msh.txt
%         required fields: np,Mx,My,Mz
%   ipe   canonical indices of primary edges
%
% output parameters
%   ipn   canonical indices of primary nodes
%
% See also cv2dv, pn2pe
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% extract from msh struct
np = msh.np;

% reshaping to column vector
ipe = ipe(:);

% local indices of the primary edges in x-, y- and z-direction
ix = find(ipe<=np);
iy = find((ipe>np)&(ipe<=2*np));
iz = find(ipe>2*np);

% reference primary nodes
ipnx = ipe(ix);
ipny = ipe(iy)-np;
ipnz = ipe(iz)-2*np;

% indices of the associated primary nodes
ipn(ix,:) = [ ipnx -(ipnx+msh.Mx)   ];
ipn(iy,:) = [ ipny -(ipny+msh.My)  ];
ipn(iz,:) = [ ipnz -(ipnz+msh.Mz) ];

ipn = unique(abs(ipn(:)));

end