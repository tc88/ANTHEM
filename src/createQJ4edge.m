function QJiSubString = createQJ4edge(m,msh,materials)
% CREATEQJ4EDGE returns the string representing the value of QJ for one
% edge. Inhomogeneous domains and nonlinear conductivities are supported.
%
% Input:
%   m               index of considered edge
%   msh             struct as defined by src/msh.txt
%                   required fields: np,Mx,My,Mz,DS
%   materials       struct as defined by src/materials.txt
%                   required fields: Msigma
%
% Output:
%   QJiSubString    string representing QJ for edge ipe
%
% See also createQJ4point, createResistor, createSPICEsignal, computeQJ
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% find indices of primary nodes that are connected by edge of index ipe
temp = ipe2ipn(msh,m);
ipn1 = temp(1); ipn2 = temp(2);
clear temp

% extract nonlinear information
if isnonlin(materials)
    Tref = materials.nonlin.T0;
    Mrho0Diag = 1./diag(materials.nonlin.Msigma0);
    DalphaDiag = full(diag(materials.nonlin.Dalpha));
end

% create QJiSubString for inhomogeneous domain and nonlinear materials
if isnonlin(materials)
    % check whether this edge is not conducting
    if (msh.DS(m,m) == 0) || isinf(Mrho0Diag(m))
        sigmam = '0';
    else
        sigmam = sprintf('%s%d%s%d%s%d%s%d%s%d%s','1/(',Mrho0Diag(m),'*(1+',DalphaDiag(m), ...
                         '*(((V(',ipn1,'T)+V(',ipn2,'T))*0.5)-',Tref,')))');
    end
    QJiSubString = sprintf('%s%d%s%d%s%s','0.5*pow(V(',ipn1,',',ipn2,'),2)*',sigmam);
else
    sigmam = num2str(full(materials.Msigma(m,m)));
    QJiSubString = sprintf('%s%d%s%d%s%s','0.5*pow(V(',ipn1,',',ipn2,'),2)*',sigmam);
end

end