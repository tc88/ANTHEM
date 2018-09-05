function GmString = createResistor(m,i,j,materials,matType)
% CREATERESISTOR returns the SPICE netlist entry for an electric resistor.
%
% Input:
%   m           canonical index of considered primary edge
%   i           canonical index of first primary node that touches edge m
%   j           canonical index of second primary node that touches edge m
%   materials   struct as defined by src/materials.txt
%               required fields: Msigma or Mlambda (depending on matType)
%   matType     type of material for which resistor shall be created.
%               Must be either 'sigma' or 'lambda'
%               (optional, default: 'sigma')
%
% Output:
%   GmString    Gm as SPICE netlist compatible string
%
% See also createQJ4edge, createQJ4point, createSPICEsignal
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 5, matType = 'sigma'; end

% check whether material is nonlinear
if isnonlin(materials)
    % extract nonlinear material information
    Tref = materials.nonlin.T0;
    switch matType
        case 'sigma'
            Mres0Diag = 1./diag(materials.nonlin.Msigma0);
            DtempCoeffDiag = full(diag(materials.nonlin.Dalpha));
        case 'lambda'
            Mres0Diag = 1./diag(materials.nonlin.Mlambda0);
            DtempCoeffDiag = full(diag(materials.nonlin.Dbeta));
    end

    % Gm for a nonlinear conductivity
    if isinf(Mres0Diag(m))
        GmString = '0';
    else
        GmString = sprintf('%s%d%s%d%s%d%s%d%s%d%s','1/(',Mres0Diag(m),'*(1+',DtempCoeffDiag(m),...
                         '*(((V(',i,'T)+V(',j,'T))*0.5)-',Tref,')))');
    end
else
    % Gm for a linear conductivity
    switch matType
        case 'sigma', GmString = sprintf('%d',full(materials.Msigma(m,m)));
        case 'lambda',GmString = sprintf('%d',full(materials.Mlambda(m,m)));
    end
end

end