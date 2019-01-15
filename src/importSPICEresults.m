function data = importSPICEresults(filename,simulator,verbose)
% IMPORTSPICERESULTS imports results from circuit simulator output.
%
% Input:
%   filename        name of file that contains the results
%   simulator       specifies which simulator is used. Must be one of
%                   {'LTspice','LTspiceXVII','LTspiceASCII','Xyce','Qucs'},
%                   where the 'ASCII' postfix requires LTspice results to
%                   be in ASCII format. By standard, LTspice output files
%                   are in binary format.
%   verbose         triggers console outputs and plots
%                   (optional, default: [1 0])
%
% Output:
%   data    output data structured as defined by LTspice2Matlab
%
% See also LTspice2Matlab, LTspiceASCII2Matlab, Xyce2Matlab
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if nargin < 3, verbose = [1 0]; end

switch simulator
    case 'LTspice'
        data = LTspice2Matlab(filename);
    case 'LTspiceXVII'
        data = LTspice2Matlab(filename,'all',1,'XVII');
    case 'LTspiceASCII'
        data = LTspiceASCII2Matlab(filename,verbose);
    case 'Xyce'
        data = Xyce2Matlab(filename,verbose);
    otherwise
        error('The provided simulator is not supported');
end

end
