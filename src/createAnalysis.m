function analysisString = createAnalysis(analysis,varargin)
% CREATEANALYSIS generates valid SPICE analysis options as a netlist string
%
% Input:
%   filename   name of file to save netlist in
%   analysis   required analysis for the electric circuit. Must be one of
%              the following:
%              AC:   frequency (AC) analysis
%              tran: transient analysis
%   varargin   cell array containing extra variables for setting up the
%              analysis options. The different supported analyses expect
%              the following options:
%              AC:   fstart   starting frequency
%                    fstop    final frequency
%                    np       the number of frequency samples
%                    type     type of frequency sweep {lin,dec,oct}
%              tran: tend     end time
%
% Output:
%   analysisString  SPICE compatible string that describes the desired
%                   analysis to be run
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if strcmp(analysis,'AC')
    if (nargin-1) == 4
        fstart = varargin{1};
        fstop  = varargin{2};
        np     = varargin{3};
        type   = varargin{4};
        % returning AC analysis command as a string
        analysisString = sprintf('%s\t%s\t%d\t%-25.16e\t%-25.16e','.ac',type,np,fstart,fstop);
    else
        disp('too few arguments for AC analysis');
    end
elseif strcmp(analysis,'tran')
    tend = varargin{1};
    % returning transient analysis command as a string
    analysisString = sprintf('%s\t%-25.16e','.tran',tend);
else
    error('analysis type not implemented yet');
end

end