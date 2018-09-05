function signal = createSPICEsignal(tEnd,excitation)
% CREATESPICESIGNAL creates the netlist string that represents a signal
% over time.
%
% Currently supported signal types are: 'pulse','exp','sine' and 'constant'
%
% Input:
%   tEnd         end time (scalar)
%   excitation   struct as defined by src/excitation.txt
%                required fields: amplitude,delay,freq,phi,tau,type,t_rise
%
% Output:
%   signal      netlist string representing the desired signal
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

if strcmp(excitation.type,'pulse')
    signal = sprintf('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s','PULSE(0 ',excitation.amplitude,excitation.delay,excitation.t_rise,0,tEnd,tEnd,' 1)');
elseif strcmp(excitation.type,'exp')
    signal = sprintf('%s\t%d\t%d\t%d\t%d\t%d\t%s','EXP(0 ',excitation.amplitude,excitation.delay,excitation.tau,tEnd,excitation.tau,')');
elseif strcmp(excitation.type,'sine')
    offset = excitation.amplitude;
    signal = sprintf('%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s','SINE(',offset,excitation.amplitude,excitation.freq,0,0,excitation.phi,' 1)');
elseif strcmp(excitation.type,'constant')
    signal = sprintf('%d%s',excitation.amplitude,'Vdc');
else
    error('excitation type not supported');
end

end