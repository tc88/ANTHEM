function [hnew,enew] = leapfrog(hold,eold,js,Mnu,Meps,Msigma,C,dt)
% LEAPFROG runs the Leapfrog algorithm to compute efit and hfit in time
% domain.
%
% Input:
%   hold    hfit at previous time step (3np-by-1)
%   eold    efit at previous time step (3np-by-1)
%   js      current excitation (3np-by-1)
%   Mnu     reluctance matrix (3np-by-3np)
%   Meps    capacitance matrix (3np-by-3np)
%   Msigma  conductance matrix (3np-by-3np)
%   C       curl matrix (3np-by-3np)
%   dt      current time step
%
% Output:
%   hnew    hfit at current time step
%   enew    efit at current time step
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

% compute new magnetic grid voltage
hnew = hold - dt*Mnu*C*eold;

% compute new electric grid voltage
enew = nullInv(Msigma+1/dt*Meps)*(1/dt*Meps*eold + C'*hnew - js);

end