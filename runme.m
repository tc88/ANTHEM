% this script calls the different test cases to reproduce all numerical
% results from the following paper:
%
% T. Casper, D. Duque, S. Schöps, H. De Gersem. Automated Netlist 
% Generation for 3D Electrothermal and Electromagnetic Field Problems.
% Journal of Computational Electronics, submitted. arXiv:1809.08588.
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

clear;
close all;
addpath(genpath('src'));
verbose = [0 0];                                                           % triggers console outputs and plots
modelsize = 'large';                                                       % {'small','large'}
wineEXE = '/usr/bin/wine';                                                 % path to wine (only required for UNIX systems)
LTspiceEXE = '~/.wine/drive_c/Program Files (x86)/LTC/LTspiceIV/scad3.exe';% path to LTspice (may include whitespaces and parantheses)

% if Octave is used, add compatOctave to path for compatibility
if isOctave, addpath('compatOctave'); end

% find LTspice version from LTspice path
if ~isempty(strfind(LTspiceEXE,'LTspiceIV'))
    versionLTspice = 'IV';
elseif ~isempty(strfind(LTspiceEXE,'LTspiceXVII'))
    versionLTspice = 'XVII';
else
    error('expected ''LTspiceEXE'' to contain either ''LTspiceIV'' or ''LTspiceXVII''');
end

% run individual test cases
testcaseET(modelsize,LTspiceEXE,wineEXE,versionLTspice,verbose);
testcaseConvergence(LTspiceEXE,wineEXE,verbose);                           % only in official journal version
testcaseChip(LTspiceEXE,wineEXE,versionLTspice,verbose);
testcaseEM(modelsize,LTspiceEXE,wineEXE,versionLTspice,verbose);
testcaseABC(modelsize,LTspiceEXE,wineEXE,versionLTspice,verbose);