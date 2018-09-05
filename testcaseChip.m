function testcaseChip(LTspiceEXE,wineEXE,versionLTspice,verbose)
% TESTCASECHIP runs chip model testcase for circuit extraction
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running chip test case ...\n');

modelname = 'testcaseChip';                                                 % name of the model that is also used for output files

% settings
forceNewNetlistGen = false;                                                 % as a default, netlist generation is only done if .cir-file does not exist. To force new netlist generation, set this parameter to true
forceNewCirSim     = false;                                                 % as a default, circuit simulation is only done if .raw-file does not exist. To force new circuit simulation, set this parameter to true
refineFactorTime4Fit = 1;                                                   % ntFIT/ntCir, thus the factor the circuit time axis is refined to be used as FIT time axis

% load chip model
load(['./',modelname,'.mat']);

% extract netlist from linear electrothermal chip model
if ~exist(['./',modelname,'.cir'],'file') || forceNewNetlistGen
    extractNetlist4FITET(modelname,msh,idx,excitation,materials,time.t_end,T.init);
    % insert electric and thermal connection to represent the wire
    insertLumpedElement(modelname,wire.idx(find(wire.select),:),wire.Gel(find(wire.select)));
    insertLumpedElement(modelname,wire.idx(find(wire.select),:),wire.Gth(find(wire.select)),'lumpedT','T');
else
    warning('Netlist found, netlist generation skipped!');
end

% solve problem using circuit solver (results written in .raw-file)
if ~exist(['./',modelname,'.raw'],'file') || forceNewCirSim
    runLTspice(['./',modelname '.cir'],'-b',LTspiceEXE,wineEXE);
else
    warning('Circuit results found, circuit simulation skipped!');
end

% extract time axis and circuit results
data = LTspice2Matlab(['./',modelname '.raw'],'all',1,versionLTspice);
timeCir = data.time_vect;
phiCir  = spiceOrderOutput(data       );
tempCir = spiceOrderOutput(data,'','t');

% clean initial ramping from circuit result data
assert(all(T.init(1) == T.init))
[timeCir,tempCir,idxTime2keep] = spiceCleanInitRamp(timeCir,tempCir,T.init(1));
phiCir = phiCir(:,idxTime2keep);
timeFIT = refineAxis(timeCir,refineFactorTime4Fit);

% number of time steps
ntFIT = length(timeFIT);

% define excitation (Dirichlet conditions) for field problem
phiDir = NaN*ones(msh.np,ntFIT);
for i = 1:length(idx.elect.excitation)
    phiDirAnomFun = phiDirFun{idx.elect.excitation(i)};
    phiDir(idx.elect.excitation(i),:) = phiDirAnomFun(timeFIT);
    clear phiDirAnomFun
end
for i = 1:length(idx.elect.gnd)
    phiDir(idx.elect.gnd(i),:) = phiDirFun{idx.elect.gnd(i)};
end

% solve problem using field solver (results saved in variables phiFIT and tempFIT)
[phiFIT,tempFIT] = solveCoupledET(msh,materials,idx,phiDir,T,timeFIT,wire,verbose);

%% calculate errors
phiCirInterpol = spline(timeCir,phiCir,timeFIT);
tempCirInterpol = spline(timeCir,tempCir,timeFIT);

normedErrorPhi = zeros(1,ntFIT);
normedErrorTemp = zeros(1,ntFIT);
normedFITphi = zeros(1,ntFIT);
normedFITtemp = zeros(1,ntFIT);
for i = 1:ntFIT
        normedErrorPhi(i) = norm(phiCirInterpol(:,i)-phiFIT(:,i));
        normedFITphi(i) = norm(phiFIT(:,i));
        normedErrorTemp(i) = norm(tempCirInterpol(:,i)-tempFIT(:,i));
        normedFITtemp(i) = norm(tempFIT(:,i));
end
errPhi = max(normedErrorPhi)/max(normedFITphi);
fprintf('rel. error in phi: %f\n',errPhi);
errTemp = max(normedErrorTemp)/max(normedFITtemp);
fprintf('rel. error in temperature: %f\n',errTemp);

%% write .vtr file for Paraview
% write to VTK file
fit_write_vtk(msh.x,msh.y,msh.z,['./',modelname,'Results.vtr'], ...
    {'PotentialFIT',phiFIT(:,ntFIT); ...
     'TemperatureFIT',tempFIT(:,ntFIT); ...
     'PotentialCir',phiCirInterpol(:,ntFIT); ...
     'TemperatureCir',tempCirInterpol(:,ntFIT)});

%% plot hottest point in insulator
[~,idx_hottestFIT] = max(sum(tempFIT,2));
idx_hottest = idx_hottestFIT; % force to use FIT's hottest

figure(19); clf; hold all;
plot(timeFIT,tempFIT(idx_hottest,:));
plot(timeCir,tempCir(idx_hottest,:),'--');
xlabel('Time in s');
ylabel('Temperature in K');
legend('FIT','Circuit','Location','NorthWest');
print(['./',modelname,'.pdf'],'-dpdf');

fprintf('finished chip test case after %d seconds.\n',toc(tstart));