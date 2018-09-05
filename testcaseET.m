function testcaseET(modelsize,LTspiceEXE,wineEXE,versionLTspice,verbose)
% TESTCASEET implements the electrothermal circuit validation considering
% linear and nonlinear conductivities.
%
% Corresponds to the results presented in Section 7.1 of the paper.
% The model used here is shown in Fig. 17.
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running electrothermal test case ...\n');

% name of the model that is also used for output files
modelname = 'testcaseET';

% capitalize modelsize
modelsize = strcat(upper(modelsize(1)),modelsize(2:end));

% refineFactorTime4Fit = ntFIT/ntCir, where nt{FIT,Cir} is the number of time steps used for the {FIT,Cir} simulation
refineFactorTime4Fit = 3;

%% linear model

tstartLin = tic;
fprintf('running linear electrothermal problem ...\n');

% filename for linear case
filenameLin = ['./',modelname,'lin',modelsize];

% load linear model
load([filenameLin,'.mat']);

% extract netlist from linear electrothermal model
extractNetlist4FITET(filenameLin,msh,idx,excitation,materials,time.t_end,T.ref);

% solve linear problem using circuit solver (results written in .raw-file)
runLTspice(filenameLin,'-b',LTspiceEXE,wineEXE);

% extract time axis and circuit results
spiceDataLin = LTspice2Matlab([filenameLin,'.raw'],'all',1,versionLTspice);
timeCirLin = spiceDataLin.time_vect;
timeFITlin = refineAxis(timeCirLin,refineFactorTime4Fit);
ntFITlin   = length(timeFITlin);
phiCirLin  = spiceOrderOutput(spiceDataLin);
tempCirLin = spiceOrderOutput(spiceDataLin,'','t');

% define excitation (Dirichlet conditions) for field problem
vpotential = NaN*ones(msh.np,ntFITlin);
vpotential(idx.elect.excitation,:) = repmat(excitation.signal(timeFITlin),numel(idx.elect.excitation),1);
vpotential(idx.elect.gnd,:) = 0;

% solve linear problem using field solver (results saved in variables phiFITlin and tempFITlin)
[phiFITlin,tempFITlin] = solveCoupledET(msh,materials,idx,vpotential,T,timeFITlin,[],verbose);

fprintf('finished linear electrothermal problem after %d seconds.\n',toc(tstartLin));
clear msh materials idx excitation time

%% nonlinear model

tstartNonlin = tic;
fprintf('running nonlinear electrothermal problem ...\n');

% filename of generated nonlinear netlist
filenameNonlin = [modelname,'nonlin',modelsize];

% load nonlinear model
load([filenameNonlin,'.mat']);

% extract netlist from model
extractNetlist4FITET(filenameNonlin,msh,idx,excitation,materials,time.t_end,T.ref);

% solve nonlinear problem using circuit solver (results written in .raw-file)
runLTspice(filenameNonlin,'-b',LTspiceEXE,wineEXE);

% extract time axis and circuit results
spiceDataNonlin = LTspice2Matlab([filenameNonlin,'.raw'],'all',1,versionLTspice);
timeCirNonlin = spiceDataNonlin.time_vect;
timeFITnonlin = refineAxis(timeCirNonlin,refineFactorTime4Fit);
ntFITnonlin = length(timeFITnonlin);
phiCirNonlin = spiceOrderOutput(spiceDataNonlin);
tempCirNonlin = spiceOrderOutput(spiceDataNonlin,'','t');

% define excitation (Dirichlet conditions) for field problem
vpotential = NaN*ones(msh.np,ntFITnonlin);
vpotential(idx.elect.excitation,:) = repmat(excitation.signal(timeFITnonlin),numel(idx.elect.excitation),1);
vpotential(idx.elect.gnd,:) = 0;

% solve nonlinear problem using field solver (results saved in variables phiFITnonlin and tempFITnonlin)
[phiFITnonlin,tempFITnonlin] = solveCoupledET(msh,materials,idx,vpotential,T,timeFITnonlin,[],verbose);

fprintf('finished linear electrothermal problem after %d seconds.\n',toc(tstartNonlin));

%% calculate linear errors

% clean ramping to initial temperature from circuit results (as not present in field results)
[timeCirLin,tempCirLin,idxTime2keepLin] = spiceCleanInitRamp(timeCirLin,tempCirLin,T.ref);
phiCirLin = phiCirLin(:,idxTime2keepLin);

% interpolate circuit results to finer time axis of the field solver
phiCirLinInterpol = spline(timeCirLin,phiCirLin,timeFITlin);
tempCirLinInterpol = spline(timeCirLin,tempCirLin,timeFITlin);

% calculate error norms as desribed in paper
normTempFITlin = zeros(1,ntFITlin);
normTempCirLinInterpol = zeros(1,ntFITlin);
normPhiFITlin = zeros(1,ntFITlin);
normPhiDiffLin = zeros(1,ntFITlin);
normTempDiffLin = zeros(1,ntFITlin);
for i = 1:ntFITlin
    normPhiFITlin(i) = norm(phiFITlin(:,i));
    normPhiDiffLin(i) = norm(phiCirLinInterpol(:,i)-phiFITlin(:,i));
    normTempFITlin(i) = norm(tempFITlin(:,i));
    normTempDiffLin(i) = norm(tempCirLinInterpol(:,i)-tempFITlin(:,i));
    normTempCirLinInterpol(i) = norm(tempCirLinInterpol);
end
deltaPhiLin = max(normPhiDiffLin)/max(normPhiFITlin);
deltaTempLin = max(normTempDiffLin)/max(normTempFITlin);
display(['linear el. potential error: ' num2str(deltaPhiLin)]);
display(['linear temperature error: ' num2str(deltaTempLin)]);

%% calculate nonlinear errors

% clean ramping to initial temperature from circuit results (as not present in field results)
[timeCirNonlin,tempCirNonlin,idxTime2keepNonlin] = spiceCleanInitRamp(timeCirNonlin,tempCirNonlin,T.ref);
phiCirNonlin = phiCirNonlin(:,idxTime2keepNonlin);

% interpolate circuit results to finer time axis of the field solver
phiCirNonlinInterpol = spline(timeCirNonlin,phiCirNonlin,timeFITnonlin);
tempCirNonlinInterpol = spline(timeCirNonlin,tempCirNonlin,timeFITnonlin);

% calculate error norms as desribed in paper
normPhiFITnonlin = zeros(1,ntFITnonlin);
normPhiDiffNonlin = zeros(1,ntFITnonlin);
normTempFITnonlin = zeros(1,ntFITnonlin);
normTempDiffNonlin = zeros(1,ntFITnonlin);
for i = 1:ntFITnonlin
    normPhiFITnonlin(i) = norm(phiFITnonlin(:,i));
    normPhiDiffNonlin(i) = norm(phiCirNonlinInterpol(:,i)-phiFITnonlin(:,i));
    normTempFITnonlin(i) = norm(tempFITnonlin(:,i));
    normTempDiffNonlin(i) = norm(tempCirNonlinInterpol(:,i)-tempFITnonlin(:,i));
end
deltaPhiNonlin = max(normPhiDiffNonlin)/max(normPhiFITnonlin);
deltaTempNonlin = max(normTempDiffNonlin)/max(normTempFITnonlin);
display(['nonlinear el. potential error: ' num2str(deltaPhiNonlin)]);
display(['nonlinear temperature error: ' num2str(deltaTempNonlin)]);

%% plot comparison of field and circuit results
idx2cmp = ceil(Nxyz(1)/2);
figure(18); clf; subplot(1,2,1);
plot(1e6*timeFITlin,phiFITlin(idx2cmp,:),'*'); hold on;
plot(1e6*timeFITnonlin,phiFITnonlin(idx2cmp,:),'x');
plot(1e6*timeCirLin,phiCirLin(idx2cmp,:));
plot(1e6*timeCirNonlin,phiCirNonlin(idx2cmp,:));
xlabel('time in {\mu}s');
ylabel('Potential $$\varphi|_{\mathbf{x}_{0}}$$ in V','Interpreter','Latex');
legend('FIT linear','FIT non-linear','Circuit linear','Circuit non-linear','Location','southeast');

subplot(1,2,2);
plot(1e6*timeFITlin,tempFITlin(idx2cmp,:),'*'); hold on;
plot(1e6*timeFITnonlin,tempFITnonlin(idx2cmp,:),'x');
plot(1e6*timeCirLin,tempCirLin(idx2cmp,:));
plot(1e6*timeCirNonlin,tempCirNonlin(idx2cmp,:));
xlabel('time in {\mu}s');
ylabel('Temperature $$\mathbf{T}|_{\mathbf{x}_{0}}$$ in K','Interpreter','Latex');
legend('FIT linear','FIT non-linear','Circuit linear','Circuit non-linear','Location','southeast');

print([modelname,lower(modelsize),'.pdf'],'-dpdf');

fprintf('finished electrothermal test case after %d seconds.\n',toc(tstart));