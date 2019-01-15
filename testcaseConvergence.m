function testcaseConvergence(LTspiceEXE,wineEXE,verbose)
% TESTCASECONVERGENCE runs convergence model testcase for circuit extraction
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

tstart = tic;
fprintf('running convergence test case ...\n');
modelname = 'testcaseConvergence';                                         % name of the model that is also used for output files

% settings
Neval = 1e5;                                                               % number of points to calculate error

% load convergence model
data = load(['./',modelname,'.mat']);

% extract from loaded data
msh4refine = data.msh4refine;                                              % []       : mesh as defined by src/msh.txt for different levels of refinement
idx4refine = data.idx4refine;                                              % []       : indices as defined by src/idx.txt for different levels of refinement
materials4refine = data.materials4refine;                                  % []       : materials as defined by src/materials.txt for different levels of refinement
phiDir4refine = data.phiDir4refine;                                        % [V]      : Dirichlet voltages for different levels of refinement
T4refine = data.T4refine;                                                  % [K]      : Dirichlet temperatures for different levels of refinement
Nvec = data.Nvec(1:3); % Nvec = data.Nvec;                                 % [#]      : number of 1D points to consider
excitation = data.excitation;                                              % []       : excitation as defined by src/excitation.txt
c = data.c;                                                                % [S/m]    : electric conductivity for 0 < x < ell
C = data.C;                                                                % [S/m]    : electric conductivity for ell < x < L
d = data.d;                                                                % [W/(Km)] : thermal conductivity for 0 < x < ell
D = data.D;                                                                % [W/(Km)] : thermal conductivity for ell < x < L
ell = data.ell; L = data.L;                                                % [m]      : lengths of regions
V = data.V;                                                                % [V]      : electric Dirichlect conditions
Tl = data.Tl; Tr = data.Tr;                                                % [K]      : thermal Dirichlet conditions

if Neval/Nvec(end) < 10
    error('number of evaluation points should be at least one magnitude greater than the finest resolution');
end

DoFfit    = zeros(1,length(Nvec));
DoFcir    = zeros(1,length(Nvec));
errFIT    = zeros(1,length(Nvec));
errLT     = zeros(1,length(Nvec));
errXyce   = zeros(1,length(Nvec));
for i = 1:length(Nvec)
    N = Nvec(i);
    msh = msh4refine{i};
    idx = idx4refine{i};
    materials = materials4refine{i};
    phiDir = phiDir4refine{i};
    T = T4refine{i};

    % generate netlist
    if ~all(unique(phiDir(idx.elect.dir))==unique(T.dir(idx.therm.dir)))
        error('netlist extraction only supports electric and thermal Dirichlet conditions to be the same');
    end
    extractNetlist4FITET(sprintf('%sN%d',modelname,N),msh,idx,excitation,materials,0,293);

    % run LTspice and import results
    runLTspice(sprintf('%sN%d',modelname,N),'-ascii -b',LTspiceEXE,wineEXE);
    if exist(sprintf('./%sN%d.raw',modelname,N),'file')
        dataLT = importSPICEresults(sprintf('%sN%d.raw',modelname,N),'LTspiceASCII',verbose);
        phiLT = spiceOrderOutput(dataLT,'','',verbose);
        TLT = spiceOrderOutput(dataLT,'','t',verbose);
        % take last value since one is forced to use a DC ramping in SPICE
        phiLT = phiLT(:,end);
        TLT = TLT(:,end);
    else
        phiLT = NaN*ones(msh.np,1);
        TLT = NaN*ones(msh.np,1);
    end

    % import Xyce results
    if exist(sprintf('./%sN%dXyce.prn',modelname,N),'file')
        dataXyce = importSPICEresults(sprintf('%sN%dXyce.prn',modelname,N),'Xyce',verbose);
        Txyce = spiceOrderOutput(dataXyce,'','T,0',verbose);
        % take last value since one is forced to use a DC ramping in SPICE
        Txyce = Txyce(:,end);
    else
        Txyce = NaN*ones(msh.np,1);
    end

    % run field solver
    [phiFIT,TFIT] = solveCoupledETstatic(msh,materials,idx,phiDir,T,verbose);

    % calculate error
    xEval = linspace(0,max(msh.x),Neval)';
    [phiAna,Tana] = anaSolET(xEval,c,C,d,D,ell,L,V,Tl,Tr);
    TFITeval  = interp1(msh.x, TFIT(1:msh.nx),xEval);
    TLTeval   = interp1(msh.x,  TLT(1:msh.nx),xEval);
    TxyceEval = interp1(msh.x,Txyce(1:msh.nx),xEval);
    DSd = spdiags([xEval(1)/2;(xEval(1:end-1)+xEval(2:end))/2;xEval(end)/2],0,Neval,Neval);
    errFIT(i)  = (TFITeval-Tana)'*DSd*(TFITeval-Tana)/(Tana'*DSd*Tana);
    errLT(i)   = (TLTeval-Tana)'*DSd*(TLTeval-Tana)/(Tana'*DSd*Tana);
    errXyce(i) = (TxyceEval-Tana)'*DSd*(TxyceEval-Tana)/(Tana'*DSd*Tana);

    % save DoFs
    DoFfit(i) = 2*msh.np;
    DoFcir(i) = dataLT.logdata.matrixSize;
end

% plots
figure(1); clf;
subplot(2,2,1); hold on;
plot(xEval,phiAna);
plot(msh.x,phiFIT(1:msh.nx),'-x');
plot(msh.x,phiLT(1:msh.nx),'-x');
xlabel('x in m');
ylabel('\Phi(x) in V');
legend('exact','FIT','LTspice','Location','SouthEast');

subplot(2,2,2); hold on;
plot(xEval,Tana);
plot(msh.x,TFIT(1:msh.nx),'-x');
plot(msh.x,TLT(1:msh.nx),'-x');
xlabel('x in m');
ylabel('T(x) in K');
legend('exact','FIT','LTspice','Location','SouthEast');

subplot(2,2,3);
loglog(DoFfit,errFIT,'-x'); hold on;
xlabel('DoF (FIT)');
ylabel('rel. error of field simulation');

subplot(2,2,4);
loglog(DoFcir,errLT,'-x'); hold on;
loglog(DoFcir,errXyce,'-x'); hold on;
legend('LTspice','Xyce');
xlabel('DoF (Cir)');
ylabel('rel. error of circuit simulation');

fprintf('finished convergence test case after %d seconds.\n',toc(tstart));

end

function [phi,T] = anaSolET(x,c,C,d,D,ell,L,V,Tl,Tr)
% ANASOLET returns the analytical solution of the 1D static ET problem
    N = length(x);
    idx1 = find(x<ell);
    idx2 = idx1(end)+1:N;
    phi1 = C*V*x    /((C-c)*ell+c*L);
    phi2 = c*V*(x-L)/((C-c)*ell+c*L) + V;
    q = c*C^2*V^2/((C-c)*ell+c*L)^2;
    Q = c^2*C*V^2/((C-c)*ell+c*L)^2;
    T1 = -q/(2*d)*x.^2+(-2*q*ell^2*d+2*q*ell*L*d+Q*ell^2*d-2*Q*ell*L*d+Q*L^2*d+q*ell^2*D-2*Tl*d*D+2*Tr*d*D)/(2*ell*d*D-2*ell*d^2+2*L*d^2)*x+Tl;
    T2 = -Q/(2*D)*(x.^2-L*x)+(-q*ell^2*D+2*Q*ell^2*D-Q*ell^2*d+Q*ell*L*d-2*Tl*d*D+2*Tr*d*D-Q*ell*L*D)/(2*ell*D^2-2*ell*d*D+2*L*d*D)*(x-L)+Tr;
    phi = zeros(N,1);
    phi(idx1) = phi1(idx1);
    phi(idx2) = phi2(idx2);
    T = zeros(N,1);
    T(idx1) = T1(idx1);
    T(idx2) = T2(idx2);
end
