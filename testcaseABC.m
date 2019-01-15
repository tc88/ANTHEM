function testcaseABC(modelsize,LTspiceEXE,wineEXE,versionLTspice)
% TESTCASEABC is a test case that implements absorbing boundary
% conditions (ABC) in an exemplary setting. A coaxial cable of rectangular
% cross section is used to simulate wave propagation along the cable from
% port 1 to port 2. To avoid reflections at port 2, the cable is terminated
% with its wave impedance.
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

modelname = 'testcaseABC';                                                 % name of the model that is also used for output files
terminationCell = {'NoReflect','Reflect'};

for i = 1:length(terminationCell)
    
    termination = terminationCell{i};
    tstart = tic;
    fprintf('running absorbing boundary condition test case with %s termination ...\n',termination);

    % settings
    forceNewCirSim = false;                               % as a default, circuit simulation is only done if .raw-file does not exist. To force new circuit simulation, set this parameter to true
    refineFactorTime4Fit = 1;                             % ntFIT/ntCir, thus the factor the circuit time axis is refined to be used as FIT time axis
    t2plot = 1e-9*[3 5 7 8];                              % [s]  : time points for which U_oi shall be plotted
    tend = 1e-8;                                          % [s]  : end time
    dt4excite = 2e-11;                                    % [s]  : time spacing used for excitation signal
    fmax_gauss = 1e9;                                     % [Hz] : maximum frequency in Gauss signal
    ipe2measureP1 = 5;                                    % [#]  : primary edge index at port 1 used for current and voltage measurements
    fileNetlist = ['./',modelname,modelsize,termination]; % [#]  : filename of generated netlist

    % load problem
    load(['./',modelname,modelsize,'.mat']);
    np = msh.np;

    % calculate characteristic impedance of line in Ohm
    Z0 = (max(max(Mnu))*max(max(Meps)))^(-0.5)/8;

    % value of termination resistance in Ohm
    if strcmp(termination,'Reflect'), R = Inf;
    elseif strcmp(termination,'NoReflect'), R = Z0;
    else, error('non-supported termination chosen');
    end

    % terminate port 2 of the transmission line with R
    ipe4RatP1 = [[5,7,9,11] [2,3,10,11]+np];
    ipe4RatP2 = ipe4RatP1 + (msh.nz-1)*msh.Mz;
    Msigma = nullInv(sparse(ipe4RatP2,ipe4RatP2,8*R*ones(numel(ipe4RatP2),1),3*np,3*np));

    % primary edge index at port 2 used for current and voltage measurements
    ipe2measureP2 = ipe2measureP1+(msh.nz-1)*msh.Mz;

    % define indices of dual facets used for excitation
    idfExciteMinus = [5,9,np+2,np+3];
    idfExcitePlus = [7,11,np+10,np+11];

    % find all indices of primary edges for voltage measurement
    ipe2measure = ipe2measureP1:msh.Mz:ipe2measureP2;

    % define excitation signal
    time4excite = 0:dt4excite:tend;
    omega = 2*pi*fmax_gauss;
    sigmaG = sqrt(2*log(100)/omega^2);
    tm = time4excite - sqrt(2*sigmaG^2 * log(1000));
    Ii = exp(-tm.^2 ./ (2*sigmaG^2));

    % export current source to ASCII file
    filename = ['./',modelname,'currentSourcePos.csv'];
    dlmwrite(filename,[time4excite',Ii'/8]);
    filename = ['./',modelname,'currentSourceNeg.csv'];
    dlmwrite(filename,[time4excite',-Ii'/8]);

    % use ASCII file as PWL excitation for circuit simulation
    Isetting.type = 'pwl';
    Isetting.amp = cell(3*msh.np,1);
    for idx = idfExcitePlus
        Isetting.amp{idx} = ['./',modelname,'currentSourcePos.csv'];
    end
    for idx = idfExciteMinus
        Isetting.amp{idx} = ['./',modelname,'currentSourceNeg.csv'];
    end

    % extract netlist from electromagnetic field model
    analysis = 'tran';
    extractNetlist4FITEM(fileNetlist,msh,Meps,Msigma,Mnu,Isetting,analysis,ipePEC,tend);

    % solve electromagnetic problem using circuit solver (results written in .raw file)
    if ~exist([fileNetlist,'.raw'],'file') || forceNewCirSim
        runLTspice(fileNetlist,'-b',LTspiceEXE,wineEXE);
    else
        warning('Circuit results found, circuit simulation skipped!')
    end

    % extract time axis and circuit results
    spiceData = LTspice2Matlab([fileNetlist '.raw'],'all',1,versionLTspice);
    timeCir = spiceData.time_vect;
    timeFIT = refineAxis(timeCir,refineFactorTime4Fit);
    ntFIT = length(timeFIT);
    espice = spiceOrderOutput(spiceData);
    V2Cir  = -espice(ipe2measureP2,:);
    VoiCir = -espice(ipe2measure,:);

    % find time indices to plot
    if isOctave
        [~,it2plotFIT] = ismembertol(t2plot,timeFIT,'tol',1e-3);
        [~,it2plotCir] = ismembertol(t2plot,timeCir,'tol',1e-3);
    else 
        [~,it2plotFIT] = ismembertol(t2plot,timeFIT,1e-3);
        [~,it2plotCir] = ismembertol(t2plot,timeCir,1e-3);
    end
    fprintf('time to plot: %e s\n',timeFIT(it2plotFIT));

    % define excitation signal for time axis used by field solver
    Ii = interp1(time4excite,Ii,timeFIT);

    % set current excitation such that it is distributed over the entire face of port 1
    jifit = zeros(3*np,ntFIT);
    jifit(idfExciteMinus,:) = repmat(-Ii/8,numel(idfExciteMinus),1);
    jifit(idfExcitePlus,:)  = repmat( Ii/8,numel(idfExcitePlus) ,1);

    % initializations
    hfit=zeros(3*np,ntFIT);
    efit=zeros(3*np,ntFIT);
    V1FIT = zeros(1,ntFIT);
    I1FIT = zeros(1,ntFIT);
    V2FIT = zeros(1,ntFIT);
    I2FIT = zeros(1,ntFIT);

    % solve electromagnetic problem using field solver in time domain (leapfrog)
    tstartFIT = tic;
    fprintf('running field simulation using FIT ...\n');
    VoiFIT = zeros(msh.nz,ntFIT);
    for k=2:ntFIT
        % set time step and excitation for current time step
        dt = timeFIT(k)-timeFIT(k-1);
        j = jifit(:,k);

        % leapfrog method
        [hfit(:,k),efit(:,k)] = leapfrog(hfit(:,k-1),efit(:,k-1),j,Mnu,Meps,Msigma,msh.C,dt);

        % voltage and current for input and output
        V1FIT(k) = efit(ipe2measureP1,k);                                      % left lower x-edge
        I1FIT(k) = (sum(abs(j))-V1FIT(k)/R);                                   % imposed current + conduction current
        V2FIT(k) = efit(ipe2measureP2,k);                                      % voltage at end of line
        I2FIT(k) = -V2FIT(k)/R;                                                % current at end of line

        % voltage along the transmission line (for each z-layer one value)
        VoiFIT(:,k) = efit(ipe2measure,k);
    end
    fprintf('finished field simulation using FIT after %d seconds.\n',toc(tstartFIT));

    % calculate error of V2
    V2CirInterpol = spline(timeCir,V2Cir,timeFIT);
    relerr = norm(V2CirInterpol-V2FIT)/norm(V2FIT);
    fprintf('relative error = %f\n',relerr);

    %% plot results

    % only plot results for large model, since this is the model used in the paper
    if strcmp(modelsize,'large')
        % plot Voi over z
        figure(22);
        switch termination
            case 'Reflect', subplot(2,1,1);
            case 'NoReflect', subplot(2,1,2);
        end
        hold all;
        plot(msh.z,VoiFIT(:,it2plotFIT));
        plot(msh.z,VoiCir(:,it2plotCir),'--');
        xlabel('z in m');
        ylabel('V_{oi}(z,t) in V');
        legend('t1','t2','t3','t4');
        switch termination
            case 'Reflect', title('total reflection case. FIT: solid, Circuit: dashed');
            case 'NoReflect', title('reflection-free case. FIT: solid, Circuit: dashed');
        end
        
        print(['./',modelname,modelsize,'Voi.pdf'],'-dpdf');

        % plot V2 over time
        figure(23);
        switch termination
            case 'Reflect', subplot(2,1,1);
            case 'NoReflect', subplot(2,1,2);
        end
        hold all;
        plot(1e9*timeFIT,V2FIT);
        plot(1e9*timeCir,V2Cir,'--');
        xlabel('time in ns');
        ylabel('V_{2}(t) in V');
        legend('FIT','Circuit');
        switch termination
            case 'Reflect', title('total reflection case. FIT: solid, Circuit: dashed');
            case 'NoReflect', title('reflection-free case. FIT: solid, Circuit: dashed');
        end
        print(['./',modelname,modelsize,'V2.pdf'],'-dpdf');
    end
    fprintf('finished absorbing boundary condition test case with %s termination after %d seconds.\n',termination,toc(tstart));
end