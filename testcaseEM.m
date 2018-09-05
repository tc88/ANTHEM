function testcaseEM(modelsize,LTspiceEXE,wineEXE,versionLTspice,verbose)
% TESTCASEEM implements a cuboidal cavity with centred excitation that can
% either be a voltage or a current. The resulting TE and TM modes are
% compared with the analytical reference solution
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

modelname = 'testcaseEM';                                                  % name of the model that is also used for output files

% settings
problemSet = {'TE','TM'};
forceNewFITsim = false;   % as a default, FIT     simulation is only done if .raw-file does not exist. To force new FIT     simulation, set this parameter to true
forceNewCirSim = false;   % as a default, circuit simulation is only done if .raw-file does not exist. To force new circuit simulation, set this parameter to true

for i = 1:length(problemSet)

    problem = problemSet{i};
    tstart = tic;
    fprintf('running electromagnetic %s test case...\n',problem);
    fileCir = ['./',modelname,'results',problem,modelsize,'Cir'];
    fileFit = ['./',modelname,'results',problem,modelsize,'Fit'];

    % sanitary checks
    if ~strcmp(problem,'TM') && ~strcmp(problem,'TE')
        error('wrong problem type defined, should be ''TE'' or ''TM''!');
    end

    % load problem
    load(['./',modelname,modelsize,'.mat'])

    %% solve eigenvalue problem

    % analytical solution for eigenfrequencies
    c = 1/sqrt(mu0*mur*eps0*epsr);
    fAna = @(m,n,p) 0.5*c*sqrt((m/a)^2+(n/b)^2+(p/d)^2);

    % choose modes to calculate
    modes2calc = {'011','110','111','121','013','122'};

    fr = zeros(length(modes2calc),1);
    frEfit = zeros(length(modes2calc),1);
    for i = 1:length(modes2calc)
        currentMode = modes2calc{i};

        % calculate analytical resonant frequencies fr and expected eigenvalue
        m = str2double(currentMode(1));
        n = str2double(currentMode(2));
        p = str2double(currentMode(3));
        fr(i) = fAna(m,n,p);
        lambdaExpect = 4*pi^2*fr(i)^2;

        % find edges that are not PEC nor ghost edges and thus the inner edges
        ipeInner = find(diag(Meps)~=0);

        % solve generalized eigenvalue problem
        A = msh.C'*Mnu*msh.C;
        A = A(ipeInner,ipeInner);
        B = Meps(ipeInner,ipeInner);
        lambda1 = eigs(A,B,5,lambdaExpect);
        [~,idx2pick1] = min(abs(lambda1-lambdaExpect));
        frEfit(i) = sqrt(lambda1(idx2pick1))/2/pi;
    end

    %% solve excited problem using FIT
    fstart = 5e8; fstop = 3e9;

    % place excitation current js in center (single edge for TM, loop for TE)
    ipnCenter = idx2canonical(msh,ceil(msh.nx/2),ceil(msh.ny/2),ceil(msh.nz/2));
    js = zeros(3*msh.np,1);
    switch problem
        case 'TM', Nfreq = 2000;
                   ipeEx = 2*msh.np + ipnCenter;
                   js(ipeEx) = 1;
        case 'TE', Nfreq = 3000;
                   ipeEx = [ipnCenter       ,msh.np+ipnCenter+msh.Mx, ...
                            ipnCenter+msh.My,msh.np+ipnCenter];
                   js(ipeEx(1:2)) = -1;
                   js(ipeEx(3:4)) = 1;
    end
    omegaVec = 2*pi*linspace(5e8,3e9,Nfreq);
    freqFit = omegaVec/2/pi;
    % solve electromagnetic problem using field solver
    if ~exist(['./' fileFit '.mat'],'file') || forceNewFITsim
        tstartFIT = tic;
        fprintf('running field simulation using FIT ...\n');
        efit = zeros(3*msh.np,length(omegaVec));
        for i = 1:length(omegaVec)
            omega = omegaVec(i);
            if verbose(1), fprintf('solving for f=%e...\n',omega/2/pi); end
            Aex = msh.C'*Mnu*msh.C - omega^2*Meps;
            b = -1j*omega*js;
            efit(ipeInner,i) = Aex(ipeInner,ipeInner)\b(ipeInner);
        end
        save([fileFit,'.mat'],'N','efit');
        fprintf('finished field simulation using FIT after %d seconds.\n',toc(tstartFIT));
    else
        warning('FIT results found, FIT simulation skipped!');
    end

    %% create circuit and solve excited problem using circuit solver

    % use current source as excitation in circuit model
    Isetting.pos = ipeEx;
    if     strcmp(problem,'TM'), Isetting.type = 'e'; Isetting.amp = 1; options.np = Nfreq;
    elseif strcmp(problem,'TE'), Isetting.type = 'e'; Isetting.amp = [-1 -1 1 1]; options.np = 3000;
    end

    % extract netlist from electromagnetic field model
    analysis = 'AC';
    options.fstart = fstart;
    options.fstop = fstop;
    options.type = 'dec';
    extractNetlist4FITEM(fileCir,msh,Meps,Msigma,Mnu,Isetting,analysis,msh.ipePEC,options.fstart,options.fstop,options.np,options.type);

    % solve electromagnetic problem using circuit solver
    if ~exist([fileCir,'.raw'],'file') || forceNewCirSim
        runLTspice([fileCir,'.cir'],'-b',LTspiceEXE,wineEXE);
    else
        warning('Circuit results found, circuit simulation skipped!')
    end

    %% plot results

    % only plot results for large model, since this is the model used in the paper
    if strcmp(modelsize,'large')

        % load results from file
        if ~exist('efit','var'), load([fileFit '.mat']); end
        if ~exist('spiceData','var'), spiceData = LTspice2Matlab([fileCir '.raw'],'all',1,versionLTspice); end

        % choose edge to plot depending on the problem
        switch problem
            case 'TM', ipe2plot = 6*msh.Mx + 5*msh.My + 5*msh.Mz;
                       ipe2save = ipe2plot + 2*msh.np;
            case 'TE', ipe2plot = 6*msh.Mx + 3*msh.My + 5*msh.Mz;
                       ipe2save = ipe2plot;
        end

        % get voltages in dB
        freqCir = spiceData.freq_vect;
        eCir = spiceOrderOutput(spiceData);
        eFitInDeb = 20*log10(abs(efit));
        eCirInDeb = 20*log10(abs(eCir));

        % identify resonant frequencies and calculate error between fit and circuit
        if strcmp(problem,'TE')
            figure(20); subplot(1,2,1); hold all;
            plot(freqFit,eFitInDeb(ipe2save,:));
            plot(freqCir,eCirInDeb(ipe2save,:),'--');
            legend('FIT','Circuit');
            title('TE mode');
            print(['./',modelname,'results',upper(modelsize(1)),modelsize(2:end),'.pdf'],'-dpdf');

            if ~isOctave
                [~,frTEfit] = findpeaks(eFitInDeb(ipe2save,:),freqFit);
                [~,frTEcir] = findpeaks(eCirInDeb(ipe2save,:),freqCir);
                epsrTE = abs(frTEfit(1:end-1)-frTEcir)./frTEfit(1:end-1);
                fprintf('epsrTE = \n');
                display(epsrTE(1:6)');
            end
        elseif strcmp(problem,'TM')
            figure(20); subplot(1,2,2); hold all;
            plot(freqFit,eFitInDeb(ipe2save,:));
            plot(freqCir,eCirInDeb(ipe2save,:),'--');
            legend('FIT','Circuit');
            title('TM mode');
            print(['./',modelname,'results',upper(modelsize(1)),modelsize(2:end),'.pdf'],'-dpdf');

            if ~isOctave
                [~,frTMfit] = findpeaks(eFitInDeb(ipe2save,:),freqFit);
                [~,frTMcir] = findpeaks(eCirInDeb(ipe2save,:),freqCir);
                epsrTM = abs(frTMfit-frTMcir)./frTMfit;
                fprintf('epsrTM = \n');
                display(epsrTM(1:4)');
            end
        end
    end
    
    fprintf('finished electromagnetic %s test case after %d seconds...\n',problem,toc(tstart));
    clear efit spiceData
end