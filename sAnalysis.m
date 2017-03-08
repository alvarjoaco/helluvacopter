%% Sensitivity Analysis
clear
close all
clc 

%% Analysis Set Up
% set to 1 if you want to run this sensitivity analysis for the parameters
isSol = 1; % solidity
isPhi = 0; % empty weight fraction keep 0
isDL = 1; % disk loading
isSFC = 1; % sfc
isCDO = 1; % drag coeff
isNR = 1; % number of rotors
isDut = 0; % Duct or No Duct
isCox = 1; % Coaxial or No Coaxial
isPlot = 1; % plot on or off
isSubP = 0; % individual plots on or off
res = 70; % samples to probe or resolution
SetUp = 'vehicle_define.m';

%% Sensitivity Parameters
% min:increments:max

% Momemtum parameters
DiskLoadingArr = linspace(1,10,res); % lb/ft^2
sigmaArr = linspace(.05,.15,res); % solidity
nRotorsArr = 1:1:4;
cd0Arr = linspace(.007,.01,res);
sfc_baselineArr = linspace(.25,.55,res); % lb/hp-hr
structureWeightFractionArr = linspace(.22,.4,res);%.285; % Structure Weight/ Gross Weight
driveSystemWeightFractionArr = linspace(.02,.1,res); % Drive System Weight/ Gross Weight
engineWeightFractionArr = linspace(.03,.1,res); %0.05; % Engine Weight/ Gross Weight

%CBEMT parameters
nBladesArr = 2:1:8;
VtArr = 450:10:800; % ft/s rotor tip speed
cruiseSpeedArr = 1:2:50; % kts, cruise speed in knots
hoverAArr = [22, 1, 8, 1]; % hour, hover time at 1st station
hoverBArr = [1, 22, 8, 1]; % hour, hover time at 2nd station
hoverCArr = 24 - hoverAArr - hoverBArr;

% Cofiguration parameters
isDuctArr = [0 1]; % Duct is there = 1, no duct = 0 ++++
isCoaxialArr = [0 1]; % is Coaxial =1, no coaxial = 0 *** 
fudgeCoaxDiskArr = 1:.01:1.15;%1.05; %disk area fudge ##change if isCoaxial = 1##

time = tic;
%% Disk Loading Sensitivity
tic
designGWArr = zeros(1,res);
if isDL
    
    for a = 1:length(DiskLoadingArr)
        %% Vehicle Define & Variable
        run(SetUp)
        DiskLoading = DiskLoadingArr(a);

        %%  Intermediate Calculation parameters
        intermediate_calcu_file

        %% Mission Definition
        mission_loading

        %% RF Method
        rf_method

        %% Ignore non-converge solution
        if designGw == Gw1
            DiskLoadingArr(a);
            designGw = NaN;
        end

        %% Saving GW
        designGWArr(a) = designGw;
    end
    if isSubP
        figure
        plot(DiskLoadingArr,designGWArr)
        title('Disk Loading vs. Gross Weight')
        xlabel('Disk Loading (lb/ft^2)')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_DL, locMinDL] = min(designGWArr);
    minDL = DiskLoadingArr(locMinDL);
    DiskLoading = NaN;
    fprintf('=========== Disk Loading Sensitivity Analysis ===========\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal Disk Loading (lbs/ft^2): %.2f\n', minDL)
    fprintf('Optimal DL GW (lbs):             %d\n',minGW_DL)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_dl = designGWArr;
else
    gwArr_dl = zeros(1,res);
end

%%  SFC Sensitivity
tic
designGWArr = zeros(1,res);
if isSFC
    
    for a = 1:length(sfc_baselineArr)
        %% Vehicle Define & Variable
        run(SetUp)
        sfc_baseline = sfc_baselineArr(a);

        %%  Intermediate Calculation parameters
        intermediate_calcu_file

        %% Mission Definition
        mission_loading

        %% RF Method
        rf_method

        %% Ignore non-converge solution
        if designGw == Gw1
            sfc_baselineArr(a);
            designGw = NaN;
        end
        %% Saving GW
        designGWArr(a) = designGw;
    end
    if isSubP
        figure
        plot(sfc_baselineArr,designGWArr)
        title('SFC vs. Gross Weight')
        xlabel('SFC (lb/hp-hr)')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_SFC, locMinSFC] = min(designGWArr);
    minSFC = sfc_baselineArr(locMinSFC);
    sfc_baseline = NaN;
    fprintf('=============== SFC Sensitivity Analysis ===============\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal SFC (lbs/hp-hr):     %.3f\n', minSFC)
    fprintf('Optimal SFC GW (lbs):        %d\n',minGW_SFC)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_sfc = designGWArr;
else
    gwArr_sfc = zeros(1,res);
end

%%  Number of Rotor Sensitivity
tic
designGWArr = zeros(1,length(nRotorsArr));
if isNR
    for a = 1:length(nRotorsArr)
        %% Vehicle Define & Variable
        run(SetUp)
        nRotors = nRotorsArr(a);

        %%  Intermediate Calculation parameters
        intermediate_calcu_file

        %% Mission Definition
        mission_loading

        %% RF Method
        rf_method

        %% Ignore non-converge solution
        if designGw == Gw1
            nRotorsArr(a);
            designGw = NaN;
        end
        %% Saving GW
       designGWArr(a) = designGw;

end
    if isSubP
        figure
        plot(nRotorsArr,designGWArr)
        title('Number of Rotors vs. Gross Weight')
        xlabel('Number of Rotors')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_NR, locMinNR] = min(designGWArr);
    minNR = nRotorsArr(locMinNR);
    nRotors = NaN;
    fprintf('=============== Num. Rotors Sensitivity Analysis =============\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal Number of Rotors:         %d\n', minNR)
    fprintf('Optimal Num. Rotors GW (lbs):     %d\n',minGW_NR)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_nr = designGWArr;
else
    gwArr_nr = zeros(1,length(nRotorsArr));
end

%% Solidity Sensitivity
tic
designGWArr = zeros(1,res);
if isSol
    for a = 1:length(sigmaArr)
        %% Vehicle Define & Variable
        run(SetUp)
        sigma = sigmaArr(a);

        %%  Intermediate Calculation parameters
        intermediate_calcu_file

        %% Mission Definition
        mission_loading

        %% RF Method
        rf_method

        %% Ignore non-converge solution
        if designGw == Gw1
            sigmaArr(a);
            designGw = NaN;
        end
        %% Saving GW
       designGWArr(a) = designGw;
    end
    if isSubP
        figure
        plot(sigmaArr,designGWArr)
        title('Solidity vs. Gross Weight')
        xlabel('Solidity (lb/ft^2)')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_Sigma, locMinSigma] = min(designGWArr);
    minSigma = sigmaArr(locMinSigma);
    sigma = NaN;
    fprintf('=============== Solidity Sensitivity Analysis ===============\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal Sigma (lbs/ft^2):       %.3f\n', minSigma)
    fprintf('Optimal Sigma GW (lbs):         %d\n',minGW_Sigma)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_sigma = designGWArr;
else
    gwArr_sigma = zeros(1,res);
end

%%  Empty Weight Sensitivity
tic
designGWArr = zeros(1,res);
if isPhi
    emptyWF = [];
    for a = 1:length(structureWeightFractionArr)
        % Vehicle Define & Variable
        run(SetUp)
        
        %  Intermediate Calculation parameters
        intermediate_calcu_file
        
        % Mission Definition
        mission_loading
        
        % Rf Required Calculation
        rf_required
        
        % For sensitivity phi change
        structureWeightFraction = structureWeightFractionArr(a);
        engineWeightFraction = engineWeightFractionArr(a);
        driveSystemWeightFraction = driveSystemWeightFractionArr(a);
        emptyWf = structureWeightFraction...
            +engineWeightFraction + driveSystemWeightFraction;
        emptyWF(end+1) = emptyWf;
        
        % Rf Available Calculation
        wCrew = nCrew .* wPerson;
        crewFraction = wCrew ./ GW;
        payloadFraction = wPayload ./ GW;
        emptyWeightFraction = structureWeightFraction ...
            + engineWeightFraction + driveSystemWeightFraction ...
            + miscellaneousWeightFraction;
        
        fuelFractionAvail = 1 - emptyWeightFraction - crewFraction...
            - payloadFraction;
        
        % Rf Method
        % Picks minimum GW with enough fuel fraction
        differentFuelFraction = fuelFractionAvail - fuelFractionReq;
        [~, indi] = find(differentFuelFraction>=0);
        if isempty(indi)
            warning('There is NO converge Gross Weight for a given configuration');
            indi = 1;
        end
        baseMin = indi(1);
        Rf = fuelFractionReq(baseMin);
        designGw = GW(baseMin);
        
        % Ignore non-converge solution
        if designGw == Gw1
            emptyWF(end);
            designGw = NaN;
        end
        % Saving GW
       designGWArr(a) = designGw;
    end
    if isSubP
        figure
        plot(emptyWF,designGWArr)
        title('Empty Weight Fraction vs. Gross Weight')
        xlabel('Empty Weight Fraction')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_WF, locMinWF] = min(designGWArr);
    minWF = emptyWF(locMinWF);
    engineWeightFraction = NaN;
    driveSystemWeightFraction = NaN;
    structureWeightFraction = NaN;
    emptyWeightFraction = NaN;
    fprintf('=============== Phi Sensitivity Analysis ===============\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal Empty Weight Fraction :             %.3f\n', minWF)
    fprintf('Optimal Empty Weight Fraction GW (lbs):     %d\n',minGW_WF)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_phi = designGWArr;
else
    gwArr_phi = zeros(1,res);
end

%% Cd0 Sensitivity
tic
designGWArr = zeros(1,res);
if isCDO
    for a = 1:length(cd0Arr)
        %% Vehicle Define & Variable
        run(SetUp)
        cd0 = cd0Arr(a);

        %%  Intermediate Calculation parameters
        intermediate_calcu_file

        %% Mission Definition
        mission_loading

        %% RF Method
        rf_method

        %% Ignore non-converge solution
        if designGw == Gw1
            cd0Arr(a);
            designGw = NaN;
        end
        %% Saving GW
        designGWArr(a) = designGw;
    end
    if isSubP
        figure
        plot(cd0Arr,designGWArr)
        title('Cd_0 vs. Gross Weight')
        xlabel('Cd_0')
        ylabel('Design Gross Weight (lbs)')
    end
    [minGW_cd0, locMincd0] = min(designGWArr);
    mincd0 = cd0Arr(locMincd0);
    cd0 = NaN;
    fprintf('=============== Cd0 Sensitivity Analysis ===============\n')
    printVehDefine
    fprintf('----------- Optimized Parameters ----------\n')
    fprintf('Optimal Cd0:       %.3f\n', mincd0)
    fprintf('Optimal Cd0 GW (lbs):         %d\n',minGW_cd0)
    fprintf('-------------------------------------------\n')
    toc
    fprintf('======================== End ===========================\n\n')
    gwArr_cd0 = designGWArr;
else
    gwArr_cd0 = zeros(1,res);
end
  

%% All Sensitivity on 1 plot
if isPlot
    figure
    dl_a = linspace(0,1,length(DiskLoadingArr));
    sg_a = linspace(0,1,length(sigmaArr));
    nr_a =  linspace(0,1,length(nRotorsArr));
    cd_a =  linspace(0,1,length(cd0Arr));
    sfc_a = linspace(0,1,length(sfc_baselineArr));
    ph_a =   linspace(0,1,length(engineWeightFractionArr));
    plot(dl_a,gwArr_dl,sg_a,gwArr_sigma,nr_a,gwArr_nr,...
        cd_a,gwArr_cd0,sfc_a,gwArr_sfc,ph_a,gwArr_phi)
    xlabel('Low to High: Parameter Values')
    xlim([0 1])
    ylabel('Gross Weight (lb)')
    title('Sensitivity Analysis')
    legend('Disk Loading','Solidity','# of Rotors','Cd_0','SFC','\phi'...
        ,'Location','best')
end

fprintf('================ Sensitivity Analysis Completed ===============\n\n')

toc(time)
