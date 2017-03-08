%% Nested Sensitivity Analysis to find optimum point
clear
close all
clc
    
 %% Optimization Parameters   
DiskLoadingArr = 3:1:5; % lb/ft^2
sigmaArr = .07:.01:0.09; % solidity
nRotorsArr = 1:1:2;
cd0Arr = 0.008:0.001:0.009;
sfc_baselineArr = .325:.005:0.35; % lb/hp-hr
isDuctArr = [0 1]; % Duct is there = 1, no duct = 0 ++++
isCoaxialArr = [0 1]; % is Coaxial =1, no coaxial = 0 *** 

%% Set Constraints
LimitBladeLoading = 0.14; % Max Blade Loading
LimitFigureMerit = 0.75; % min Fig of Merit
LimitPhi = 0.25; % Min Phi

%% Calculation
nDesigns = length(DiskLoadingArr).*length(nRotorsArr).*length(sigmaArr)...
    .*length(cd0Arr).*length(sfc_baselineArr).*length(isDuctArr).*...
    length(isCoaxialArr);

% Calculate minimum gross weights according to mission in sizing.m
fprintf('Sizing %d combinations using sizing.m ...\n\n', nDesigns')
designNumber = 1;
tic;
nondesignNumber = [];
for iisCoaxial = 1:length(isCoaxialArr)
    for iisDuct = 1:length(isDuctArr)
        for inRotors = 1:length(nRotorsArr)
            for isfc = 1:length(sfc_baselineArr)
                for icd0 = 1:length(cd0Arr)
                    for isigma = 1:length(sigmaArr)
                        for iDL = 1:length(DiskLoadingArr)
                            % saving design concept
                            design(designNumber).isCoaxial = isCoaxialArr(iisCoaxial);
                            design(designNumber).isDuct = isDuctArr(iisDuct);
                            design(designNumber).nRotors = nRotorsArr(inRotors);
                            design(designNumber).sfc = sfc_baselineArr(isfc);
                            design(designNumber).cd0 = cd0Arr(icd0);
                            design(designNumber).sigma = sigmaArr(isigma);
                            design(designNumber).DiskLoading = DiskLoadingArr(iDL);
                            
                            % Vehicle Define & Variable
                            vehicle_define
                            
                            % Initializing design
                            isCoaxial = isCoaxialArr(iisCoaxial);
                            isDuct = isDuctArr(iisDuct);
                            nRotors = nRotorsArr(inRotors);
                            sfc_baseline = sfc_baselineArr(isfc);
                            cd0 = cd0Arr(icd0);
                            sigma = sigmaArr(isigma);
                            DiskLoading = DiskLoadingArr(iDL);
                            
                            %  Intermediate Calculation parameters
                            intermediate_calcu_file
                            
                            % Mission Definition
                            mission_loading
                            
                            % RF Method
                            rf_method
                            
                            % Ignore non-converge solution
                            if designGw == Gw1
                                nondesignNumber(end+1) = designNumber;
                                designGw = 1000000;
                            end
                            bladeLoading = max(takeOffblLoad,designBladeLoading);
                            
                            % Saving GW
                            design(designNumber).GrossWeight = designGw;
                            design(designNumber).BladeLoading = bladeLoading;
                            design(designNumber).FigureMerit = FigMrt;
                            design(designNumber).phi = designEmptyWeightFraction;
                            design(designNumber).Radius = Roptimum;
                            designNumber = designNumber + 1;
                        end
                    end
                end
            end
        end
    end
end

time = toc;
fprintf('Sized %d combinations in %0.2f min.\n\n', nDesigns, time./60)
fprintf('%d combinations of designs did NOT work\n\n',...
    length(nondesignNumber))

%% Optimum design Selection
range = 1:designNumber;
rangeAcceptable = range([design.BladeLoading] <= LimitBladeLoading);
rangeAcceptable = range([design(rangeAcceptable).FigureMerit] >= LimitFigureMerit);
rangeAcceptable = range([design(rangeAcceptable).phi] >= LimitPhi);
fprintf(['%d combinations pass all constraint functions.\n\n'], ...
    length(rangeAcceptable))

% Select optimal design
[~, index] = min([design(rangeAcceptable).GrossWeight]);
optDesign = design(rangeAcceptable(index));

%% Run Sizing on optimul Design Configuration

vehicle_define

% Initializing design
isCoaxial = optDesign.isCoaxial;
isDuct = optDesign.isDuct;
nRotors = optDesign.nRotors;
sfc_baseline = optDesign.sfc;
cd0 = optDesign.cd0;
sigma = optDesign.sigma;
DiskLoading = optDesign.DiskLoading;

%  Intermediate Calculation parameters
intermediate_calcu_file

% Mission Definition
mission_loading

% RF Method
rf_method

% Plotting
isPlotting = true;
plot_results

% Result Printing
print_results
