%% Rf method And finding optimum parameters

% Rf Required Calculation

rf_required

% Rf Available Calculation

rf_available



% Picks minimum GW with enough fuel fraction
differentFuelFraction = fuelFractionAvail - fuelFractionReq;
[~, indi] = find(differentFuelFraction>=0);
%[~, indi] = min(abs(differentFuelFraction));
if isempty(indi)
    %warning('There is NO converge Gross Weight for a given configuration')
    indi = 1;
end
baseMin = indi(1);
Rf = fuelFractionReq(baseMin);
designGw = GW(baseMin);

% Optimum parameters
designMcp = maxContPowerReq(baseMin);
designBladeLoading = maxBladeLoading(baseMin);
[~,colDesgn] = size(emptyWeightFraction);
if colDesgn < baseMin
    designstructureWeightFraction = structureWeightFraction;
    designEmptyWeightFraction = emptyWeightFraction;
    designEngineWeightFraction = engineWeightFraction;
    designdriveSystemWeightFraction = driveSystemWeightFraction;
else
    designEngineWeightFraction = engineWeightFraction(baseMin);
    designdriveSystemWeightFraction = driveSystemWeightFraction(baseMin);
    designstructureWeightFraction = structureWeightFraction(baseMin);
    designEmptyWeightFraction = emptyWeightFraction(baseMin);
end
    
Roptimum = R(baseMin);
Coptimum = chord(baseMin);
omegaOptimum = omega(baseMin);
omegaOptimum = convangvel(omegaOptimum,'rad/s','rpm');

% Figure of Merit of optimum configuration
f = 0.035.* designGw.^0.67; % ft^2, flat plate area
h_sea = 0; % ft seal level altitutde 
V = 0; % ft/s
[takeoffP, takeOffblLoad,FigMrt] = hpReq(h_sea, designGw, Roptimum, ... 
                            sigma, Vt, cd0,f, V,nRotors,isCoaxial,isDuct);
