%% Rf available calculation

wCrew = nCrew .* wPerson;
crewFraction = wCrew ./ GW;
payloadFraction = wPayload ./ GW;

[structureWeightFraction,engineWeightFraction,driveSystemWeightFraction]...
    = stWeight(R,nBlades,chord,GrossInitial,isDuct,maxContPowerReq,...
    nRotors, isCoaxial); 

emptyWeightFraction = structureWeightFraction ...
    + engineWeightFraction + driveSystemWeightFraction ... 
    + miscellaneousWeightFraction; 

fuelFractionAvail = 1 - emptyWeightFraction - crewFraction...
    - payloadFraction;
