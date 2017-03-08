%% Printing Vehicle Defining
fprintf('------------ Base-Line Values ------------\n')
fprintf('Gross Weight Increment (lbs.):  %d\n', GW(2)-GW(1)); 
fprintf('Disk Loading (lbs/ft^2):        %.2f\n', DiskLoading);
fprintf('Sigma (lbs/ft^2):               %.3f\n', sigma);
fprintf('Number of Blades:               %d\n', nBlades);
fprintf('Blade Tip Speed (ft/s):         %d\n', Vt);
fprintf('SFC (lbs/hr-hp):                %.3f\n', sfc_baseline);
fprintf('Total Empty Weight Fraction:    %.3f\n', designEmptyWeightFraction);
sprintf('Structure Weight phi:           %.3f\n', structureWeightFraction);
sprintf('Drive System Weight phi:        %.3f\n', driveSystemWeightFraction);
sprintf('Engine Weight phi:              %.3f\n', engineWeightFraction);
fprintf('Number of Rotors:               %d\n', nRotors);
fprintf('Payload Weight (lbs.):          %d\n', wPayload);
fprintf('Forward Flight Speed (kts.):    %d\n', cruiseSpeed);
fprintf('Hover Segment A (hr.):          %d\n', hoverTimeAtStationA);
fprintf('Hover Segment B (hr.):          %d\n', hoverTimeAtStationB);
fprintf('Hover Segment C (hr.):          %d\n', 24-hoverTimeAtStationA-hoverTimeAtStationB);
fprintf('-------------------------------------------\n')
if isDuct && isCoaxial
    fprintf('Ducted Coaxial Tail-Sitter Configuration\n');
elseif isDuct && ~isCoaxial
    fprintf('Ducted SMR Tail-Sitter Configuration\n');
elseif ~isDuct && isCoaxial
    fprintf('Coaxial Tail-Sitter Configuration\n');
else
    fprintf('SMR Tail-Sitter Configuration\n');
end