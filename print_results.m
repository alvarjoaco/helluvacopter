%% Printing Results
%% Result Printing
fprintf('=========== Multi Rotor Rotorcraft ==========================\n')
fprintf('=========== Optimized Parameters ==============\n')
fprintf('Gross Weight                       = %.1f [lb]\n', designGw)
fprintf('Radius of main rotor               = %.2f  [ft]\n', Roptimum)
fprintf('Average chord                      = %.3f  [ft]\n', Coptimum)
fprintf('Rotor rpm                          = %.1f  [rpm]\n', omegaOptimum)
fprintf('MCP required                       = %.2f [hp]\n', designMcp)
fprintf('Take off Power                     = %.2f [hp]\n', takeoffP)
fprintf('Weighted Average Power Req.        = %.4f [hp]\n',...
  sum(table2array(tab1(:,3))'.*Power_consmp)./sum(table2array(tab1(:,3))));
fprintf('Figure of Merit                    = %.4f \n', FigMrt)
fprintf('Take off Blade loading             = %.4f \n', takeOffblLoad)
fprintf('Maximum Blade loading              = %.4f\n', designBladeLoading)
fprintf('Fuel Fraction required             = %.4f\n', Rf)
fprintf('Empty Weight fraction (Phi)        = %.4f\n', ...
    designEmptyWeightFraction)
fprintf('Engine weight fraction             = %.4f\n', ... 
    designEngineWeightFraction)
fprintf('Drive system weight fraction       = %.4f\n', ...
    designdriveSystemWeightFraction)
fprintf('Airframe structure fraction        = %.4f\n', ...
    designstructureWeightFraction)

fprintf('========== Input Parameters ================\n')
fprintf('Disk loading           = %.3f [lb/ft^2]\n', DiskLoading)
fprintf('sfc                    = %.3f [lb/hr-hp]\n',sfc_baseline)
fprintf('Number of blades       = %d \n',nBlades)
fprintf('Solidity               = %.3f \n', sigma)
fprintf('Cd0                    = %.3f \n',cd0)
fprintf('Number of Rotors       = %d \n', nRotors)
fprintf('Payload                = %d\n', wPayload)
if isCoaxial
    fprintf('Coaxial Rotor\n')
end
if isDuct
    fprintf('Ducted Fan Configuration\n')
end
tab1