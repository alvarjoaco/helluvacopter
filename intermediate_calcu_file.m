%% Intermediate Calculation file

GrossInitial = GW; % just assigning initial GW for Rf_available calc.
Gw1 = GW(1); % just for ploting results no use in the code

area = (GW./DiskLoading)./nRotors; % ft^2, each main rotor area

if isCoaxial
    nRotors = 2.*nRotors;
    areaCoax = fudgeCoaxDisk.*area./2; % ft^2, individual rotor area
    R = sqrt(areaCoax./pi);
else
    R = sqrt(area./pi); % ft, main rotor radius
end
chord = (pi .* sigma.* R)./nBlades; % solidity
omega = Vt./R; % rad/s rotor rpm
f = 0.035.* GW.^0.67; % ft^2, flat plate area
[~, ~, ~, rho_sealevel] = atmoscoesa(convlength(0,'ft','m'));
rho_sealevel = convdensity(rho_sealevel, 'kg/m^3', 'slug/ft^3');


