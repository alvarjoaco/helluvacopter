function [hpTotal, blLoad,FigMrt] = hpReq(h, GW, R, sigma, VT, ...
                              cd0,f, V,nRotors,isCoaxial,isDuct, varargin)
%{
Author:
 Karan Sutaria
 Last Modified: 2017-03-07 (yyyy-mm-dd)

    ALL INPUTS & CALCULATION DONE IN STANDARD ENGLISH UNITS

   Input variables:
    h           = altitude,                                         (ft)
    GW          = gross weight,                                     (lb)
    R           = radius of the blade,                              (ft)
    sigma       = blade solidity
    VT          = Rotor tip speed                                   (ft/s)
    cd0         = airfoil's cd0,
    f           = flat plate area,                                  (ft^2)
    V           = forward speed,                                    (ft/s)
    nRotors     = number of rotors    
    isCoaxial   = Coaxial or not i.e. 1 or 0
    isDuct      = Duct or not i.e. 1 or 0
                                            
   Optional variables:
    varagin = [Vc, isIdle, isIGE] in that order variables,
    Vc      = climb speed,                                      (ft/s)
    isIdle  = ideal 5% power required of engine condition,      
    isIGE   = in ground effects 10% power required condition,   

   Output variable:
    hpReq   = power required for one instance of the segment,   (hp)
    blLoad  = blade loading for each segment
    FigMrt  = figure of merit based on hover flight condition at sea level
    
   
    NOTES: 
        1. For climb Vc = +Vc, and for descent Vc = -Vc already assumed as
            given in the inputs.
        2. Momentum theory is used to estimate power requiared with 
            correction factors. 
        
        Cp = k*Ct*lambda_h*Ku+sigma*cd0/8*(1+4.65*mu^2)+f.*mu.^3./(2.*A).

        3. For tail rotor, in ground effects, and ideal condition
            correction factors were applied.
        4. Figure of merit equation:
                                  Power Ideal(Ct*lambda_h)
    Figure of Merit =   ------------------------------------------
                          Total Power (  k*Ct*lambda_h*Ku+sigma*cd0/8 )
%}


%% Assumptions 

% Power correction factors
tailRotorFactor = .1; % 10% power consumed by counter-torque, tailRotorFactor = 0.1
groundEffectFactor = 0.9; % for IGE, only 90% of power required.
idleFactor = 0.05; % for Ideal Engine only 5% of power required
kappa = 1.15;
ductFactor = 0.8; % 20% power saving
CoaxialFactor = 0.22; % 28.08% power increment to single rotor.

if isCoaxial
    tailRotorFactor = 0;
    ductFactor = 0.90;
end

if mod(nRotors,2) == 0
    tailRotorFactor = 0;
end
%% Calculation 

% Check number of inputs
numvarargs = length(varargin);
if numvarargs > 3
    error('hpReq:TooManyInputs', ...
        'requires at most 3 optional inputs');
end

% Set defaults for optional inputs
optargs = {0, false,false};
optargs(1:numvarargs) = varargin;
[Vc, isIdle,isIGE] = optargs{:};
A = pi .* R.^2; % ft^2 for one rotor
[~, ~, ~, rho] = atmoscoesa(convlength(h, 'ft','m'));
rho = convdensity(rho, 'kg/m^3', 'slug/ft^3');

% Power calculation

D = 0.5 .* rho .* V.^2 .* f; % parasite drag
a_tpp = atan(D ./ GW);
CT = GW ./ (rho .* A .* VT.^2 .*cos(a_tpp) .*nRotors); % for each rotor Ct
blLoad = CT./sigma; % blade loading for each segment for each GW
lambda_h = sqrt(CT ./ 2); % hover induced inflow ratio
mu = V .* cos(a_tpp) ./ VT; % advance ratio
% Glauert inflow correction factor
Ku = sqrt((-(mu ./ lambda_h).^2 + sqrt((mu ./ lambda_h).^4 +4)) ./ 2); 
% power coeffs. for each rotor
CPi = kappa .* CT .* lambda_h .* Ku; % induced power
CP0 = sigma .* cd0 ./ 8 .*(1 + 4.65 .* mu.^2); % profile power
CPp = f .*mu.^3 ./ (2 .* A); % parasite power
%CPp = D .*V ./ (rho .* A .* VT.^3) % parasite power
CP =  (CPi + CP0 + CPp); % total power equired by each rotor

hpTotal = CP .*(rho .* A .* VT.^3) ./ 550 ; % hp, of all rotors


% Power for n Rotors
hpTotal = nRotors.*hpTotal;

if isCoaxial
    hpTotal = hpTotal .* (1 + CoaxialFactor);
end

if isDuct
    hpTotal = ductFactor.*hpTotal;
end

% Tail Rotor include
hpTotal = hpTotal .* (1 + tailRotorFactor); % total power

% Climb power required
hpClimb = Vc .* GW ./ 550; % hp main rotor power with climb
if h == 0
    hpClimb = 0;
end
hpTotal = hpTotal + hpClimb;

% HIGE 
if isIGE
    hpTotal = hpTotal .* groundEffectFactor;
end

% Idle
if isIdle
    hpTotal = hpTotal .* idleFactor;
end

% Figure of Merit calculation
CT = GW ./ (rho .* A .* VT.^2);
lambda_h = sqrt(CT ./ 2); % hover induced inflow ratio
Pideal = CT.*lambda_h.* (rho .* A .* VT.^3) ./ 550 ; % hp
Pprofile = sigma.*(cd0./8).*(rho .* A .* VT.^3) ./ 550 ; % hp
Pactual = kappa.*Pideal + Pprofile;
FigMrt = Pideal./Pactual;

end

