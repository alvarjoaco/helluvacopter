%% Vehicle Defining

GW = 500:2:10000; % lb desire range of weight sizing analysis range ...
    ... for the gross weight

% Baseline input parameters

DiskLoading = 5;% lb/ft^2 ****
% ****NOTE: Disk loading will be each rotor's disk loading, i.e. single
% rotor or multi rotor or coaxial all's each rotor has same disk loading,
% so for multi and coaxial, area will be half the area of single rotor.
sigma = 0.08; %0.0757; % solidity, is also of each rotor
nBlades = 4;
Vt = 550; % ft/s rotor tip speed
cd0 = 0.01; 
sfc_baseline = 0.35; % lb/hp-hr
nRotors = 1; % keep at 1 for SMR & Coaxial
nCrew = 0; % fuck the crew
wPerson = 200; % lb, single crew weight
wPayload = 180; % lb, payload
isDuct = 0; % Duct is there = 1, no duct = 0 ++++
isCoaxial = 1; % is Coaxial = 1, no coaxial = 0 *** 
fudgeCoaxDisk = 1.05; %disk area fudge correction, add extra 5% area
miscellaneousWeightFraction = 0; % Miscellaneous Weight/ Gross Weight

% for sensitivity of Phi ONLY, otherwise code by itself do weight build up
structureWeightFraction = 0.3; % Structure Weight/ Gross Weight
driveSystemWeightFraction = 0.03; % Drive System Weight/ Gross Weight
engineWeightFraction = 0.05; % Engine Weight/ Gross Weight


%+++NOTE: Duct will help reduce 20% power very optimistic, but realistic 
%***NOTE: you need to add extra structure weight for coaxial & Duct, and
%   coaxial will increase power by 28.08% to single main rotor and has same
%   rotor configuration i.e chord, blades etc as of single rotor.

%% Mission Defining

baseAltitude = 0; % ft, helipad base altitude
warmUpTime = 2/60; % hour, warm up time
hoverTaxi = 1/16; % hour, hover taxi - Hover IGE - in ground effects
altitudeToHover = 500; % ft, hover altitude at each hover station
vrocForClimb = convvel(500,'ft/min','ft/s'); % ft/s, vertical rate of ...
 ... climb for take off.
segmentLength = 1; % nm, segment distance between each station in ... 
                       ... nautical miles
cruiseSpeed = 8; % kts, cruise speed in knots
hoverTimeAtStationA = 8; % hour, hover time at 1st station
hoverTimeAtStationB = 8; % hour, hover time at 2nd station
totalHoverTime = 24; % hour, total hour time
hoverSpeed = 0; % kts, hovering speed in the special case of drons ...
                ... making circles on one place
shutDownTime = 2/60; % hour, shut down time
maxSegmentTime = 0.1; % hour, maximum time allow before updating ... 
                        ... gross weight for fuel reduction. So, if any ...
                        ... segment is greater than 10 min, it will ...
                        ... decomose into sub segments for more accurate...
                        ... results.