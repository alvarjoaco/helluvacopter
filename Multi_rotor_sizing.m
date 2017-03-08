%% Sizing script
%{
Author:
 Karan Sutaria
 Last Modified: 2017-03-07 (yyyy-mm-dd)
 Version 2.7

    ALL INPUTS & CALCULATION DONE IN STANDARD ENGLISH UNITS

   Engine weight, specific fuel consumption, and drive train weight are
   based on curve fits scaled to work with the baseline specifications of
   the A-160 hummingbird. The structure weight fraction is also based on
   A-160 hummingbird and is assumed to be constant, barring advances in 
   airframe technology.

   Only the minimum information required need be entered for the mission
   segments. Unknows will be calculated from the given information. During
   climb and descent, the final altitude of the segment is used to
   calculate the power requirement over the entire segment. This
   overpredicts power required for climb and underpredicts power required
   for descent.

    Input paramters:
        - GW                            = Gross Weight ,            (lb)
        - DiskLoading                   = Disk Loading ,          (lb/ft^2)
        - sigma                         = Blade solidity ,
        - nBlades                       = Number of blades
        - Vt                            = Rotor Tip speed ,        (ft/s)
        - cd0                           = Cd0
        - sfc_baseline                  = sfc @ baseline engine, (lb/hp-hr)
        - structureWeightFraction       = Structure weight fraction
        - sfc                           = sfc for instantaneous segment**
        - engineWeightFraction          = Engine weight fraction**
        - driveSystemWeightFraction     = Drive system weight fraction**
        - f                             = flat plate area** ,       (ft^2)        
        - mcp_baseline                  = Baseline MCP ,            (hp)
        - wPayload                      = Payload ,                 (lb)
        - nCrew                         = number of Crew
        - wPerson                       = weight of the person ,    (lb)
        - nRotors                       = Number of Rotors++ 
        - isCoaxial                     = Coaxial Configuration = 1 
        - isDuct                        = if Duct = 1, no duct = 0
        
    Mission segments parameters:
        - Time              - hover, IGE, Ideal, climb, decent, warm up
        - Forward speed     - cruise i.e. max endurance speed
        - Altitude          - hover altitude
        - Distance          - length of each segment
        - Climb/Decent speed
        - Ideal engine time
        - In ground effect hover time 

    **NOTE:
        Flate plate area, engine weight fraction, and drive system 
        weight fracction are estimated based on curve fit equations, so 
        MUST MODIFY for specfic rotorcraft configuration.
    ++NOTE:
        IF even numbers of rotor used, then tail rotor does not reqire, so
        modify hpReq.m file accordingly for tail rotor. 
  
   Dependencies:
       hpReq.m
       mission_loading.m
       rf_required.m
       rf_available.m
       rf_method.m
       vehicle_define.m
       plot_results.m
       print_results.m
       intermediate_calc_file.m
       stWeight.m

    NOTES:
        1. Read hpReq.m Introduction for how that function works, and what
           assumptions were made to estimate power required. 
        2. The disk loading was equally distributed between all rotors.
        3. All rotors have same chord and radius.
        4. Each rotor had same solidity as of input.
        5. Disk loading was same for each rotor as of input loading.
        6. Read vehicle_define.m file
        7. Coaxial will add 28.08% power to single main rotor system and
           assume same rotor configuration as of top rotor.
        8. Duct will help reduce 15% power

%}
clear;
close all;
clc;

%% Inputs
% Draw plot for Rf method, Optimum design power consumption over different
% segments, Weight Change over segments for optimum design

isPlotting = true;

%% Vehicle Defining
vehicle_define

%%  Intermediate Calculation parameters
intermediate_calcu_file

%% Mission Definition
mission_loading

%% RF Method
rf_method

if designGw ~= Gw1
    % Plotting
    plot_results
    
    % Result Printing
    print_results
else
    printVehDefine
end
