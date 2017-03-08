%%  Rf-required calculation

for iSegment = 1:length(mission.segment)
    
    %  for j = 1:2 % WARNING: need to study why Max put this here
    % Distance = speed * time
    % Calculate segment unknowns - time,and speed if distance known
    if ~isempty(mission.segment(iSegment).distance)
        if ~isempty(mission.segment(iSegment).speed) ...
                && isempty(mission.segment(iSegment).time)
            mission.segment(iSegment).time = ...
                mission.segment(iSegment).distance ...
                ./ mission.segment(iSegment).speed ...
                ./ 3600;
        elseif ~isempty(mission.segment(iSegment).time) ...
                && isempty(mission.segment(iSegment).speed)
            mission.segment(iSegment).speed = ...
                mission.segment(iSegment).distance ...
                ./ mission.segment(iSegment).time ...
                .* 3600;
        end
    end
    % Altitude change = rate of climb * time
    % Calculation of time requried to acheive desire altitude
    if ~isempty(mission.segment(iSegment).altitudeChange)
        if ~isempty(mission.segment(iSegment).roc) ...
                && isempty(mission.segment(iSegment).time)
            mission.segment(iSegment).time = ...
                mission.segment(iSegment).altitudeChange ...
                ./ mission.segment(iSegment).roc ...
                ./ 3600;
        elseif ~isempty(mission.segment(iSegment).time) ...
                && isempty(mission.segment(iSegment).roc)
            mission.segment(iSegment).roc = ...
                mission.segment(iSegment).altitudeChange ...
                ./ mission.segment(iSegment).time ...
                .* 3600;
        end
    end
    % end % WARNING: need to study why Max put this here , for loop warning.
    
    % Set distance if none is specified
    if isempty(mission.segment(iSegment).distance)
        mission.segment(iSegment).distance = 0;
    end
    % Set velocity if none is specified
    if isempty(mission.segment(iSegment).speed)
        mission.segment(iSegment).speed = 0;
    end
    % Set altitude change if none is specified
    if isempty(mission.segment(iSegment).altitudeChange)
        mission.segment(iSegment).altitudeChange = 0;
    end
    % Set roc if none is specified
    if isempty(mission.segment(iSegment).roc)
        mission.segment(iSegment).roc = 0;
    end
    % Set time if none is specified
    if isempty(mission.segment(iSegment).time)
        mission.segment(iSegment).time = 0;
    end
    % Set isIdle if none is specified
    if isfield(mission.segment(iSegment), 'isIdle')
        if isempty(mission.segment(iSegment).isIdle)
            mission.segment(iSegment).isIdle = false;
        end
    else
        mission.segment(iSegment).isIdle = false;
    end
    % Set isIGE if none is specified
    if isfield(mission.segment(iSegment), 'isIGE')
        if isempty(mission.segment(iSegment).isIGE)
            mission.segment(iSegment).isIGE = false;
        end
    else
        mission.segment(iSegment).isIGE = false;
    end
    
    if iSegment ~=1
        % Calculate altitude to account for altitude change
        mission.segment(iSegment).altitude = ...
            mission.segment(iSegment - 1).altitude ...
            + mission.segment(iSegment).altitudeChange;
    end
    % Update gross weight to account for fuel burned in previous segment
    if iSegment == 1
        mission.segment(iSegment).gw = GW;
        
    else
        mission.segment(iSegment).gw = ...
            mission.segment(iSegment - 1).gw ...
            - mission.segment(iSegment - 1).wFuel;
    end
    
    % For small time segment updates
    if mission.segment(iSegment).time >= maxSegmentTime
        numIteration = floor(mission.segment(iSegment).time ...
            ./ maxSegmentTime);
        lastIterTime = mission.segment(iSegment).time - numIteration ...
            .*maxSegmentTime;
        fuelIncrement = [];
        hpReqIncrement = [];
        blLoadIncrement = [];
        gwIncrement = mission.segment(iSegment).gw;
        for ii = 1:numIteration
            
            % Calculate hp for the segment 
            [hpReqIncrement(end+1,:),blLoadIncrement(end+1,:)]...
                = hpReq(mission.segment(iSegment).altitude, ...
                gwIncrement, R, sigma, Vt, cd0, f, ...
                mission.segment(iSegment).speed, nRotors,isCoaxial,...
                isDuct, mission.segment(iSegment).roc, ...
                mission.segment(iSegment).isIdle,...
                mission.segment(iSegment).isIGE);
            % Calculate fuel 
            fuelIncrement(end+1,:) = ...
                sfc_baseline .* ...
                hpReqIncrement(end,:) .* maxSegmentTime;
            gwIncrement = gwIncrement - fuelIncrement(end,:);
        end
        if lastIterTime > 0
            
            % Calculate hp for the segment's last leftover time
            [hpReqIncrement(end+1,:),blLoadIncrement(end+1,:)] ...
                = hpReq(mission.segment(iSegment).altitude, ...
                gwIncrement, R, sigma, Vt, cd0, f, ...
                mission.segment(iSegment).speed, nRotors, isCoaxial,...
                isDuct, mission.segment(iSegment).roc, ...
                mission.segment(iSegment).isIdle,...
                mission.segment(iSegment).isIGE);
            % Calculate fuel 
            fuelIncrement(end+1,:) = ...
                sfc_baseline .* ...
                hpReqIncrement(end,:) .* lastIterTime;
            gwIncrement = gwIncrement - fuelIncrement(end,:);
            
        end
        mission.segment(iSegment).wFuel = sum(fuelIncrement);
        mission.segment(iSegment).hpReq = max(hpReqIncrement);
        mission.segment(iSegment).blLoad = max(blLoadIncrement);
       
    else
        % regular mission segment where time is less than 10 min
        [mission.segment(iSegment).hpReq,...
            mission.segment(iSegment).blLoad]...
            = hpReq(mission.segment(iSegment).altitude, ...
            mission.segment(iSegment).gw , R, sigma, Vt, cd0, f, ...
            mission.segment(iSegment).speed, nRotors, isCoaxial,isDuct,...
            mission.segment(iSegment).roc, ...
            mission.segment(iSegment).isIdle,...
            mission.segment(iSegment).isIGE);
        % Calculate fuel 
        mission.segment(iSegment).wFuel = ...
            sfc_baseline .* ...
            mission.segment(iSegment).hpReq .* ...
            mission.segment(iSegment).time;
        
    end
    
end

% Finds total fuel for each GW to calculate fuel fraction
wFuelTotal = zeros(1,length(GW));
for iSegment = 1:length(mission.segment)
    wFuelTotal = wFuelTotal + mission.segment(iSegment).wFuel;
end
fuelFractionReq = wFuelTotal ./ GW;

% Finds maximum continuous power required at any point in the mission for
% each GW
maxContPowerReq = zeros(1,length(GW));
for iSegment = 1:length(mission.segment)
    h = mission.segment(iSegment).altitude;
    [~, ~, ~, rho_h] = atmoscoesa(convlength(h,'ft','m'));
    rho_h = convdensity(rho_h, 'kg/m^3', 'slug/ft^3');
    mission.segment(iSegment).hpReq = (rho_sealevel ./ rho_h)...
        .* mission.segment(iSegment).hpReq; % @ sea-level power req
    mask =  mission.segment(iSegment).hpReq > maxContPowerReq;
    maxContPowerReq(mask) = mission.segment(iSegment).hpReq(mask);
end

% Finds max blade loading any point in the mission
maxBladeLoading = zeros(1,length(GW));
for iSegment = 1:length(mission.segment)
    mask = mission.segment(iSegment).blLoad > maxBladeLoading;
    maxBladeLoading(mask) = mission.segment(iSegment).blLoad(mask);
end