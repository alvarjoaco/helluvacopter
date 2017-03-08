%% Mission Defination

% Mission Segments

% T1 - Warmup
mission.segment(1).altitude = baseAltitude;
mission.segment(1).isIdle = true;
mission.segment(1).time = warmUpTime;
mission.segment(2).speed = 0;

% T2 - Hover-taxi
mission.segment(2).speed = 0;
mission.segment(2).isIGE = true;
mission.segment(2).time = hoverTaxi;

% T3 - Take-off
mission.segment(3).altitudeChange = altitudeToHover - baseAltitude;
mission.segment(3).speed = convvel(0, 'kts', 'ft/s');
mission.segment(3).roc = vrocForClimb;

% T4 - Cruise flight
mission.segment(4).distance = convlength(segmentLength, 'naut mi', 'ft');
mission.segment(4).speed = convvel(cruiseSpeed, 'kts', 'ft/s');

% T5 - Hover - HOGE - hover station 1
mission.segment(5).speed = convvel( hoverSpeed, 'kts', 'ft/s');
mission.segment(5).time = hoverTimeAtStationA;

% T6 - Cruise flight
mission.segment(6).distance = convlength(segmentLength, 'naut mi', 'ft');
mission.segment(6).speed = convvel(cruiseSpeed, 'kts', 'ft/s');

% T7 - Hover - HOGE - hover station 2
mission.segment(7).speed = convvel( hoverSpeed, 'kts', 'ft/s');
mission.segment(7).time = hoverTimeAtStationB;

% T8 - Cruise flight
mission.segment(8).distance = convlength(segmentLength, 'naut mi', 'ft');
mission.segment(8).speed = convvel(cruiseSpeed, 'kts', 'ft/s');

% T9 - Hover - HOGE - hover station 3
mission.segment(9).speed = convvel( hoverSpeed, 'kts', 'ft/s');
mission.segment(9).time = totalHoverTime -hoverTimeAtStationA...
                          - hoverTimeAtStationB;

% T10 - Landing to base
mission.segment(10).altitudeChange = -altitudeToHover;
mission.segment(10).speed = convvel(0, 'kts', 'ft/s');
mission.segment(10).roc = -vrocForClimb;

% T11 - Hover-taxi
mission.segment(11).speed = 0;
mission.segment(11).isIGE = true;
mission.segment(11).time = hoverTaxi;

% T12 - Shutdown
mission.segment(12).speed = 0;
mission.segment(12).isIdle = true;
mission.segment(12).time = shutDownTime;

