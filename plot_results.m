%% Plotting results

% Printing units
plotCounter = 1;
units = 'inches';
width = 6.25;
height = 5;


if isPlotting
    figure(1)
    plot(GW(fuelFractionAvail > 0), ...
        fuelFractionAvail(fuelFractionAvail > 0),'k--', ...
        'DisplayName', 'R_f avail.')
    hold on
    plot(GW, fuelFractionReq, 'k-','DisplayName' , 'R_f req.')
    plot(designGw, Rf, 'ks', ...
        'DisplayName',  sprintf('GW = %d lb\nR_f = %0.4f', designGw, Rf));
    title('Weight Balance Sizing for Baseline A-160 Hummingbird')
    xlabel('Gross Weight, GW (lb)')
    ylabel('Fuel Weight Fraction, R_f')
    legend('show', 'location', 'best')
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    hold off
    
    
    time_vec = [];
    Power_consmp = [];
    GrossWeight_Change = [];
    time_seg = [];
    for ii = 1:length(mission.segment)
        if ii == 1
            time_vec(end+1) = mission.segment(ii).time;
        else
            time_vec(end+1) = time_vec(end)+ mission.segment(ii).time;
        end
        Power_consmp(end+1) = mission.segment(ii).hpReq(baseMin);
        GrossWeight_Change(end+1) = mission.segment(ii).gw(baseMin);
        time_seg(end+1) = mission.segment(ii).time;
    end
    
    
    figure
    title('Optimum design Power consumption')
    subplot(1,2,1)
    bar(Power_consmp)
    title('Optimum Design Power consumption')
    xlabel('Segments')
    ylabel('Power [hp]')
    axis('xy')
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    segmts = 1:12;
    str = {'WarmUp';'Hover-taxi';'Take-off';'Cruise Flight';...
        'Hover Station 1';'Cruise Flight';'Hover Station 2';...
      'Cruise Flight';'Hover Station 3';'Landing';'Hover-Taxi';'Shutdown'};
    tab1 = table(segmts',str,time_seg','VariableNames',{'Segments',...
    'Station','Time'});
    
    %legend(str,'Location','best')
    colormap parula
    subplot(1,2,2)
    bar(GrossWeight_Change)
    title('Optimum Design Gross Weight Change')
    xlabel('Segments')
    ylabel('Gross Weight [lb]')
    ax = gca;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    str = {'1 WarmUp';'2 Hover-taxi';'3 Take-off';'4 Cruise Flight';...
        '5 Hover Station 1';'6 Cruise Flight';'7 Hover Station 2';...
        '8 Cruise Flight';'9 Hover Station 3';'10 Landing';...
        '11 Hover-Taxi';'12 Shutdown'};
    text(11,designGw-500,str)
   
end
