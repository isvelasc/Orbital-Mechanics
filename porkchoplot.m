function [dVshort,dVlong,bestdv,C3short,C3long,Vinfshort,Vinflong,TOF] = porkchoplot(departloc, arriveloc, departtime, arrivetime, opts)

    %{
        This function calculates and plots pork chop plots for orbit
        transfer between any two planets.

        Input:
            departloc   - string denoting departure location
            arriveloc   - string denoting arrival location
                          Both must be from one of the following:
                            Mercury
                            Venus
                            Earth
                            Mars
                            Jupiter
                            Saturn
                            Uranus
                            Neptune
                            Pluto
            departtime - 2D datetime vector with date window
            arrivetime - 2D datetime vector with date window
            opts       - optional struct input with the following optional
                         parameters:
                            lowerC3    - lower bound for C3
                            upperC3    - upper bound for C3
                            lowerVinf  - lower bound for Vinf
                            upperVinf  - upper bound for Vinf
                            C3levels   - specifies C3 contour lines to display
                            Vinflevels - specifies Vinf contour lines to display
                            TOFlevels  - specifies time of flight contour lines to display
                            optimalt   - plots the best dV for departure
                                         and arrival

        Example:
            opts.lowerC3    = 4;
            opts.upperC3    = 110;
            opts.lowerVinf  = 2;
            opts.upperVinf  = 10;
            opts.C3levels   = [16,17,18,22,28,36,45,56,70,100];
            opts.Vinflevels = [2.5,3,4,5,7,10];
            opts.TOFlevels  = 50:50:500;
            opts.optimalt   = 1;

            porkchoplot('Earth', 'Mars', [datetime('06-Jun-2005 12:00:00'),...
                                          datetime('07-Nov-2005 12:00:00')],...
                                         [datetime('01-Dec-2005 12:00:00'),...
                                          datetime('24-Feb-2007 12:00:00')], opts)
    %}

    planets = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'}; 

    % check inputs
    
    options.lowerC3    = 10;
    options.upperC3    = 120;
    options.lowerVinf  = 0;
    options.upperVinf  = 15;
    options.C3levels   = 10:5:120;
    options.Vinflevels = 0:15;
    options.TOFlevels  = 50:50:500;
    options.optimalt   = 0;
    if isstruct(opts)
        if isfield(opts,'lowerC3') options.lowerC3 = opts.lowerC3; end
        if isfield(opts,'upperC3') options.upperC3 = opts.upperC3; end
        if isfield(opts,'lowerVinf') options.lowerVinf = opts.lowerVinf; end
        if isfield(opts,'upperVinf') options.upperVinf = opts.upperVinf; end
        if isfield(opts,'C3levels')   && isvector(opts.C3levels)   options.C3levels   = opts.C3levels;   end
        if isfield(opts,'Vinflevels') && isvector(opts.Vinflevels) options.Vinflevels = opts.Vinflevels; end
        if isfield(opts,'TOFlevels')  && isvector(opts.TOFlevels)  options.TOFlevels  = opts.TOFlevels;  end
        if isfield(opts,'optimalt') options.optimalt = opts.optimalt; end
    else
        warning('opts is not a struct. opts will be ignored.')
    end
    if ~any(strcmp(departloc,planets))
        error('Departing location is beyond the scope of this function')
    end
    if ~any(strcmp(arriveloc,planets))
        error('Arrival location is beyond the scope of this function')
    end
    if arrivetime(1) > arrivetime(2)
        error('Arrival time window must be incremental.')
    end
    if departtime(1) > departtime(2)
        error('Departure time window must be incremental.')
    end
    if strcmp(departloc,arriveloc)
        fprintf('Departure and arrival coincide\n');
        return
    end

    % acquire preliminary data
    departID  = find(strcmp(planets,departloc));
    arriveID  = find(strcmp(planets,arriveloc));
    departlen = days(departtime(2) - departtime(1));
    arrivelen = days(arrivetime(2) - arrivetime(1));
    departdt  = linspace(departtime(1), departtime(2), departlen);
    arrivedt  = linspace(arrivetime(1), arrivetime(2), arrivelen);
    departdv  = datevec(departdt);
    arrivedv  = datevec(arrivedt);


    dVshort   = zeros(arrivelen,departlen);
    dVlong    = zeros(arrivelen,departlen);
    Vinfshort = zeros(arrivelen,departlen);
    C3short   = zeros(arrivelen,departlen);
    Vinflong  = zeros(arrivelen,departlen);
    C3long    = zeros(arrivelen,departlen);
    TOF       = zeros(arrivelen,departlen);

    % find every combination
    for adt = 1:arrivelen
        for ddt = 1:departlen
            % planet ephemeris
            [~, depart_planet_r, depart_planet_v, ~] = AERO557planetcoe_and_sv(departID, departdv(ddt,1), departdv(ddt,2), departdv(ddt,3), departdv(ddt,4), departdv(ddt,5), departdv(ddt,6));
            [~, arrive_planet_r, arrive_planet_v, ~] = AERO557planetcoe_and_sv(arriveID, arrivedv(adt,1), arrivedv(adt,2), arrivedv(adt,3), arrivedv(adt,4), arrivedv(adt,5), arrivedv(adt,6));
            
            % sc velocities
            time = days(seconds(etime(arrivedv(adt,:),departdv(ddt,:))));
            [depart_sc_short_v, arrive_sc_short_v,~,~] = lambert(depart_planet_r', arrive_planet_r', time, 0, 1.327124e11);
            [depart_sc_long_v, arrive_sc_long_v,~,~]   = lambert(depart_planet_r', arrive_planet_r', -time, 0, 1.327124e11);

            C3s   = norm(depart_sc_short_v' - depart_planet_v)^2;
            C3l   = norm(depart_sc_long_v' - depart_planet_v)^2;
            Vinfs = norm(arrive_sc_short_v' - arrive_planet_v);
            Vinfl = norm(arrive_sc_long_v' - arrive_planet_v);

            % check data at cutoff 
            if C3s < options.lowerC3 C3short(adt,ddt) = options.lowerC3;
            elseif C3s > options.upperC3 C3short(adt,ddt) = options.upperC3;
            else C3short(adt,ddt) = C3s;
            end

            if C3l < options.lowerC3 C3long(adt,ddt) = options.lowerC3;
            elseif C3l > options.upperC3 C3long(adt,ddt) = options.upperC3;
            else C3long(adt,ddt) = C3l;
            end

            if Vinfs < options.lowerVinf Vinfshort(adt,ddt) = options.lowerVinf;
            elseif Vinfs > options.upperC3 Vinfshort(adt,ddt) = options.upperVinf;
            else Vinfshort(adt,ddt) = Vinfs;
            end

            if Vinfl < options.lowerVinf Vinflong(adt,ddt) = options.lowerVinf;
            elseif Vinfl > options.upperC3 Vinflong(adt,ddt) = options.upperVinf;
            else Vinflong(adt,ddt) = Vinfl;
            end            

            dVshort(adt,ddt) = Vinfshort(adt,ddt) + sqrt(C3short(adt,ddt));
            dVlong(adt,ddt)  = Vinflong(adt,ddt) + sqrt(C3long(adt,ddt));
            TOF(adt,ddt)     = time;
        end
    end

    % constrain any errors from plotting
    warnID     = 'MATLAB:contour:ConstantData';
    warnStruct = warning('off',warnID);

    % get array for ploting optimal time given best delta v
    bsv = min(min(dVshort));
    blv = min(min(dVlong));
    if bsv < blv bestdv = bsv; else bestdv = blv; end
    if bestdv == bsv
        [r,c] = find(dVshort==bestdv);
    else
        [r,c] = find(dVlong==bestdv);
    end
    locbdv = strjoin({strcat('(',num2str(r),','), strcat(num2str(c),')')});


    figure
    hold on
    contour(C3long,options.C3levels,'r','ShowText','on')
    contour(C3short,options.C3levels,'r','ShowText','on')
    contour(Vinflong,options.Vinflevels,'b','ShowText','on')
    contour(Vinfshort,options.Vinflevels,'b','ShowText','on')
    contour(TOF,options.TOFlevels,'k','ShowText','on')
    if options.optimalt
        scatter(r,c,'MarkerFaceColor','#7E2F8E','MarkerEdgeColor','#7E2F8E')
        line([0,r],[c,c],'Color','#7E2F8E','LineStyle','--')
        line([r,r],[0,c],'Color','#7E2F8E','LineStyle','--')
        text(r,c,locbdv,'FontWeight','bold','Color','#7E2F8E',...
            'FontSize',10,'Position',[r c+15],'HorizontalAlignment','center')
    end
    labelx = strjoin(['Days Past', string(month(departtime(1),'shortname')),...
             strcat(num2str(day(departtime(1))), ','),...
             num2str(year(departtime(1),'gregorian'))]);
    labely = strjoin(['Days Past', string(month(arrivetime(1),'shortname')),...
             strcat(num2str(day(arrivetime(1))), ','),...
             num2str(year(arrivetime(1),'gregorian'))]);
    xlabel(labelx)
    ylabel(labely)
    if options.optimalt
        legend('C3 $\left(\frac{km^2}{s^2}\right)$','',...
            arrivelabel(arriveID),'','Time of Flight (days)',...
            'Best $\Delta$V','Interpreter','latex','Orientation','horizontal',...
            'FontSize',7,'Location','northoutside')
    else
    legend('C3 $\left(\frac{km^2}{s^2}\right)$','',...
            arrivelabel(arriveID),'','Time of Flight (days)',...
            'Interpreter','latex','NumColumns',12,'FontSize',7,'Location','northoutside')
    end
    hold off
    

    function label = arrivelabel(ID)
        labels = {'V$_\infty$ at Mercury $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Venus $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Earth $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Mars $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Jupiter $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Saturn $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Uranus $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Neptune $\left(\frac{km}{s}\right)$',...
                  'V$_\infty$ at Pluto $\left(\frac{km}{s}\right)$'};
        label = labels{ID};
    end

end