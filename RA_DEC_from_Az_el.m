function [alpha, delta] = RA_DEC_from_Az_el(Lat, Az, el, lst, option)
    %{
        This function calculates the declination and right ascension
            from azimuth and elevation
        
        Input:
            Lat    - Latitude
            Az     - Azimuth
            el     - Elevation
            lst    - Local Sidereal Time -> ALWAYS in degrees
            option - 'd' for degrees 'r' for radians
        Output:
            alpha  - right ascension (degrees)
            delta  - declination (degrees)
    %}
    
    if option == 'd'
        delta = asind(cosd(Lat)*cosd(Az)*cosd(el) + sind(Lat)*sind(el));
        if Az >= 180
            HA = acosd((cosd(Lat)*sind(el) - sind(Lat)*cosd(Az)*cosd(el))/cosd(delta));
        else
            HA = 360 - acosd((cosd(Lat)*sind(el) - sind(Lat)*cosd(Az)*cosd(el))/cosd(delta));
        end
        alpha = lst - HA;
    elseif option == 'r'
        delta = rad2deg(asin(cos(Lat)*cos(Az)*cos(el) + sin(Lat)*sin(el)));
        if Az >= pi
            HA = rad2deg(acos((cos(Lat)*sin(el) - sin(Lat)*cos(Az)*cos(el))/cosd(delta)));
        else
            HA = 360 - rad2deg(acos((cos(Lat)*sin(el) - sin(Lat)*cos(Az)*cos(el))/cosd(delta)));
        end
        alpha = lst - HA;
    else
        error("option: 'd' or 'r' for degrees or radians.")
    end

    % correct direction
    if alpha < 0
        alpha = 360 + alpha;
    end

end