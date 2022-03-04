function [Az,el] = Az_el_from_RA_DEC(alpha,delta,lat,lst,option)

    %{
        This function calculates the azimuth and elevation from
            right ascension and declination

        Input:
            alpha  - Right Ascension
            delta  - Declination
            lat    - Latitude
            HA     - Hour Angle
            option - 'd' for degrees, 'r' for radians
        Output:
            Az     - Azimuth (degrees)
            el     - elevation (degrees)
    %}

    HA = lst - alpha;

    if option == 'd'
        el = asind(sind(delta)*sind(lat) + cosd(delta)*cosd(lat)*cosd(HA));
        if sind(HA) < 0
            Az = acosd((sind(delta) - sind(el)*sind(lat))/(cosd(el)*cosd(lat)));
        else
            Az = 360 - acosd((sind(delta) - sind(el)*sind(lat))/(cosd(el)*cosd(lat)));
        end
    elseif option == 'r'
        el = rad2deg(asin(sin(delta)*sin(lat) + cos(delta)*cos(lat)*cos(HA)));
        if sin(HA) < 0
            Az = rad2deg(acos((sin(delta) - sind(el)*sin(lat))/(cosd(el)*cos(lat))));
        else
            Az = 360 - rad2deg(acos((sin(delta) - sind(el)*sin(lat))/(cosd(el)*cos(lat))));
        end
    else
         error("option: 'd' or 'r' for degrees or radians.")
    end

    % correct direction
    if Az < 0
        Az = 360 + Az;
    end

end