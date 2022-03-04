function st = LST(DT,L)
    %{
        This function calculates the local sidereal time
            
        Input:
            DT - Datetime vector in the form of [Y M D UT]
                    Y  - Year
                    M  - Month
                    D  - Day
                    UT - Universal Time H M S (Hours Minutes Seconds)
                         or decimal format

                Example: [2016 07 29 16 05 24]
                Example: [2016 07 29 16.09]

            L  - Longitude vector in the form of [D M S] or decimal format
                    D - Degrees     -359 <= D <= 359  -> negative connotes
                    west, positive connotes east
                    M - Minutes     0 <= M <= 59
                    S - Seconds     0 <= S <= 59

                Example: [144 58 19]
                Example: 144.9719

        Output:
            st - local sidereal time in degrees
    %}
    if isempty(DT) | isempty(L)
        error('Check for empty input arguments.')
    end
    if length(DT) < 4
        error('Datetime does not have enough arguments.')
    end

    if DT(2) < 1 || DT(2) > 12
        error('Datetime Month beyond scope.')
    end
    if DT(3) < 1 || DT(3) > 31 || (DT(1) == 2 && DT(3) > 29)
        error('Datetime Day beyond scope.')
    end

    Y  = DT(1);
    M  = DT(2);
    D  = DT(3);
    UT = DT(4:end);

    if UT(1) > 23 || UT(1) < 0
        error('UT Hours beyond scope.')
    end
    if length(UT) > 1 && (UT(2) > 59 || UT(2) < 0)
        error('UT Minutes beyond scope.')
    end
    if length(UT) > 2 && (UT(3) > 59.999999 || UT(3) < 0)
        error('UT Seconds beyond scope.')
    end
    
    if length(UT) == 2
        UT = UT(1) + UT(2)/60;
    elseif length(UT) == 3
        UT = UT(1) + UT(2)/60 + UT(3)/3600;
    end

    if L(1) >= 360
        error('L Degrees beyond scope.')
    end
    if length(L) > 1 && L(2) > 59
        error('L Minutes beyond scope.')
    end
    if length(L) > 2 && L(3) > 59.999999
        error('L Seconds beyond scope.')
    end

    if L(1) < 0
        if length(L) == 2
            el = 360 + L(1) - L(2)/60;
        elseif length(L) == 3
            el = 360 + L(1) - L(2)/60 - L(3)/3600;
        else
            el = 360 + L(1);
        end
    else
        if length(L) == 2
            el = L(1) + L(2)/60;
        elseif length(L) == 3
            el = L(1) + L(2)/60 + L(3)/3600;
        else
            el = L(1);
        end
    end

    
    J0 = 367*Y - floor((7*(Y+floor((M+9)/12)))/4) + floor((275*M)/9) + D + 1721013.5;
    T0 = (J0 - 2451545)/36525;

    Theta_G0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0^2 - 2.58*(10^-8)*T0^3;
    Theta_G = Theta_G0 + 360.98564724*(UT/24);

    Theta = Theta_G + el;
    st = mod(Theta,360);
end