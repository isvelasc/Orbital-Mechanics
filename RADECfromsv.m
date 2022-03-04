function [alpha,delta] = RADECfromsv(r,opts)

    %{
        This function calculates right ascension and declination from the
        position vector.
        Input:
            r - position vector (km)
            opts - 'd' for degrees, 'r' for radians
        Output:
            alpha - right ascension
            delta - declination
    %}


    rnorm = r/norm(r);
    l     = rnorm(1);
    m     = rnorm(2);
    n     = rnorm(3);

    if opts == 'r'
        delta = asin(n);
        if m > 0
            alpha = acos(l/cos(delta));
        else
            alpha = deg2rad(360 - acosd(l/cos(delta)));
        end
    elseif opts == 'd'
        delta = asind(n);
        if m > 0
            alpha = acosd(l/cosd(delta));
        else
            alpha = 360 - acosd(l/cos(delta));
        end
    else
        error('Unknown option.')
    end
end