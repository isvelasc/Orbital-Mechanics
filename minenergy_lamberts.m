function [dtm,dtp,v0] = minenergy_lamberts(r0,r1,tm,mu)

    %{
        This function computes the minimum time required
            between to two points in an orbit.
            It assumes that the orbit occurs within 1 period.
       
        Input:
            r0 - first position vector (km)
            r1 - second position vector (km)
            tm - transfer method
                        for short way -> 1   Δθ < 180 deg
                        for long way -> -1   Δθ > 180 deg
            mu - gravitational parameter (km^3/s^2)

        Output:
            dtm - time transfer for minimum energy (s)
            dtp - time for parabolic transfer for minimum energy (s)
            v0  - initial velocity vector (km/s)
    %}

    nr0   = norm(r0);
    nr1   = norm(r1);
    cdt   = dot(r0,r1)/(nr0*nr1); % cos(Δθ)
    sdt   = tm*sqrt(1 - cdt^2); % sin(Δθ)
    c     = sqrt(nr0^2 + nr1^2 - 2*nr0*nr1*cdt); % chord length
    s     = (nr0 + nr1 + c)/2; % semi-parameter
    am    = s/2; % semi major axis for minimum energy
    pm    = nr0*nr1/c*(1 - cdt); % min semi-parameter
    alpha = pi;
    betam = asin(sqrt((s-c)/s))*2;

    if tm == 1
        dtm = (am^1.5*(alpha - (betam - sin(betam))))/sqrt(mu);
    elseif tm == -1
        dtm = (am^1.5*(alpha + (betam - sin(betam))))/sqrt(mu);
    else
        error('tm not correct value.')
    end

    dtp = ((sqrt(2)/3)*(s^1.5 - (s - c)^1.5))/sqrt(mu);
    v0  = sqrt(mu*pm)/(nr0*nr1*sdt)*(r1 - r0*(1 - nr1/pm*(1 - cdt)));
end