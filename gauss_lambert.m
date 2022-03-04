function [v0,v1] = gauss_lambert(r0,r1,dt,tm,deg,mu)

    %{
        This function uses Gauss' trigonometric approach to solve 
            Lambert's problem.

        Input:
            r0 - first position vector (km)
            r1 - second position vector (km)
            dt - time of flight between the two vectors (s)
            tm - transfer method
                        for short way -> 1   Δθ < 180 deg
                        for long way -> -1   Δθ > 180 deg
            deg - degree of taylor series expansion
            mu  - gravitational parameter (km^3/s^2)

        Output:
            v0 - first velocity vector (km/s)
            v1 - second velocity vector (km/s)
    %}

    nr0   = norm(r0);
    nr1   = norm(r1);
    cdt   = dot(r0,r1)/(nr0*nr1); % cos(Δθ)
    dth   = acos(cdt); % Δθ
    chdt  = cos(dth/2); % cos(Δθ/2)
    sdt   = tm*sqrt(1 - cdt^2); % sin(Δθ)
    l     = (nr0 + nr1)/(4*sqrt(nr0*nr1)*chdt) - 0.5;
    m     = (mu*dt^2)/(2*sqrt(nr0*nr1)*chdt)^3;
    y0    = 1;
    ERR   = 1;
    TOL   = 1e-6;
    ts    = taylor_frac_expansion(deg);

    while ERR > TOL
        x1  = m/y0^2 - l;
        x2  = taylor_expansion(deg,x1,ts);
        y   = 1 + x2*(l + x1);
        ERR = abs(y - y0);
        y0  = y;
    end

    chde = 1 - 2*x1; % cos(ΔE/2)
    p    = (nr0*nr1*(1 - cdt))/(nr0 + nr1 - 2*sqrt(nr0*nr1)*chdt*chde); % parameter

    % find lagrange coefficients
    f    = 1 - nr1/p*(1 - cdt);
    g    = nr1*nr0*sdt/sqrt(mu*p);
    gdot = 1 - nr0/p*(1 - cdt);

    % solve for vector velocities
    v0   = (r1 - f*r0)/g;
    v1   = (gdot*r1 - r0)/g;

    function ts = taylor_frac_expansion(deg)
        ts = ones(deg+1,1);
        numer = 6;
        denom = 5;
        incr  = 2;

        % populate fraction expansion
        for step = 2:deg+1
            ts(step) = ts(step-1)*(numer/denom);
            numer = numer + incr;
            denom = denom + incr;
        end
    end

    function s = taylor_expansion(deg,x1,ts)        
        series = ones(deg+1,1);

        % apply power expansion
        for step = 2:deg+1
            series(step) = ts(step)*x1^(step-1);
        end

        s = 4/3*sum(series);
    end
end