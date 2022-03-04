function [f,g] = FG_UV(alpha, chi, r0, t, mu)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function calculates the Lagrange coefficients
    %       f and g in terms of change in the universal anomaly
    %
    %   Input:
    %       alpha - reciprocal of the semimajor axis (1/km)
    %       chi   - the universal anomaly after time t (km^0.5)
    %       r0    - the radial position (scalar) at initial time (km)
    %       t     - the time elapsed since inital state vector (s)
    %       mu    - gravitational parameter (km^3/s^2)
    %   Output:
    %       f     - the Lagrange f coefficient (none)
    %       g     - the Lagrange g coefficient (s)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    z = alpha*chi^2;
    f = 1 - chi^2/r0*stumpC(z);
    g = t - 1/sqrt(mu)*chi^3*stumpS(z);
    
end

