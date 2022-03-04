function [f_dot,g_dot] = FG_dot_UV(alpha, chi, r0, r, mu)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function calculates the time derivatives of
    %       the Lagrange coefficients f and g in terms of 
    %       change in the universal anomaly
    %
    %   Input:
    %       alpha - reciprocal of the semimajor axis (1/km)
    %       chi   - the universal anomaly after time t (km^0.5)
    %       r0    - the radial position (scalar) at initial time (km)
    %       r     - the radial position (scalar) after time t (km)
    %       mu    - gravitational parameter (km^3/s^2)
    %   Output:
    %       f_dot - time derivative of the Lagrange f coefficient (1/s)
    %       g_dot - time derivative of the Lagrange g coefficient (none)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    z = alpha*chi^2;
    f_dot = sqrt(mu)/r/r0*(z*stumpS(z) - 1)*chi;
    g_dot = 1 - chi^2/r*stumpC(z);
    
end

