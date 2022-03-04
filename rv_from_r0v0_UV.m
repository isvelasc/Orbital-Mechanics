function [r,v] = rv_from_r0v0_UV(r0_vec, v0_vec, t, mu)

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function computes the state vector from an initial
    %       state vector and an elapsed time using the solution
    %       of the universal variable equation
    %
    %   Input:
    %       r0_vec - initial position vector (km)
    %       v0_vec - initial velocity vector (km/s)
    %       t      - elapsed time (s)
    %       mu     - gravitational parameter (km^3/s^2)
    %   Output:
    %       r      - final position vector (km)
    %       v      - final velocity vector (km/s)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    [chi, alpha]  = Kepler_UV(r0_vec, v0_vec, t, mu);
    [f,g]         = FG_UV(alpha, chi, norm(r0_vec), t, mu);
    r             = f*r0_vec + g*v0_vec;
    [f_dot,g_dot] = FG_dot_UV(alpha, chi, norm(r0_vec), norm(r), mu);
    v             = f_dot*r0_vec + g_dot*v0_vec;
    
end

