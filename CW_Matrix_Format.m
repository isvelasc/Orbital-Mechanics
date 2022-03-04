function [delta_r, delta_v] = CW_Matrix_Format(delta_r0, delta_v0, n, t)

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function calculates the simplified version of
    %       Equations of Motion (EOM) using the 
    %       Clohessy-Wiltshire equations in matrix format
    %       following Curtis.
    %
    %   Input:
    %       delta_r0 - initial position vector (m, column vector)
    %       delta_v0 - initial velocity vector (m/s, column vector)
    %       n        - mean motion of target (rad/s)
    %       t        - time since last position (s)
    %
    %   Output:
    %       delta_r - final position vector (m)
    %       delta_v - final velocity vector (m/s)
    %
    %
    %   NOTE: This function only works if rho/R << 1
    %
    %       i.e.     the relative position magnitude
    %                -------------------------------  << 1
    %                   chaser position magnitude
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    phi_rr = [4 - 3*cos(n*t),        0,    0
              6 - (sin(n*t) - n*t),  1,    0
                     0,              0, cos(n*t)];
    
    phi_rv = [   1/n*sin(n*t),       2/n*(1 - cos(n*t)),          0
              2/n*(cos(n*t) - 1), 1/n*(4*sin(n*t) - 3*n*t),       0
                     0,                      0,             1/n*sin(n*t)];
    
    phi_vr = [   3*n*sin(n*t),    0,     0
              6*n*(cos(n*t) - 1), 0,     0
                       0,         0, -n*sin(n*t)];
    
    phi_vv = [  cos(n*t),    2*sin(n*t),      0
              -2*sin(n*t), 4*cos(n*t) - 3,    0
                   0,            0,        cos(n*t)];
    
    delta_r = phi_rr*delta_r0 + phi_rv*delta_v0;
    delta_v = phi_vr*delta_r0 + phi_vv*delta_v0;

end

