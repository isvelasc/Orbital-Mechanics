function [delta_r, delta_v] = CW_Equation_Format(delta_r0, delta_v0, n, t)

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function calculates the simplified version of
    %       Equations of Motion (EOM) using the 
    %       Clohessy-Wiltshire equations in equation format
    %       using Laplace Transforms, a.k.a Dr. A's format.
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

    delta_x0     = delta_r0(1);
    delta_y0     = delta_r0(2);
    delta_z0     = delta_r0(3);
    delta_x0_dot = delta_v0(1);
    delta_y0_dot = delta_v0(2);
    delta_z0_dot = delta_v0(3);
    
    delta_x     = 4*delta_x0 + 2*delta_y0_dot/n + delta_x0_dot/n*sin(n*t)...
                    - (3*delta_x0 + 2*delta_y0_dot/n)*cos(n*t);
    delta_y     = 2*delta_x0_dot/n*cos(n*t) + (6*delta_x0 + 4*delta_y0_dot/n)*sin(n*t)...
                    - (6*n*delta_x0 + 3*delta_y0_dot)*t - 2*delta_x0_dot/n + delta_y0;
    delta_z     = delta_z0*cos(n*t) + delta_z0_dot/n*sin(n*t);
    delta_x_dot = delta_x0_dot*cos(n*t) + (3*n*delta_x0 + 2*delta_y0_dot)*sin(n*t);
    delta_y_dot = (6*n*delta_x0 + 4*delta_y0_dot)*cos(n*t) - 2*delta_x0_dot*sin(n*t)...
                    - (6*n*delta_x0 + 3*delta_y0_dot);
    delta_z_dot = -delta_z0*n*sin(n*t) + delta_z0_dot*cos(n*t);
    
    delta_r = [delta_x; delta_y; delta_z];
    delta_v = [delta_x_dot; delta_y_dot; delta_z_dot];
end

