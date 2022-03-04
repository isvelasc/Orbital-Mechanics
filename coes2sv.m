function [r, v] = coes2sv(h, ecc, inc, raan, omega, theta, mu)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    %   This function computes the state vector (r, v)          
    %    from the classical orbital elements (coes).            
    %   
    %   Input:
    %       mu    - gravitational parameter (km^3/s^2)              
    %       h     - angular momentum (km^2/s)                       
    %       e     - eccentricity                                    
    %       inc   - inclination of the orbit (rad)                  
    %       raan  - right ascension of the ascending node (rad)     
    %       omega - argument of perigee (rad)                       
    %       theta - true anomaly (rad)
    %   Output:
    %       r     - position vector (km)
    %       v     - velocity vector (km/s)
    %                                  
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    

    % Rotation Matrices
    R3     = @(n) [cos(n) sin(n) 0;
                  -sin(n) cos(n) 0;
                     0      0    1];
               
    R1     = @(n) [1    0      0;
                   0  cos(n) sin(n);
                   0 -sin(n) cos(n)];
               
    % Initial Perifocal Frame
    radius_perifocal   = (h^2 / mu) * (1 / (1 + ecc * cos(theta))) .* [cos(theta); sin(theta); 0];
    velocity_perifocal = (mu / h) .* [-sin(theta); ecc + cos(theta); 0];
    
    % Matrix transformation from perifocal to ECI
    Q_p_to_ECI = transpose(R3(omega) * R1(inc) * R3(raan));
    
    % Set to column vectors
    radius   = Q_p_to_ECI * radius_perifocal;
    velocity = Q_p_to_ECI * velocity_perifocal;
    
    % Set to row vectors
    r = radius;
    v = velocity;
end

