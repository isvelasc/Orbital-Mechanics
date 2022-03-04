function Rs = RS_ECI(sigma,H,theta)
    
    %{
        This function calculates the site location in ECI coordinates
        
        Input:
            Re    - radius of Earth at equator (km)
            f     - oblateness factor
            sigma - geodetic latitude (rad)
            H     - elevation of the site (km)
            theta - local sidereal time (rad)
        Output:
            Rs    - site loction vector in ECI
    %}

    Re = 6378;
    f  = 0.003353;
    Z1 = cos(sigma)*(Re/(sqrt(1 - (2*f - f^2)*sin(sigma)^2)) + H);
    Z2 = (Re*(1 - f)^2)/(sqrt(1 - (2*f - f^2)*sin(sigma)^2)) + H;
    
    Rs = [Z1*cos(theta); Z1*sin(theta); Z2*sin(sigma)];

end