function TM = CTransform_Matrix(coes, ra_vec, ha_vec)

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function is the "Coordinate Transformation Matrix"
    %       needed to get from ECI to RSW
    %
    %   Input:
    %       coes - classical orbital elements {struct format}
    %              ω - argument of perigee (rads)
    %              θ - true anomaly (rads)
    %              Ω - right ascencion of ascending node (rads)
    %              i - inclination (rads)
    %
    %       ra_vec - target radius (row vector)
    %       ha_vec - target angular momentum (row vector)
    %
    %   Usage:
    %       'Nan' for coes if target vectors are used
    %       'NaN' for ra_vec and ha_vec if coes are used
    %
    %   Output:
    %       TM - transformation matrix
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if isnan(coes)
        
        i_hat = ra_vec/norm(ra_vec);
        k_hat = ha_vec/norm(ha_vec);
        j_hat = cross(k_hat, i_hat);
       
        TM = [i_hat;j_hat;k_hat];
    
    else
        u = coes.omega + coes.theta;
        r = coes.raan;
        i = coes.inc;

        R3_u = [cos(u), sin(u), 0
               -sin(u), cos(u), 0
                  0,      0,    1];
        R1_i = [1,    0,      0
                0,  cos(i), sin(i)
                0, -sin(i), cos(i)];
        R3_r = [cos(r), sin(r), 0
               -sin(r), cos(r), 0
                  0,      0,    1];

        TM = R3_u*R1_i*R3_r;
    end
    
end

