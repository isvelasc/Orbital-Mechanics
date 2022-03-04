function [Az,el,rho,rhosez] = sv2razel(r,Rs,lst,gd)
    
    %{
        This function converts state vector in ECI to azimuth, elevation
        and slant range in SEZ.

        Input:
            r   - row positiion vector (km)
            Rs  - row site vector (km)
            lst - local sidereal time (deg)
            gd  - geodedict latitude (deg)
        Output:
            az  - azimuth (deg)
            el  - elevation (deg)
            rho - slant range (km)
    %}

    rhoeci  = r - Rs;
    tempvec = rotz(lst)\rhoeci';
    rhosez  = roty(90-gd)\tempvec;
    rho     = norm(rhosez);
    l       = rhosez(1)/rho;
    m       = rhosez(2)/rho;
    n       = rhosez(3)/rho;
    el      = asind(n);
    Az      = atan2d(m,-l);
    if Az < 0
        Az = Az + 360;
    end
    
end