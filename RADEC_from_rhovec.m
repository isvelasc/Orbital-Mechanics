function [alpha, delta] = RADEC_from_rhovec(rhovec)

    %{
        This function computes right ascension and declination from slant
        range vector. 
        Input:
            rhovec - slant range vector (km)
        Output:
            alpha - right ascension (rad)
            delta - declination (rad)
    %}

    unitrho = rhovec/norm(rhovec);
    I       = unitrho(1);
    J       = unitrho(2);
    K       = unitrho(3);
    delta   = asin(K);
    alpha   = acos(I/cos(delta));

end