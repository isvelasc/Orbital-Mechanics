function coes = coesfromTLE(TLEinfo)
    
    %{
        This function returns coes from TLE struct. Assumes 0<e<1 orbit.

        Input:
            TLEinfo:
                    name - Satellite Name
                  satnum - Satellite catalog number
                   class - Classification (U: unclassified, C: classified, S: secret)
              designator - International Designator (last two digits of launch year, launch number of the year, piece of the launch)
                   epoch - format [yy mm dd H M S]
                  motion - [First derivative of mean motion (the ballistic coefficient), Second derivative of mean motion]
                    drag - B*, the drag term, or radiation pressure coefficient
               ephemeris - Ephemeris type
                 element - Element set number
                checksum - Checksum (modulo 10)
                     inc - Inclination (deg)
                    raan - Right Ascension of Ascendion Node (deg)
                     ecc - Eccentricity
                    argp - Argument of Perigee (deg)
                      ma - Mean Anomaly (deg)
                       n - Mean Motion (rev/day)
                     rev - Revolution number at epoch

        Output:
            coes:
                     inc - Inclination (rad)
                    raan - Right Ascension of Ascendion Node (rad)
                     ecc - Eccentricity
                       w - Argument of Perigee (rad)
                      Me - Mean Anomaly (rad)
                       N - Mean Motion (rad/s)
                       a - Semi-major Axis (km)
                       h - Angular Momentum km^2/s
                      TA - True Anomaly (rad)
    %}

    mu        = 398600; % km^3/s^2
    coes.inc  = deg2rad(TLEinfo.inc);
    coes.raan = deg2rad(TLEinfo.raan);
    coes.ecc  = TLEinfo.ecc;
    coes.w    = deg2rad(TLEinfo.argp);
    coes.Me   = deg2rad(TLEinfo.me);
    coes.N    = TLEinfo.n*2*pi/24/3600;
    coes.a    = (mu/coes.N^2)^(1/3);
    coes.h    = sqrt(coes.a*mu*(1 - coes.ecc^2));
    coes.TA   = KeplersEq(coes.Me, coes.ecc);

end