function TLEinfo = readTLE(TLEfname)

    %{

        This function parses a Two Line Element data format encoding
            a list of orbital elements of an Earth-orbiting object 
            for a given point in time.

        Input:
            TLEfname: TLE filename

        Output:
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
                     ecc - Eccentricity (deg)
                    argp - Argument of Perigee (deg)
                      me - Mean Anomaly (deg)
                       n - Mean Motion (rev/day)
                     rev - Revolution number at epoch
    %}
    
    TLEfile = fopen(TLEfname,'rb');
    if TLEfile < 0
        error('Error. \nTLE file %s not found.',TLEfname)
    end

    NOS = strtrim(fgetl(TLEfile));
    L1  = strtrim(fgetl(TLEfile));
    L2  = strtrim(fgetl(TLEfile));
    
    % check for correct formating
    if length(NOS) > 24
        error('Name greater than 24 characters.')
    elseif length(L1) ~= 69
        error('Line 1 format not according to specifications. Length must be 69 characters long.')
    elseif length(L2) ~= 69
        error('Line 2 format not according to specifications. Length must be 69 characters long.')
    end

    if L1(1) ~= '1'
        error('Line 1 does not contain correct line number.')
    elseif L2(1) ~= '2'
        error('Line 2 does not contain correct line number.')
    end

    hasspace = [L1(2) L1(9) L1(18) L1(33) L1(44) L1(53) L1(62) L1(64)...
                L2(2) L2(8) L2(17) L2(26) L2(34) L2(43) L2(52)];
    if ~isspace(hasspace)
        error('Correct spacing format not adhered to.')
    end

    if L1(3:7) ~= L2(3:7)
        error('Satellite Catalog Number differs.')
    end

    TLEinfo.name           = NOS;
    TLEinfo.satnum         = str2num(L1(3:7));
    TLEinfo.class          = L1(8);
    designator.year        = str2num(L1(10:11));
    designator.launch      = str2num(L1(12:14));
    designator.piece       = strtrim(L1(15:17));
    TLEinfo.designator     = designator;
    year                   = fix(str2num(L1(19:32))/1000);
    TLEinfo.epoch          = datevec(datenum(year,0,str2double(L1(19:32)) - year*1000)); 
    motion.dt1             = str2double(L1(34:43));
    if contains(L1(45:52),'-')
        motion.dt2  = str2double(strtrim(L1(45:50)))*10^-str2num(L1(52));
    elseif contains(L1(45:52),'+')
        motion.dt2  = str2double(strtrim(L1(45:50)))*10^str2num(L1(52));
    else
        error('Second derivative of mean motion symbol must be (+) or (-).')
    end
    TLEinfo.motion    = motion;
    if contains(L1(55:61),'-')
        TLEinfo.drag  = str2double(strtrim(L1(54:59)))*10^-str2num(L1(61));
    elseif contains(L1(55:61),'+')
        TLEinfo.drag  = str2double(strtrim(L1(54:59)))*10^str2num(L1(61));
    else
        error('Radiation pressure coefficient symbol must be (+) or (-).')
    end

    TLEinfo.ephemeris      = str2num(L1(63));
    TLEinfo.element        = str2num(L1(65:68));
    TLEinfo.checksum       = str2num(L1(69));
    TLEinfo.inc            = str2double(strtrim(L2(9:16)));
    TLEinfo.raan           = str2double(strtrim(L2(18:25)));
    if L2(27) == '0'
        TLEinfo.ecc        = str2double(strcat('.',L2(27:33)));
    else
        TLEinfo.ecc        = str2double(L2(27:33));
    end        
    TLEinfo.argp           = str2double(strtrim(L2(35:42)));
    TLEinfo.me             = str2double(strtrim(L2(44:51)));
    TLEinfo.n              = str2double(strtrim(L2(53:63)));
    TLEinfo.re             = str2num(L2(64:68));

end