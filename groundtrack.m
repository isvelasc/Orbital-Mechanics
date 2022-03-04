function out = groundtrack(data,time,opt)
    
    %{
        This function calculates and plots ground tracks over a given
        period.
        NOTE: The function currently does not account for any orbit
        deterioration or pertubations. 

        Inputs:
            opt  - a string describing one of three data sets the
                   function takes in.
                  'state' - state element
                  'tle'   - two line element
                  'coes'  - classical orbital elements 
            time - time to plot ground track for the given dataset
            data - 
                   state - Mx6 row verctor with position velocity components
                        OR
                   coes  - Mx6 table containing the following parameters:
                                h    - angular momentum (km^2/s)
                                ecc  - eccentricity
                                inc  - inclination (rad)
                                raan - right ascension of ascending node (rad)
                                argp - argumnent of perigee (rad)
                                ta   - true anomaly (rad)
                        OR
                   tle  - Mx1 string vector containing TLE file names
                               each file must be formatted to two line elemnet parameters, 
                               more info can be viewed here:
                               https://en.wikipedia.org/wiki/Two-line_element_set
                    
            Example:
            1 entry: 
                     data     = [−3048.0182, 0725.0827, 2759.1729, −9.2437, 1.5842, 7.6957]
                     time     = 5.3497e5
                     opt      = 'state'
                 OR
                        data.h = 3.0279e3
                      data.ecc = 0.9987
                      data.inc = 0.9035
                     data.raan = 0.5332
                     data.argp = 5.3756
                       data.ta = 3.0494
                     time      = 5.3497e5
                     opt       = 'coes'
                 OR
                     data = ['ISS TLE.txt']
                     time = 5.3497e5
                     opt  = 'tle'

            Multiple Entries
                    data = [-1984 -5348    3471   10.36  -5.763   -2.961;
                            48200 -2658   -24660   5.590  1.078   -3.484;
                            23047 -6972.4 -9219.6  6.6563 0.88638 -3.9680]
                 OR
                        data.h = 8.1590e4 8.1573e4 8.1589e4
                      data.ecc = 0.5012   0.4999   1.5012
                      data.inc = 0.6108   0.6109   0.6108
                     data.raan = 2.2691   2.2690   2.2691
                     data.argp = 2.0071   2.0071   2.0071
                       data.ta = 6.2830   2.0503   1.8077
                OR
                     data  = ["TLE 1.txt";
                              "TKE 2.txt";
                              "TLE 3.txt"]
    %}

    % check inputs
    if strcmp(opt,'state')
        if ~isa(data,'double') || length(data) ~= 6
            error('State vector contains either incorrect parameters or is incorrect length.')
        end
    elseif strcmp(opt,'coes')
        if ~istable(data) || ~ismember(dataset.Properties.VariableNames,{'h','ecc','inc','raan','argp','ta'})
            error('Coes is not a table or is missing one or more elements.')
        end
    elseif strcmp(opt,'tle')
        if ~isfile(data)
            error('TLE filename not a string.')
        end
    else
        error("Unkown option '%s'", opt)
    end

    if time == 0
        fprintf('TBD\n');
        return
    end

    


    out = nan;
   
end