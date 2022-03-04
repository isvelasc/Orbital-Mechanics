function [TLE,tstep,sstep] = TLEpropagator(TLEfname,time,opts,ret)

    %{
        NOTE: This function does not take into account any petubational
        effects and therfore should only be used for either -1 or +1
        rotation at most. Assumes orbit is elliptical.
        
        This function rudimentarily propagates a TLE.
        Input:
            TLEfname - TLE file name
            time     - Time to be propagated decimal format range
                        should be -1*period<=time<=1*period (s)
            opts     - 'f' for file format
                       's' for string format
            ret      - return propagation values 0 for no, any other
                        integer value for yes
        Output:
            TLE      - new TLE in TLE format
            tstep    - time step vector from ode45
            sstate   - state step vector from ode45
    %}

    % check input parameters
    if time == 0
        fprintf('No propagation needed.\n');
        TLE = 'No propagation';
        return
    end

    if opts ~= 's' && opts ~= 'f'
        error('Unknown option %s',opts);
    end

    if nargin == 3
        retval = 0;
    elseif nargin == 4
        retval = ret;
    end


    tleinfo = readTLE(TLEfname);
    coes    = coesfromTLE(tleinfo);
    [r, v]  = coes2sv(coes.h, coes.ecc, coes.inc, coes.raan, coes.w, coes.TA, 398600);
    coes    = sv2coes(r,v);

    % check for correct interval
    if abs(time) > coes.T
        warning('ERROR: TLEpropagator does not account for pertubations. Predictions may outside of the scope of this function.')
    end
    
    % propagate
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    dt      = [0 time];
    state   = [r' v'];
    [t,s]   = ode45(@orbitprop,dt,state,options,398600);
    r       = s(end,1:3)';
    v       = s(end,4:6)';


    % grab old TLE paramters as strings
    TLEfid  = fopen(TLEfname,'rb');
    if TLEfid < 0
        error('Error. \nTLE file %s not found.',TLEfname)
    end
    name   = strtrim(fgetl(TLEfid));
    L1     = strtrim(fgetl(TLEfid));
    L2     = strtrim(fgetl(TLEfid));
    fclose(TLEfid);

    % set new TLE parameters
    % new epoch time
    L1(19:32) = num2str(str2double(L1(19:32)) + time/24/3600, '%.8f');
    % change element set number if possible
    if str2num(strtrim(L1(65:68))) < 999
        if time == 1 || time == -1
            esn = num2str(str2num(strtrim(L1(65:68))) + time);
            if length(esn) ~= 4
                L1(65:68) = padwhspc(esn,4 - length(esn));
            else
                L1(65:68) = esn;
            end
        end
    end
    % update coes
    coes      = sv2coes(r,v);
    inc       = num2str(rad2deg(coes.inc), '%.4f');
    L2(9:16)  = padwhspc(inc, 8 - length(inc));
    raan      = num2str(rad2deg(coes.raan), '%.4f');
    L2(18:25) = padwhspc(raan, 8 - length(raan));
    L2(27:33) = delimitleft(num2str(coes.ecc,'%.7f'));
    argp      = num2str(rad2deg(coes.omega), '%.4f');
    L2(35:42) = padwhspc(argp, 8 - length(argp));
    me        = num2str(rad2deg(coes.Me), '%.4f');
    L2(44:51) = padwhspc(me, 8 - length(me));
    n         = num2str(sqrt(398600/coes.a^3)*24*3600/(2*pi), '%.8f');
    L2(53:63) = n;

    % update revolutions
    if time == 1 || time == -1
        L2(64:68) = num2str(str2num(L2(64:68)) + time);
    end

    if retval
        tstep = t;
        sstep = s;
    else
        tstep = nan;
        sstep = nan;
    end

    data = strjoin({name,L1,L2},'\n');
    if opts == 's'
        TLE = data;
    elseif opts == 'f'
        newTLE = strcat(name,' Propagated TLE.txt');
        fid = fopen(newTLE, 'w');
        fwrite(fid, data);
        fclose(fid);
        TLE = newTLE;
    end

    function n = delimitleft(num)
        while contains(num,'.')
            num = num(2:end);
        end
        n = num;
    end

    function s = padwhspc(str,len)
        s = sprintf([blanks(len),'%s'],str);
    end
end