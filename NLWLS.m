function [Xnom,P,H,W,y,bXnom,rvec] = NLWLS(dataset,variance,lat,long,alt,od)
    
    %{
        Non-Linear Weighted Least Squares
        This function applies linearization to a non-linear problem using 
        least squares method.
        
        Input:
            dataset  - Mx9 table containing the following parameters for
            any given number of observations (variable location does 
            not matter as long as all memebers are present and naming
            convetion is adhered to):
                yyyy - Year     
                mm   - Month    
                dd   - Day      
                H    - Hour     (UT)
                M    - Minute   
                S    - Second   
                Az   - Azimuth (deg)
                el   - elevation (deg)
                rho  - slant range (km)
            variance - 1x3 table with the following parameters:
                rho  - slant range variance or error (km)
                Az   - Azimuth variance or error (deg)
                el   - elevation variance or error (deg)
            lat      - latitude of site (deg)
            long     - longitude of site (deg)
            alt      - altitude of site (km)
            od       - orbit determination
                       'hg' for Herrick-Gibs orbit determination
                            NOTE - observations must be less than 1 degree
                            apart
                       'g' for Gibbs orbit determination
                            NOTE - observations must be greater than 1 
                            degree apart
        Output:
            Xnom  - nominal state vector after k-iterations [km km/s]
            P     - variance/covariance matrix
            H     - mapping matrix
            W     - weight matrix (non-normalized)
            y     - residuals vector
            NOTE: all of these outputs are the end states
            bXnom - best guess of the nominal state vector prior to
            iteration
            
    %}

    %% check inputs
    if ~ismember(dataset.Properties.VariableNames,{'yyyy','mm','dd','H','M','S','Az','el','rho'})
        error('dataset missing %s parameter,',dataset.Properties.VariableNames{~ismember(dataset.Properties.VariableNames,{'yyyy','mm','dd','H','M','S','Az','el','r'})})
    end
    if ~ismember(variance.Properties.VariableNames,{'Az','el','rho'})
        error('variance missing %s parameter',variance.Properties.VariableNames{~ismember(variance.Properties.VariableNames,{'Az','el','rho'})})
    end
    if abs(lat) >= 360 % ****** THIS MIGHT BE WORNG
        error('latitude beyond scope of function')
    end
    if abs(long) >= 360
        error('latitude beyond scope of function')
    end
    if ~strcmp(od,'hg') | od ~= 'g'
        error("od must be 'hg' or 'g'.")
    end


    %% calculate preliminary data
    dates = [dataset.yyyy dataset.mm dataset.dd dataset.H dataset.M dataset.S];

    len  = length(dataset.rho);
    lst  = zeros(len,1); % deg
    RA   = zeros(len,1); % deg
    DEC  = zeros(len,1); % deg
    rhoh = zeros(len,3); % km
    RS   = zeros(len,3); % km

    % local sidereal time
    for i = 1:len
        lst(i) = LST(dates(i,:),long);
    end

    % right ascension and declination
    for i = 1:len
        [alpha, delta] = RA_DEC_from_Az_el(lat, dataset.Az(i), dataset.el(i), lst(i), 'd');
        RA(i)  = alpha;
        DEC(i) = delta;
    end

    % slant range vector in ECI
    for i = 1:len
        rhoh(i,:) = rhohat(RA(i),DEC(i),'d')';
    end

    % site vector in ECI
    for i = 1:len
        RS(i,:) = RS_ECI(deg2rad(lat),alt,deg2rad(lst(i)))';
    end

    % distance vector in ECI
    rvec = rhoh.*dataset.rho + RS; 

    % normalized weight matrix
    W    = diag([1/variance.Az^2 1/variance.el^2 1/variance.rho^2]);
    WN   = W/norm(W);
    trudata = [dataset.Az dataset.el dataset.rho];
  
        
    %% get all 'middle' velocity vectors
    midnum = floor(len - 2);
    vmid   = zeros(midnum,3);
    ct     = 2;
    if strcmp(od,'hg')
        for i = 1:midnum
            vmid(i,:) = hgibbs(rvec(ct-1,:),rvec(ct,:),rvec(ct+1,:),juliandate(dates(ct-1,:)),juliandate(dates(ct,:)),juliandate(dates(ct+1,:)));
            ct = ct + 1;
        end
    else
        for i = 1:midnum
            vmid(i,:) = gibbs(rvec(ct-1,:),rvec(ct,:),rvec(ct+1,:));
            ct = ct + 1;
        end
    end

    % backward propagate all 'middle' vectors to t1
    ct      = 2;
    tmid    = zeros(midnum, 2);
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    for i = 1:midnum
        tmid(i,:) = [0 etime(dates(1,:),dates(ct,:))];
        ct = ct + 1;
    end

    rbackp = zeros(1,3);
    vbackp = zeros(1,3);
    ct     = 2;
    for i = 1:midnum
        [~, s1] = ode45(@orbitprop,tmid(i,:),[rvec(ct,:) vmid(i,:)],options,398600);
        rbackp(i,:) = s1(end,1:3);
        vbackp(i,:) = s1(end,4:6);
        ct = ct + 1;
    end

    % average state vectors to get nominal state at t1
    Xnom  = [sum(rbackp)/midnum sum(vbackp)/midnum];
    bXnom = Xnom;


    %% iterate nominal to get true solution
    n       = 3; % number of observing variables
    numobs  = length(rvec); % number of observations
    RMS_old = 2;
    RMS_new = 1;
    TOL     = 1.e-3;
    MAXITR  = 10e3;
    NUMITR  = 1;
    
    while abs((RMS_old - RMS_new)/RMS_old) > TOL && NUMITR < MAXITR

        H    = zeros(3,6);
        HTWH = zeros(6,6);
        HTWy = zeros(6,1);

        for obs = 1:numobs
            % get observations from state vector (nominal)
            [nomAz,nomel,nomrho] = sv2razel(Xnom(1:3),RS(obs,:),lst(obs),lat);
            % calculate residuals
            y = [trudata(obs,1) - nomAz; trudata(obs,2) - nomel; trudata(obs,3) - nomrho];

            if obs == 1
                tempXnom = Xnom;
            else
                % propagate nominal state to current time step
                [~,sv] = ode45(@orbitprop,[0  etime(dates(obs,:),dates(1,:))],Xnom,options,398600);
                % build propagated nominal state
                tempXnom = sv(end,:);
            end

            % build mapping matrix (change one element of the state
            for modf = 1:length(Xnom) 
                delta = tempXnom(modf)*0.001;
                Xmod  = tempXnom;
                Xmod(modf) = Xmod(modf) + delta;
                [modAz,model,modrho] = sv2razel(Xmod(1:3),RS(obs,:),lst(obs),lat);
                H(:,modf) = ([modAz;model;modrho] - [nomAz;nomel;nomrho])/delta;
            end
            
            % accumulate respective matrices for each time step
            HTWH = HTWH +  H'*WN*H;
            HTWy = HTWy + H'*WN*y;            
        end

        % determine if covariance matrix is singular
        if det(HTWH) <= eps % zero equivalent
            [U,S,V] = svd(HTWH);
            P       = inv(V*pinv(S)*U');
            Xhat    = P*HTWy;
        else
            P       = inv(HTWH);
            Xhat    = P*HTWy;
        end

        % update nominal state vector
        Xnom = Xnom + Xhat';
        % recalculate RMS
        RMS_old = RMS_new;
        RMS_new = sqrt(y'*W*y/(n*numobs));

        NUMITR = NUMITR + 1;
    end


    if NUMITR >= MAXITR
        warning('Maximum number of iterations reached.')
        warning('Solution did not converge.')
    end


end