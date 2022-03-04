function [r2vec,v2vec] = gauss_angles_only(rhohat1,rhohat2,rhohat3,Rs1,Rs2,Rs3,JD1,JD2,JD3,mu,R2,options)

    %{
        This function solves the initial orbit determination using
            Gauss' angles only technique. This function assumes that
            delta t between observation is short ( < 10 mins apart)
            and/or there is less than 10 degrees of separation between
            observartions.

        Input:
            ***All vector inputs must be column vectors***
            rhohat1 - Slant range vector at 1st observation
            rhohat2 - Slant range vector at 2nd observation
            rhohat3 - Slant range vector at 3rd observation
            Rs1     - Site location vector at 1st observation (km)
            Rs2     - Site location vector at 2nd observation (km)
            Rs3     - Site location vector at 3rd observation (km)
            JD1     - Juliandate at 1st observation (UTC time)
            JD2     - Juliandate at 2nd observation (UTC time)
            JD3     - Juliandate at 3rd observation (UTC time)
            options - struct comprised of 3 optional parameters
                         by convention named: tau1, tau3, e
                        tau1 = T1 - T2 (s)
                        tau3 = T3 - T2 (s)
                        e    = 1 (for extended) 0 (for non-extended)
                                **non-extended default**
            mu      - gravitational parameter (km^3/s^2)

        Output:
            r2 - (km)
            v2 - (km)
    %}

    % check for options
    if nargin > 9 && ~isempty(options)
        if isfield(options, 'tau1')
            T1 = options.tau1;
        else
            T1 = 0;
        end
        if isfield(options, 'tau3')
            T3 = options.tau3;
        else 
            T3 = 0;
        end
        if isfield(options,'e')
            e = options.e;
        else
            e = 0;
        end
    end

    if ~T1
        T1 = JD1 - JD2;
    end
    if ~T3
        T3 = JD3 - JD2;
    end

    a1    = T3/(T3 - T1);
    a1u   = (T3*((T3 - T1)^2 - T1^2))/(6*(T3 - T1));
    a3    = -T1/(T3 - T1);
    a3u   = -(T1*((T3 - T1)^2 - T1^2))/(6*(T3 - T1));
    L     = [rhohat1 rhohat2 rhohat3];
    Linv  = LMatrix(rhohat1,rhohat2,rhohat3)/norm(L);
    rsite = [Rs1 Rs2 Rs3];
    M     = Linv*rsite;
    d1    = M(2,1)*a1 - M(2,2) + M(2,3)*a3;
    d2    = M(2,1)*a1u + M(2,3)*a3u;
    C     = dot(rhohat2,Rs2);

    % solve for the real root of r2
%     syms R2
    eqn   = R2^8 - (d1^2 + 2*C*d1 + norm(Rs2)^2)*R2^6 -...
        2*mu*(C*d2 + d1*d2)*R2^3 - mu^2*d2^2 == 0;
    r2    = vpasolve(eqn,R2);
    r2    = real(double(r2(8)));

    u     = mu/(r2^3);
    c1    = a1 + a1u*u;
    c2    = -1;
    c3    = a3 + a3u*u;
    mcvec = M*[-c1;-c2;-c3];
    rho1  = mcvec(1,:)/c1;
    rho2  = mcvec(2,:)/c2;
    rho3  = mcvec(3,:)/c3;

    r2vec = norm(rho2)*rhohat2 + Rs2;

    % find lagrange coefficients
    f1    = 1 - 0.5*mu*T1^2/r2^3;
    f3    = 1 - 0.5*mu*T3^2/r2^3;
    g1    = T1 - (1/6)*mu*T1^2/r2^3;
    g3    = T3 - (1/6)*mu*T3^2/r2^3;

    v2vec = 1/(f1*g3 - f3*g1)*(-f3*rho1 + f1*rho3);

    if e % extended
        
        r2vec = 0; v2vec = 0;
        error('not implemented')
    end



    function lm = LMatrix(rho1, rho2, rho3)
        lm = [ rho2(2)*rho3(3) - rho3(2)*rho2(3) -rho1(2)*rho3(3) + rho3(2)*rho1(3)  rho1(2)*rho2(3) - rho2(2)*rho1(3);...
              -rho2(1)*rho3(3) + rho3(1)*rho2(3)  rho1(1)*rho3(3) - rho3(1)*rho1(3) -rho1(1)*rho2(3) + rho2(1)*rho1(3);...
               rho2(1)*rho3(2) - rho3(1)*rho2(2) -rho1(1)*rho3(2) + rho3(1)*rho1(2)  rho1(1)*rho2(2) - rho2(1)*rho1(2)];
    end
end