function [r2vec,v2vec] = IOD_gauss_CM(rhv1,rhv2,rhv3,Rs1,Rs2,Rs3,tau1,tau3,mu,e)

    %{
        Initial orbit determination using Gauss angle-only
            formulation with Curtis Method of implementation.
            This function assumes that delta t between observation
            is short ( < 10 mins apart) and/or there is 
            less than 10 degrees of separation between observartions.

        Input:
            rhv1 - Slant range vector at 1st observation
            rhv2 - Slant range vector at 2nd observation
            rhv3 - Slant range vector at 3rd observation
            Rs1  - Site location vector at 1st observation (km)
            Rs2  - Site location vector at 2nd observation (km)
            Rs3  - Site location vector at 3rd observation (km)
            tau1 - time difference at 1st observation (t1-t2)
            tau3 - time difference at 3rd observation (t3-t2)
            mu   - gravitational parameter (km^3/s^2)
            e    - extended option, 0 for non-extende (default)
                                    1 for extended

        Output:
            r2vec - middle position vector (km)
            v2vec - middle velocity vector (km/s)
    %}


    tau  = tau3 - tau1; 
    p1   = cross(rhv2,rhv3);
    p2   = cross(rhv1,rhv3);
    p3   = cross(rhv1,rhv2);
    D0   = dot(rhv1,p1);
    D    = DMatrix(Rs1,Rs2,Rs3,p1,p2,p3);
    A    = 1/D0*(-D(1,2)*tau3/tau + D(2,2) + D(3,2)*tau1/tau);
    B    = 1/(6*D0)*(D(1,2)*(tau3^2 - tau^2)*(tau3/tau) + D(3,2)*(tau^2 -tau1^2)*(tau1/tau));
    E    = dot(Rs2,rhv2);
    R2sq = dot(Rs2,Rs2);
    a    = -(A^2 + 2*A*E + R2sq);
    b    = -2*mu*B*(A + E);
    c    = -mu^2*B^2;

    % find roots
    rts  = roots([1 0 a 0 0 b 0 0 c]);
    rrts = rts(imag(rts)==0);
    r2   = rrts(rrts>=0);

    rho1 = 1/D0*((6*(D(3,1)*(tau1/tau3) + D(2,1)*(tau/tau3))*r2^3 +...
        mu*D(3,1)*(tau^2 - tau1^2)*(tau1/tau3))/...
        (6*r2^3 + mu*(tau^2 - tau3^2)) - D(1,1));
    rho2 = A + (mu*B)/r2^3;
    rho3 = 1/D0*((6*(D(1,3)*(tau3/tau1) - D(2,3)*(tau/tau1))*r2^3 +...
        mu*D(1,3)*(tau^2 - tau3^2)*(tau3/tau1))/...
        (6*r2^3 + mu*(tau^2 - tau1^2)) - D(3,3));

    r1vec = Rs1 + rho1*rhv1;
    r2vec = Rs2 + rho2*rhv2;
    r3vec = Rs3 + rho3*rhv3;

    % lagrange coefficients
    f1    = 1 - 0.5*mu*tau1^2/r2^3;
    g1    = tau1 - (1/6)*mu*tau1^3/r2^3;
    f3    = 1 - 0.5*mu*tau3^2/r2^3;    
    g3    = tau3 - (1/6)*mu*tau3^3/r2^3;

    v2vec = 1/(f1*g3 - f3*g1)*(-f3*r1vec + f1*r3vec);
    

    % extended version
    if e
        
        prev_rho1 = rho1;
        prev_rho2 = rho2;
        prev_rho3 = rho3;
        ERR1      = 1;
        ERR2      = 1;
        ERR3      = 1;
        iter      = 0;
        max_iter  = 1000;
        TOL       = 1.e-8;

        while (ERR1 > TOL && ERR2 > TOL && ERR3 > TOL) && (iter < max_iter)
            nr2 = norm(r2vec);

            [chi1, ~]     = Kepler_UV(r2vec, v2vec, tau1, mu);
            [chi3, alpha] = Kepler_UV(r2vec, v2vec, tau3, mu);
            [nf1,ng1]     = FG_UV(alpha, chi1, nr2, tau1, mu);
            [nf3,ng3]     = FG_UV(alpha, chi3, nr2, tau3, mu);

            f1 = (f1 + nf1)/2;
            f3 = (f3 + nf3)/2;
            g1 = (g1 + ng1)/2;
            g3 = (g3 + ng3)/2;

            c1 = g3/(f1*g3 - f3*g1);
            c3 = -g1/(f1*g3 - f3*g1);

            rho1 = 1/D0*(-D(1,1) + 1/c1*D(2,1) - c3/c1*D(3,1));
            rho2 = 1/D0*(-c1*D(1,2) + D(2,2) - c3*D(3,2));
            rho3 = 1/D0*(-c1/c3*D(1,3) + 1/c3*D(2,3) - D(3,3));

            r1vec = Rs1 + rho1*rhv1;
            r2vec = Rs2 + rho2*rhv2;
            r3vec = Rs3 + rho3*rhv3;

            v2vec = 1/(f1*g3 - f3*g1)*(-f3*r1vec + f1*r3vec);

            ERR1  = abs(rho1 - prev_rho1);
            ERR2  = abs(rho2 - prev_rho2);
            ERR3  = abs(rho3 - prev_rho3);

            prev_rho1 = rho1;
            prev_rho2 = rho2;
            prev_rho3 = rho3;
        end

        if iter >= max_iter
            fprintf('Number of max iterations exceeded\n\n');
        end
    end



    function d = DMatrix(R1,R2,R3,p1,p2,p3)
        d = [dot(R1,p1) dot(R1,p2) dot(R1,p3);...
             dot(R2,p1) dot(R2,p2) dot(R2,p3);...
             dot(R3,p1) dot(R3,p2) dot(R3,p3)];
    end

end