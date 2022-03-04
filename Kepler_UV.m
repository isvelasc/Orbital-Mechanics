function [chi, alpha] = Kepler_UV(r0_vec, v0_vec, dt, mu)

    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %   This function uses Newton's method to solve the
    %       the universal Kepler equation for the
    %       universal anomaly Chi(χ)
    %
    %   Input:
    %       r0     - initial position vector (km)
    %       vr0    - initial velocity vector (km/s)
    %       dt     - elapsed time (s)
    %       mu     - gravitational parameter (km^3/s^2)
    %   Output:
    %       chi(χ) - universal anomaly (km^0.5)
    %       alpha  - reciprocal of the semimajor axis (1/km)
    %   
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    %   Postion and veloicty magnitudes
    r0 = norm(r0_vec);
    v0 = norm(v0_vec);
    
    %   Initial radial velocity
    vr0 = dot(r0_vec,v0_vec)/r0;
    
    %   Reciprocal of the semimajor axis
    alpha = 2/r0 - v0^2/mu;
    
    %   Set tolerance and limit to iterations
    TOL = 1.e-8;
    max_iter = 1000;
    
    %   Starting value for chi(χ)
    chi = sqrt(mu)*abs(alpha)*dt;
    
    %   Number of iterations for convergence
    n = 0;
    
    %   Breaking point for tolerance
    ratio = 1;
    
    %   Iterate until convergence or error on tolerance
    while abs(ratio) > TOL && n <= max_iter
        n     = n + 1;
        z     = alpha*chi^2;
        C     = stumpC(z);
        S     = stumpS(z);
        F     = r0*vr0/sqrt(mu)*chi^2*C + (1 - alpha*r0)*chi^3*S + r0*chi - sqrt(mu)*dt;
        dF    = r0*vr0/sqrt(mu)*chi*(1 - alpha*chi^2*S) + (1 - alpha*r0)*chi^2*C + r0;
        ratio = F/dF;
        chi   = chi - ratio;
    end
    
    %   Report max iterations were reached
    if n > max_iter
        fprintf('\n **No. iterations of Kepler’s equation = %g', n);
        fprintf('\n F/dF = %g\n', F/dF);
    end

end

