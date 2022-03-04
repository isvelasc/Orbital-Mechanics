function TA = KeplersEq(Me, ecc)

% ==========================================
%
% This function uses Newton's Method to solve for Kepler's equation
%
%   x0     - Inital guess
%   TOL    - Tolerance
%   F      - Function for Newton
%   FPRIME - Derivative for Newton
%   THETA  - Function to for True Anomaly (rad)
%   Me     - Mean Anomaly (rad)
%   ecc    - eccentricity
%
% ==========================================

x0      = EGuess(Me, ecc);
TOL     = 1.e-8;
F       = @(Me, E, ecc) Me - E + ecc*sin(E);
FPRIME  = @(E, ecc) -1 + ecc*cos(E);
THETA   = @(E, e) atan(tan(E/2)*sqrt((1+e)/(1-e)))*2;

MAX_ITR = 12;
x2      = x0 - F(Me, x0, ecc)/FPRIME(x0, ecc);
itr     = 1;
error   = abs( x2 - x0 );

while error > TOL && itr < MAX_ITR
    x1    = x2;
    x2    = x1 - F(Me, x1, ecc)/FPRIME(x1, ecc);
    error = abs(x2 - x1);
    itr   = itr + 1;   
end

TA = THETA(x2, ecc);
end