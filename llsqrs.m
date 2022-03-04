function [xhat,P] = llsqrs(xo,yo)
    %{
        This funcion uses linear unweighted least squares to estimate a
        linear function through a set of given observation points.

        Input:
            xo - observed values of the independent variable
                    column vector
            yo - observed values of the dependent variable
                    column vector

        Output:
            xhat(1) - intercept of the line
            xhat(2) - slope of the line
            P       - covariance matrix
    %}

    % check input parameters
    if length(xo) < 3 || length(yo) < 3
        error('llsqrs meant to be used with more than 3 observations.')
    end
    if length(xo) ~= length(yo)
        error('xo and yo must have the same number of observations.')
    end

    H     = [ones(length(xo),1) xo];
    P     = inv(H'*H);
    xhat  = P*H'*yo;

end