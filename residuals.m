function r = residuals(xo,yo,yc)

    %{
        This function calculates the errors, residuals or difference
        between actual and calculated values for llsqrs.
        Input:
            xo - observed values of the independent variable
                    column vector
            yo - observed values of the dependent variable
                    column vector
            yc(1) - intercept of the line
            yc(2) - slope of the line

        Output:
            r     - errors
    %}
    
    r = yo - (yc(1) + yc(2).*xo);
end