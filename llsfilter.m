function [fxo,fyo,fpts] = llsfilter(xo,yo,r,RMS)

    %{
        This function filters out bad data points for the least squares
            method using the difference between the root mean square and
            rsiduals.

        Input:
            xo  - observed values of the independent variable
                    column vector
            yo  - observed values of the dependent variable
                    column vector
            r   - residuals
            RMS - root mean square

        Output:
            fxo  - filtered independent data points
            fyo  - filtered dependent data points
            fpts - filtered data points
    %}

    % check input parameters
    if length(xo) < 4 || length(yo) < 4
        error('llsfilter meant to be used with more than 4 observations.')
    end
    if length(xo) ~= length(yo) || length(xo) ~= length(r)
        error('xo, yo, and r must have the same length.')
    end
    if RMS < 0 || RMS == 0
        error('RMS must be greater than zero')
    end

    fxo  = zeros;
    fyo  = zeros;
    fpts = zeros;
    err  = 2*RMS;

    loc = 1;
    for obs = 1:length(xo)
        if abs(r(obs)) <= err
            fxo(loc) = xo(obs);
            fyo(loc) = yo(obs);
            loc = loc + 1;
        else
            fpts(end+1) = obs;
        end
    end
    fxo  = fxo';
    fyo  = fyo';
    fpts = fpts(2:end)';

end