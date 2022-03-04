function rho = rhohat(alpha,delta,option)

    %{
        This function defines the distance to an object 

        Input:
            alpha - right ascension
            delta - declination
            option - 'd' for degrees, 'r' for radians
        Output:
            rho   - slant range vector
    %}

    if option == 'd'
        rho = [cosd(delta)*cosd(alpha); cosd(delta)*sind(alpha); sind(delta)];
    elseif option == 'r'
        rho = [cos(delta)*cos(alpha); cos(delta)*sin(alpha); sin(delta)];
    else
        error("Option invalid. 'd' or 'r' for degrees or radians")
    end

end