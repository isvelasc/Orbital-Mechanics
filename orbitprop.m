function dstate = orbitprop(time, state, mu)

    %{
        This function numerically computes the state vector r, v as a
            function of time, given initial values r0, v0. This function
            should be used with a differential equation solver.
        Input:
            time   - time that the solver will use to solver for each step
            state  - 1x6 column vector
                position - 1x 2y 3z of the vector (km)
                velocity - 4x 5y 6z of the vector (km/s)
            mu     - gravitational parameter (km^3/s^2)
        Output:
            dstate - the new state of the vecotr after being propagated
    %}
    
    dstate  = [state(4:6);-mu*state(1:3)/norm(state(1:3))^3];
end