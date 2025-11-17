function dx = oszillator_linear(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    d  = 0.5; % Ns/m

    % State space representation 
    A = [0, 1;
        -c1/m, -d/m];
    B = [0;
        1/m];

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % Dynamics
    dx = A*x_vec + B*u;
end