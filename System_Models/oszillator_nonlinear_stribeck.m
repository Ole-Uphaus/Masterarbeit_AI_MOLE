function dx = oszillator_nonlinear_stribeck(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 2; % N/m^3
    d  = 0.5; % Ns/m

    % Friction Parameters (Stribeck)
    Fc   = 0.4;     % Coulomb-Friction [N]
    Fs   = 0.5;     % Static Friction [N] (Fs > Fc)
    vs   = 0.1;     % Stribeck-Velocity [m/s]
    dv   = 0.1;     % viskos damping constant [Ns/m]
    v_eps = 0.001;   % Regularisation tanh()

    % States
    x = x_vec(1);
    xp = x_vec(2);

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % Stribeck-friction
    F_fric = (Fc + (Fs - Fc)*exp(-(xp/vs)^2)) * tanh(xp/v_eps) + dv*xp;
    % F_fric = 0;

    % Dynamics
    dx = zeros(2, 1);
    dx(1) = xp;
    dx(2) = 1/m*(-c1*x - c2*x^3 - d*xp - F_fric + u);
end

