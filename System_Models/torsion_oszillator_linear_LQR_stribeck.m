function dx = torsion_oszillator_linear_LQR_stribeck(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    J1  = 0.0299;    % kgm^2
    J2  = 0.0299;    % kgm^2
    c_phi = 7.309;   % Nm/rad
    d_v1 = 0.055;    % Nms/rad
    d_v2 = 0.0064;   % Nms/rad

    % Friction Parameters
    Fc   = 3;   % Coulomb-Friction [Nm]
    Fs   = 6.5;     % Static Friction [Nm] (Fs > Fc)
    vs   = 2;     % Stribeck-Velocity [rad/s]
    v_eps = 0.001;   % Regularisation tanh()

    dv_1   = 0.06;     % viskos damping constant [Nms/rad]
    dv_2   = 0.006;     % viskos damping constant [Nms/rad]

    % Friction forces
    phi_1p = x_vec(2);
    phi_2p = x_vec(4);

    F_fric_1 = (Fc + (Fs - Fc)*exp(-(phi_1p/vs)^2)) * tanh(phi_1p/v_eps) + dv_1*phi_1p;
    F_fric_2 = dv_2*phi_2p;

    % Input
    u_in = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % State space
    A = [0, 1, 0, 0;
        -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
        0, 0, 0, 1;
        c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
    b = [0;
        1/J1;
        0;
        0];

    % Control Law
    k_T = [12.524133472585133, 1.268619349231718, -9.207508682229756, 0.314246813584626];
    u = u_in - k_T*x_vec;

    % Actuator limits
    u_max = 9;
    u_min = -9;

    % Saturation
    u = min(max(u, u_min), u_max);

    % Dynamics
    dx = A*x_vec + b*u - [0; F_fric_1/J1; 0; F_fric_2/J2];
end
