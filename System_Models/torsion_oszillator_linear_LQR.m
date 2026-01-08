function dx = torsion_oszillator_linear_LQR(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    J1  = 0.0299;    % kgm^2
    J2  = 0.0299;    % kgm^2
    c_phi = 7.309;   % Nm/rad
    d_v1 = 0.055;    % Nms/rad
    d_v2 = 0.0064;   % Nms/rad

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

    % Control law
    k_T = [12.524133472585133, 1.268619349231718, -9.207508682229756, 0.314246813584626];
    u = u_in - k_T*x_vec;
    
    % Dynamics
    dx = A*x_vec + b*u;
end

