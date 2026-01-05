function dx = torsion_oszillator_linear_LQR(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    J1  = 0.031;    % kgm^2
    J2  = 0.237;    % kgm^2
    c_phi = 9;      % Nm/rad
    d_v1 = 0.070;   % Nms/rad
    d_v2 = 0.231;   % Nms/rad

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
    k_T = [7.004636887952207, 1.129661405169407, 2.995363112047798, 1.415299352286920];
    u = u_in - k_T*x_vec;
    
    % Dynamics
    dx = A*x_vec + b*u;
end

