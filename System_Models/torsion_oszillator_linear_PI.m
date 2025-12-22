function dx = torsion_oszillator_linear_PI(t, x_vec, r_vec, t_vec)
    % Simulation parameters
    J1  = 0.031;    % kgm^2
    J2  = 0.237;    % kgm^2
    c_phi = 9;      % Nm/rad
    d_v1 = 0.070;   % Nms/rad
    d_v2 = 0.231;   % Nms/rad

    % Controller parameters
    Kp = 3;
    KI = 1;

    % Input (reference trajectory)
    r = interp1(t_vec, r_vec, t, 'previous', 'extrap');

    % State space
    A = [0, 1, 0, 0;
        -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
        0, 0, 0, 1;
        c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
    b = [0;
        1/J1;
        0;
        0];
    C = [0, 0, 1, 0];

    % State space (controlled system)
    A_cont = [A-Kp*b*C, KI*b;
        -C, 0];
    b_cont = [Kp*b;
        1];
     
    % Dynamics
    dx = A_cont*x_vec + b_cont*r;
end
