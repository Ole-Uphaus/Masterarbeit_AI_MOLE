function [sys_contin, sys_disc ,sys_disc_cont, Ts, S_gain, k_T_disc] = controller_design_torsion_oszillator()
    %% Controller design (discrete state feedback DR)
    % Sample Time
    Ts = 0.001;
    
    % Simulation parameters
    J1  = 0.0299;    % kgm^2
    J2  = 0.0299;    % kgm^2
    c_phi = 7.309;   % Nm/rad
    d_v1 = 0.055;    % Nms/rad
    d_v2 = 0.0064;   % Nms/rad
    
    % State space
    A = [0, 1, 0, 0;
        -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
        0, 0, 0, 1;
        c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
    
    b = [0;
        1/J1;
        0;
        0];
    
    c_T = [0, 0, 1, 0];
    
    d = 0;
    
    % Discrete System
    sys_contin = ss(A, b, c_T, 0);
    sys_disc = c2d(sys_contin, Ts, 'zoh');
    
    % System matrices
    Ad = sys_disc.A;
    bd = sys_disc.B;
    c_Td = sys_disc.C;
    dd = sys_disc.D;
    
    % LQR weighting matrices (as in DR)
    Q_LQR = diag([1, 1, 10, 1]);
    R_LQR = 1;
    
    % LQR gain
    k_T_disc = dlqr(Ad, bd, Q_LQR, R_LQR);
    
    % Controlled system dynamics
    Ad_cont = Ad - bd * k_T_disc;
    sys_disc_cont = ss(Ad_cont, bd, c_Td, dd, Ts);
    
    % Static gain (feedforward control signal)
    S_gain = 1 / dcgain(sys_disc_cont);
end

