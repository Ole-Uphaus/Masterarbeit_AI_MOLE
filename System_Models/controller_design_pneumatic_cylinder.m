function [S_gain] = controller_design_pneumatic_cylinder()
    %% Entwurf Optimalregelung (LQR) des mechanischen Teilsystems
    % Sample Time
    Ts = 0.001;

    % Parameter des Prüfstandes (aus Praktikumsanleitung)
    m  = 1.7;      % Masse des Zylinderkolbens [kg]
    
    % Zustandsraummodell
    A_contin = [0 1;
        0 0];
    b_contin = [0;
        1/m];
    c_T = [1, 0];
    
    % Diskretisierung
    sys_contin = ss(A_contin, b_contin, c_T, 0);
    sys_disc = c2d(sys_contin, Ts, 'zoh');
    
    % Diskretes Zustandsraummodell
    Ad = sys_disc.A;
    bd = sys_disc.B;
    c_Td = sys_disc.C;
    dd = sys_disc.D;
    
    % LQR Gewichtungsmatrizen
    Q_lqr = diag([1e5, 1]);
    r_lqr = 1;
    
    % Reglerverstärkung
    k_T = dlqr(Ad, bd, Q_lqr, r_lqr);
    
    % Dynamik des geregelten Systems
    Ad_cont = Ad - bd * k_T;
    sys_disc_cont = ss(Ad_cont, bd, c_Td, dd, Ts);
    
    % Statische Verstärkung berechnen
    S_gain = 1 / dcgain(sys_disc_cont);
end

