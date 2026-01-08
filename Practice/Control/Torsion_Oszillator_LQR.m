% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      04.01.2026
% Beschreibung:
% In diesem Skript werde ich eine Zustandsregelung für den
% Torsionsschwinger auf basis von DR entwerfen, um dies später für ILC
% nutzen zu können.
% -------------------------------------------------------------

clc
clear
close all

%% Simulation (uncontrolled system)
% Simulation parameters
x0 = [0; 0; 0; 0];
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

% Input signal
t_step = 2;
u_inp = zeros(size(t_vec));
u_inp(t_vec >= t_step) = 2;

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', 1e-8, ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear(t, x, u_inp, t_vec), t_vec, x0, opts);
phi1 = x_sim(:, 1);
phi1_p = x_sim(:, 2);
phi2 = x_sim(:, 3);
phi2_p = x_sim(:, 4);

%% Controller design (state feedback control)
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

% LQR weighting matrices (as in DR)
Q_LQR = diag([1, 1, 10, 1]);
R_LQR = 1;

% LQR gain (contin.)
k_T = lqr(A, b, Q_LQR, R_LQR);

% Eigenvalues
eig_cont = eig(A - b*k_T);
fprintf('Eigenwerte des geregelten Systems (A_cont)\n');
disp(eig_cont);

%% Controller design discrete (for comparision)
% Discrete System
sys_contin = ss(A, b, c_T, 0);
sys_disc = c2d(sys_contin, Ts, 'zoh');

% System matrices
Ad = sys_disc.A;
bd = sys_disc.B;
c_Td = sys_disc.C;
dd = sys_disc.D;

% LQR gain
k_T_disc = dlqr(Ad, bd, Q_LQR, R_LQR);

% Controlled system dynamics
Ad_cont = Ad - bd * k_T_disc;
sys_disc_cont = ss(Ad_cont, bd, c_Td, dd, Ts);

% Print
fprintf('\nVergleich der Verstärkungsmatrizen (kontinuierlich vs. diskret):\n');
disp(k_T);
disp(k_T_disc);

%% Simulation (controlled System continuously)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR_cont(t, x, u_inp, t_vec), t_vec, x0, opts);
phi1_cont = x_sim(:, 1);
phi1_p_cont = x_sim(:, 2);
phi2_cont = x_sim(:, 3);
phi2_p_cont = x_sim(:, 4);

%% Simulation (controlled System discrete)
% Simulation
[~, ~, x_sim] = lsim(sys_disc_cont, u_inp(:), t_vec(:), x0);
phi1_contd = x_sim(:, 1);
phi1_p_contd = x_sim(:, 2);
phi2_contd = x_sim(:, 3);
phi2_p_contd = x_sim(:, 4);

%% Reference trajectory
% Load trajectory file
filename = 'Trajectory_01.mat';
filepath = fullfile(pwd, '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
load(filepath);

%% Feedforward control
% Calculate stationary control signal (contin)
u_ff = - ref_traj.phi2 ./ (c_T*inv((A-b*k_T))*b);

% Calculate stationary control signal (disc)
u_ffd = ref_traj.phi2 ./ dcgain(sys_disc_cont);

%% Simulation (feedforward controlled System continuously)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR_cont(t, x, u_ff, t_vec), t_vec, x0, opts);
phi1_cont_ff = x_sim(:, 1);
phi1_p_cont_ff = x_sim(:, 2);
phi2_cont_ff = x_sim(:, 3);
phi2_p_cont_ff = x_sim(:, 4);

%% Simulation (feedforward controlled System discrete)
% Simulation
[~, ~, x_sim] = lsim(sys_disc_cont, u_ffd(:), t_vec(:), x0);
phi1_contd_ff = x_sim(:, 1);
phi1_p_contd_ff = x_sim(:, 2);
phi2_contd_ff = x_sim(:, 3);
phi2_p_contd_ff = x_sim(:, 4);

%% Dynamic feedforward control (continuously)
% Controled dynamics
A_cont = A - b*k_T;
b_cont = b;
sys_contin_cont = ss(A_cont, b_cont, c_T, 0);

% Transfer function
G_contin_cont = tf(sys_contin_cont);

% Control input (stop at second derivative)
u_dff = (G_contin_cont.Denominator{1}(3).*ref_traj.phi2_pp + G_contin_cont.Denominator{1}(4).*ref_traj.phi2_p + G_contin_cont.Denominator{1}(5).*ref_traj.phi2) ./ G_contin_cont.Numerator{1}(5);

%% Simulation (dynamic feedforward controlled System continuously)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR_cont(t, x, u_dff, t_vec), t_vec, x0, opts);
phi1_cont_dff = x_sim(:, 1);
phi1_p_cont_dff = x_sim(:, 2);
phi2_cont_dff = x_sim(:, 3);
phi2_p_cont_dff = x_sim(:, 4);

%% Plot 1
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   
plot(t_vec, phi1, LineWidth=1, DisplayName='phi1 contin'); hold on;
plot(t_vec, phi2, LineWidth=1, DisplayName='phi2 contin');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (uncontrolled System)');
legend('Location', 'best');

subplot(2,2,2);  
plot(t_vec, phi1_p, LineWidth=1, DisplayName='phi1p contin'); hold on;
plot(t_vec, phi2_p, LineWidth=1, DisplayName='phi2p contin');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (uncontrolled System)');
legend('Location', 'best');

subplot(2,2,3);   
% plot(t_vec, phi1_cont, LineWidth=1, DisplayName='phi1'); hold on;
plot(t_vec, phi2_cont, LineWidth=1, DisplayName='phi2 contin'); hold on;
plot(t_vec, phi2_contd, LineWidth=1, DisplayName='phi2 disc'); 
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (controlled system)');
legend('Location', 'best');

subplot(2,2,4);  
% plot(t_vec, phi1_p_cont, LineWidth=1, DisplayName='phi1p'); hold on;
plot(t_vec, phi2_p_cont, LineWidth=1, DisplayName='phi2p contin'); hold on;
plot(t_vec, phi2_p_contd, LineWidth=1, DisplayName='phi2p disc');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (controlled system)');
legend('Location', 'best');

%% Plot 2
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   
plot(t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_cont_ff, LineWidth=1, DisplayName='phi2 contin');
plot(t_vec, phi2_contd_ff, LineWidth=1, DisplayName='phi2 disc');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (static feedforward control)');
legend('Location', 'best');

subplot(2,2,2);  
plot(t_vec, ref_traj.phi2_p, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_p_cont_ff, LineWidth=1, DisplayName='phi2p contin');
plot(t_vec, phi2_p_contd_ff, LineWidth=1, DisplayName='phi2p disc');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (static feedforward control)');
legend('Location', 'best');

subplot(2,2,3);   
plot(t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_cont_dff, LineWidth=1, DisplayName='phi2 contin');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (dynamic feedforward control)');
legend('Location', 'best');

subplot(2,2,4);  
plot(t_vec, ref_traj.phi2_p, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_p_cont_dff, LineWidth=1, DisplayName='phi2p contin');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (dynamic feedforward control)');
legend('Location', 'best');

%% Local functions
function dx = torsion_oszillator_linear(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    J1  = 0.0299;    % kgm^2
    J2  = 0.0299;    % kgm^2
    c_phi = 7.309;   % Nm/rad
    d_v1 = 0.055;    % Nms/rad
    d_v2 = 0.0064;   % Nms/rad

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % State space
    A = [0, 1, 0, 0;
        -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
        0, 0, 0, 1;
        c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
    b = [0;
        1/J1;
        0;
        0];
    
    % Dynamics
    dx = A*x_vec + b*u;
end

function dx = torsion_oszillator_linear_LQR_cont(t, x_vec, u_vec, t_vec)
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

    % Control law (cont.)
    k_T = [12.524133472585133, 1.268619349231718, -9.207508682229756, 0.314246813584626];
    u = u_in - k_T*x_vec;
    
    % Dynamics
    dx = A*x_vec + b*u;
end