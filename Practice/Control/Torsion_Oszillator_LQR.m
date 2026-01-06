% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      04.01.2026
% Beschreibung:
% In diesem Skript werde ich eine Zustandsregelung für den
% Torsionsschwinger auf basis von ASK entwerfen, um dies später für ILC
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
J1  = 0.031;    % kgm^2
J2  = 0.237;    % kgm^2
c_phi = 9;      % Nm/rad
d_v1 = 0.070;   % Nms/rad
d_v2 = 0.231;   % Nms/rad

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

% LQR weighting matrices (as in ASK)
Q_LQR = diag([50, 1, 50, 1]);
R_LQR = 1;

% LQR gain
k_T = lqr(A, b, Q_LQR, R_LQR);

% Eigenvalues
eig_cont = eig(A - b*k_T);
fprintf('Eigenwerte des geregelten Systems (A_cont)\n');
disp(eig_cont);

%% Controller design discrete (for comparision)
% Discrete System
Ts = 0.001;
sys_contin = ss(A, b, c_T, 0);
sys_disc = c2d(sys_contin, Ts, 'zoh');

% System matrices
Ad = sys_disc.A;
bd = sys_disc.B;

% LQR gain
k_T_disc = dlqr(Ad, bd, Q_LQR, R_LQR);

% Print
fprintf('\nVergleich der Verstärkungsmatrizen (kontinuierlich vs. diskret):\n');
disp(k_T);
disp(k_T_disc);

%% Simulation (controlled System)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR(t, x, u_inp, t_vec), t_vec, x0, opts);
phi1_cont = x_sim(:, 1);
phi1_p_cont = x_sim(:, 2);
phi2_cont = x_sim(:, 3);
phi2_p_cont = x_sim(:, 4);

%% Reference trajectory
% Load trajectory file
filename = 'Trajectory_01.mat';
filepath = fullfile(pwd, '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
load(filepath);

%% Feedforward control
% Calculate stationary control signal
u_ff = - ref_traj.phi2 ./ (c_T*inv((A-b*k_T))*b);

%% Simulation (feedforward controlled System)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR(t, x, u_ff, t_vec), t_vec, x0, opts);
phi1_cont_ff = x_sim(:, 1);
phi1_p_cont_ff = x_sim(:, 2);
phi2_cont_ff = x_sim(:, 3);
phi2_p_cont_ff = x_sim(:, 4);

%% Dynamic feedforward control
% Controled dynamics
A_cont = A - b*k_T;
b_cont = b;
sys_cont = ss(A_cont, b_cont, c_T, 0);

% Transfer function
G_cont = tf(sys_cont);

% Control input (stop at second derivative)
u_dff = (G_cont.Denominator{1}(3).*ref_traj.phi2_pp + G_cont.Denominator{1}(4).*ref_traj.phi2_p + G_cont.Denominator{1}(5).*ref_traj.phi2) ./ G_cont.Numerator{1}(5);

%% Simulation (dynamic feedforward controlled System)
% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_LQR(t, x, u_dff, t_vec), t_vec, x0, opts);
phi1_cont_dff = x_sim(:, 1);
phi1_p_cont_dff = x_sim(:, 2);
phi2_cont_dff = x_sim(:, 3);
phi2_p_cont_dff = x_sim(:, 4);

%% Plot 1
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   
plot(t_vec, phi1, LineWidth=1, DisplayName='phi1'); hold on;
plot(t_vec, phi2, LineWidth=1, DisplayName='phi2');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (uncontrolled System)');
legend('Location', 'best');

subplot(2,2,2);  
plot(t_vec, phi1_p, LineWidth=1, DisplayName='phi1p'); hold on;
plot(t_vec, phi2_p, LineWidth=1, DisplayName='phi2p');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (uncontrolled System)');
legend('Location', 'best');

subplot(2,2,3);   
% plot(t_vec, phi1_cont, LineWidth=1, DisplayName='phi1'); hold on;
plot(t_vec, phi2_cont, LineWidth=1, DisplayName='phi2'); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (controlled system)');
legend('Location', 'best');

subplot(2,2,4);  
% plot(t_vec, phi1_p_cont, LineWidth=1, DisplayName='phi1p'); hold on;
plot(t_vec, phi2_p_cont, LineWidth=1, DisplayName='phi2p'); hold on;
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
plot(t_vec, phi2_cont_ff, LineWidth=1, DisplayName='phi2');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (static feedforward control)');
legend('Location', 'best');

subplot(2,2,2);  
plot(t_vec, ref_traj.phi2_p, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_p_cont_ff, LineWidth=1, DisplayName='phi2p');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (static feedforward control)');
legend('Location', 'best');

subplot(2,2,3);   
plot(t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_cont_dff, LineWidth=1, DisplayName='phi2');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (dynamic feedforward control)');
legend('Location', 'best');

subplot(2,2,4);  
plot(t_vec, ref_traj.phi2_p, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, phi2_p_cont_dff, LineWidth=1, DisplayName='phi2p');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (dynamic feedforward control)');
legend('Location', 'best');

%% Local functions
function dx = torsion_oszillator_linear(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    J1  = 0.031;    % kgm^2
    J2  = 0.237;    % kgm^2
    c_phi = 9;      % Nm/rad
    d_v1 = 0.070;   % Nms/rad
    d_v2 = 0.231;   % Nms/rad

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