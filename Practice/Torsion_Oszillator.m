% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      18.11.2025
% Beschreibung:
% In diesem Skript werde ich das Systemmodell des Torsionsschwingers
% untersuchen, sodass ich es später für AI-MOLE nutzen kann. Dabei werde
% ich sowohl das ungeregelte als auch das geregelte System untersuchen.
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

%% Reference trajectory
% Load trajectory file
filename = 'Trajectory_01.mat';
filepath = fullfile(pwd, '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oscillator', 'Reference_Trajectories', filename);
load(filepath);

r_vec = ref_traj.phi2;

%% Simulation (controlled system)
% New state vector (because of the additional integrator state)
x0 = [0; 0; 0; 0; 0];

% Simulation
[~, x_sim] = ode45(@(t,x) torsion_oszillator_linear_PI(t, x, r_vec, t_vec), t_vec, x0, opts);
phi1_cont = x_sim(:, 1);
phi1_p_cont = x_sim(:, 2);
phi2_cont = x_sim(:, 3);
phi2_p_cont = x_sim(:, 4);

%% Stability analysis
% Simulation parameters
J1  = 0.031;    % kgm^2
J2  = 0.237;    % kgm^2
c_phi = 9;      % Nm/rad
d_v1 = 0.070;   % Nms/rad
d_v2 = 0.231;   % Nms/rad

% Controller parameters
Kp = 3;
KI = 1;

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

% Eigenvalues
eig_uncont = eig(A);
eig_cont = eig(A_cont);

% Print
fprintf('Eigenwerte des ungeregelten Systems (A)\n');
disp(eig_uncont);
fprintf('Eigenwerte des geregelten Systems (A_cont)\n');
disp(eig_cont);

%% Plots
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
plot(t_vec, r_vec, LineWidth=1, DisplayName='reference');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi [rad]');
title('Simulation results (feedback controller)');
legend('Location', 'best');

subplot(2,2,4);  
% plot(t_vec, phi1_p_cont, LineWidth=1, DisplayName='phi1p'); hold on;
plot(t_vec, phi2_p_cont, LineWidth=1, DisplayName='phi2p');
grid on;
xlabel('Zeit [s]'); 
ylabel('phip [rad/s]');
title('Simulation results (feedback controller)');
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