% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      28.10.2025
% Beschreibung:
% In diesem skript werde ich ILC Regelung für eine den linearen Oszillator
% nochmal zum Vergelich erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

%% System Dynamics
% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
d  = 0.5; % Ns/m
m_delay = 1;

% Noise Parameters
sigma_w_rep = 0;  % Repeating process Noise (trial invariant)
sigma_w = 0;     % Process Noise 0.05
sigma_v = 0;      % Measurement Noise 0.1
fc_w = 0.1;
fc_v = 10;

% State space representation 
A = [0, 1;
    -c1/m, -d/m];
B = [0;
    1/m];
C = [1, 0];
D = 0;

sys_cont = ss(A,B,C,D);
G_cont = tf(sys_cont);

% Discrete state space
Ts = 1.e-2;
sys_disc = c2d(sys_cont, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sys_disc);

%% Reference Trajectory
% Parameters
x_max = 0.5;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
u_init = zeros(size(r_vec, 1), 1);

%% ILC quadratic optimal design
% Lifted system dynamics
N = size(t_vec, 2);
P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);

% Parameters
N_iter = 10;
x0 = [0;
    0]; 
W = eye(size(P));
S = 0.001*eye(size(P));

% Init repeating process Noise
w_rep_vec = Gen_noise_Butter(t_vec, sigma_w_rep, fc_w);

% Q-Filter
Q_order = 2;
Q_fc = 2;

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init);
ILC_Quadr.init_Quadr_type(W, S, P)
% ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v);
[t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_sim, t_vec, (w_vec + w_rep_vec)), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v);
    [t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_sim, t_vec, (w_vec + w_rep_vec)), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
end
y_vec_Quadr = y_sim;
u_sim_Quadr = u_sim;

%% ILC PD-Type
% Parameters
kp = 0.1;
kd = 100;

% Initialisation
ILC_PD = ILC_SISO(r_vec, m_delay, u_init);
ILC_PD.init_PD_type(kp, kd);
% ILC_PD.init_Q_lowpass(Q_fc, Q_order, Ts);

% Update Loop
u_sim = [ILC_PD.u_vec; 0];
[w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v);
[t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_sim, t_vec, (w_vec + w_rep_vec)), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
for i = 1:N_iter
    % Update input
    u_sim = [ILC_PD.PD_update(y_sim); 0];

    % Simulate the system
    [w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v);
    [t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_sim, t_vec, (w_vec + w_rep_vec)), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
end
y_vec_PD = y_sim;
u_sim_PD = u_sim;

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, y_vec_Quadr, LineWidth=1, DisplayName='ILC Quadr');
plot(t_vec, y_vec_PD, LineWidth=1, DisplayName='ILC PD');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(2,2,3);
plot(1:length(ILC_Quadr.RMSE_log), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
plot(1:length(ILC_PD.RMSE_log), ILC_PD.RMSE_log, LineWidth=1, DisplayName='ILC PD');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,2);
plot(t_vec, w_vec, LineWidth=1, DisplayName='proc'); hold on;
plot(t_vec, w_rep_vec, LineWidth=1, DisplayName='proc rep');
plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare process and measurement noise');
legend()

subplot(2,2,4);
plot(t_vec, u_sim_Quadr, LineWidth=1, DisplayName='u quadr'); hold on;
% plot(t_vec, u_sim_PD, LineWidth=1, DisplayName='u pd');
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Input Signal');
legend()

%% Compare lifted represetation
% Parameters
x_sim = zeros(N, 2);

% Lifted Dynamics
P_nonlin = Lifted_dynamics_nonlinear_SISO(@(x) linear_discrete_system(x, Ts), N, m_delay, x_sim);

% Compare Results
delta_P = P - P_nonlin;
max_error_P = max(abs(P(:) - P_nonlin(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von P: %.3e\n', max_error_P);

%% Local Functions
function dx = oszillator_linear(t, x_vec, u_vec, t_vec, w_vec)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    d  = 0.5; % Ns/m

    % State space representation 
    A = [0, 1;
        -c1/m, -d/m];
    B = [0;
        1/m];

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');
    w = interp1(t_vec, w_vec, t, 'previous', 'extrap');

    % Dynamics
    dx = A*x_vec + B*(u + w);
end

function [Ad, Bd, Cd, Dd] = linear_discrete_system(x_star, Ts)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 0; % N/m^3
    d  = 0.5; % Ns/m

    % States
    x = x_star(1);
    xp = x_star(2);

    % Linearisation
    A_lin = [0, 1;
        (-c1/m - 3*c2/m*x^2), -d/m];
    B_lin = [0;
        1/m];
    C_lin = [1, 0];
    D_lin = 0;
    
    sys_cont = ss(A_lin, B_lin, C_lin, D_lin);
    
    % Discrete
    sys_disc = c2d(sys_cont, Ts, 'zoh');
    [Ad, Bd, Cd, Dd] = ssdata(sys_disc);
end

function [w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v)
% Calculate Noises
w_vec = Gen_noise_Butter(t_vec, sigma_w, fc_w);
v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v);
end