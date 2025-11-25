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

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', 'System_Models');
addpath(Model_Path);

%% System Dynamics
% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
d  = 0.5; % Ns/m
m_delay = 1;

% Noise Parameters
sigma_w_rep = 0;  % Repeating process Noise (trial invariant)
sigma_w = 0;     % Process Noise 0.05
sigma_v = 0.0;      % Measurement Noise 0.1
fc_w = 0.1;
fc_v = 10;
white = true;       % if white == true -> white noise is sampled - no filter

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
R = 0.00001*eye(size(P));

% Track Results
y_cell_quadr = cell(N_iter+1, 1);
u_cell_quadr = cell(N_iter+1, 1);

% Init repeating process Noise
w_rep_vec = Gen_noise_Butter(t_vec, sigma_w_rep, fc_w, white);

% Q-Filter
Q_order = 2;
Q_fc = 2;

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init);
ILC_Quadr.init_Quadr_type(W, S, R, P)
% ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v, white);
[t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim + w_vec + w_rep_vec), t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
y_cell_quadr{1} = y_sim;
u_cell_quadr{1} = u_sim;
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v, white);
    [t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim + w_vec + w_rep_vec), t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
    y_cell_quadr{i+1} = y_sim;
    u_cell_quadr{i+1} = u_sim;
end
% Calculate and log final error
ILC_Quadr.calculate_final_error(y_sim);

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:N_iter
    plot(t_vec, y_cell_quadr{i}, LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, y_cell_quadr{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(2,2,3);
plot(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
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
hold on;
for i = 1:N_iter
    plot(t_vec, u_cell_quadr{i}, LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, u_cell_quadr{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
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
function [Ad, Bd, Cd, Dd] = linear_discrete_system(x_star, Ts)
%linear_discrete_system this function is just used for comparison of the
%lifted matrices. It is used as a local function, because the nonlinearity
%is set to zero to compare to linear system. Otherwise we could have used
%the global function.

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

function [w_vec, v_vec] = proc_meas_noise(t_vec, fc_w, fc_v, sigma_w, sigma_v, white)
% Calculate Noises
w_vec = Gen_noise_Butter(t_vec, sigma_w, fc_w, white);
v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
end