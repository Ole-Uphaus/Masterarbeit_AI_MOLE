% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      26.10.2025
% Beschreibung:
% In diesem skript werde ich ILC Regelung für eine einfaches nichtlineares
% System entwerfen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', 'System_Models');
addpath(Model_Path);

%% Reference Trajectory
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Noise Parameters
sigma_v = 0.0;      % Measurement Noise 0.1
fc_v = 20;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
u_init = zeros(size(r_vec, 1), 1);

%% ILC PD-Type design
% Parameters
kp = 0.1;
kd = 100;
m_delay = 1;
N_iter = 10;
x0 = [0;
    0];

% Initialisation
ILC_PD = ILC_SISO(r_vec, m_delay, u_init);
ILC_PD.init_PD_type(kp, kd);

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_PD.u_vec; 0];
v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v);
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
for i = 1:N_iter
    % Update input
    u_sim = [ILC_PD.PD_update(y_sim); 0];

    % Simulate the system
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v);
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
end
% Calculate and log final error
ILC_PD.calculate_final_error(y_sim);
y_sim_pd = y_sim;

%% ILC quadratic optimal design
% Parameters
N = size(t_vec, 2);

W = eye(N-m_delay);
S = 0.5*eye(N-m_delay);

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init);
ILC_Quadr.init_Quadr_type(W, S)

% Track Results
y_cell_quadr = cell(N_iter+1, 1);
u_cell_quadr = cell(N_iter+1, 1);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v);
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
y_cell_quadr{1} = y_sim;
u_cell_quadr{1} = u_sim;
for i = 1:N_iter
    % Update input
    P = Lifted_dynamics_nonlinear_SISO(@(x) oszillator_linearized_discrete(x, Ts), N, m_delay, x_sim);
    u_sim = [ILC_Quadr.Quadr_update(y_sim, P); 0];

    % Simulate the system
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v);
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
    y_cell_quadr{i+1} = y_sim;
    u_cell_quadr{i+1} = u_sim;
end
% Calculate and log final error
ILC_Quadr.calculate_final_error(y_sim);

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
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

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
plot(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
% plot(0:(length(ILC_PD.RMSE_log)-1), ILC_PD.RMSE_log, LineWidth=1, DisplayName='ILC PD');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,2);
plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Noise');
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