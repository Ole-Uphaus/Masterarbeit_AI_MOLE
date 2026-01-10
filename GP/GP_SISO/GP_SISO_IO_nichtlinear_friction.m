% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      18.11.2025
% Beschreibung:
% In diesem skript werde ich untersuchen, wie der Gauß Prozess das
% nichtlineare Systemmodell inklusive einer ebenfalls nichtlinearen
% Stribeck Reibkraft erlernen kann.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', 'ILC', 'Simulation', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', 'System_Models');
addpath(ILC_path);
addpath(Model_Path);

%% Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

N = length(t_vec);
m_delay = 1;

% Input trajectorys
u_scale_train = [1, 3];
N_traj = length(u_scale_train);
u_scale_test = 0;

u_vec_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    u_vec_train_cell{i} = u_scale_train(i)*sin(2*pi/T_end.*t_vec');
end

% u_vec_test= u_scale_test*sin(2*pi/T_end.*t_vec');
u_vec_test = u_scale_test * (t_vec' >= 1);

%% Data generation
% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation train/test without friction
x0 = [0;
    0];

y_sim_train_cell = cell(N_traj, 1);
x_sim_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    [~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_cell{i} = x_sim(:, 1);
    x_sim_train_cell{i} = x_sim;
end

[~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_vec_test, t_vec), t_vec, x0, opts);
y_sim_test = x_sim(:, 1);

% Simulation train/test with friction
y_sim_train_fric_cell = cell(N_traj, 1);
x_sim_train_fric_cell = cell(N_traj, 1);
for i = 1:N_traj
    [~, x_sim] = ode45(@(t,x) oszillator_nonlinear_stribeck(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_fric_cell{i} = x_sim(:, 1);
    x_sim_train_fric_cell{i} = x_sim;
end

[~, x_sim] = ode45(@(t,x) oszillator_nonlinear_stribeck(t, x, u_vec_test, t_vec), t_vec, x0, opts);
y_sim_test_fric = x_sim(:, 1);

%% Predict System Dynamics
% Init Gaussian Processes
GP_IO = GP_SISO_IO();
GP_IO_fric = GP_SISO_IO();

% Train Gaussian Process
GP_IO.train_GP_model(y_sim_train_cell, u_vec_train_cell);
GP_IO_fric.train_GP_model(y_sim_train_fric_cell, u_vec_train_cell);

% Predict new Trajectory
[y_pred_test, y_std_test] = GP_IO.predict_trajectory_measurement(u_vec_test);
[y_pred_test_fric, y_std_test_fric] = GP_IO_fric.predict_trajectory_measurement(u_vec_test);

% Calculate prediction error
error_y = abs(y_pred_test - y_sim_test);
error_y_fric = abs(y_pred_test_fric - y_sim_test_fric);

%% Plots
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test, LineWidth=1, DisplayName='y-pred-test');
for i = 1:N_traj
    plot(t_vec, y_sim_train_cell{i}, LineWidth=1, DisplayName=sprintf('y-sim-train%d', i));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and Predicted Results (no friction)');
legend()

subplot(2,2,2);
plot(t_vec, y_sim_test_fric, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test_fric, LineWidth=1, DisplayName='y-pred-test');
for i = 1:N_traj
    plot(t_vec, y_sim_train_fric_cell{i}, LineWidth=1, DisplayName=sprintf('y-sim-train%d', i));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and Predicted Results (friction)');
legend()

subplot(2,2,3);
plot(t_vec, error_y, LineWidth=1, DisplayName='y-sim - y-pred'); hold on;
plot(t_vec, error_y_fric, LineWidth=1, DisplayName='y-sim-fric - y-pred-fric');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare prediction Error');
legend()

subplot(2,2,4);
plot(t_vec, u_vec_test, LineWidth=1, DisplayName='u-test'); hold on;
for i = 1:N_traj
    plot(t_vec, u_vec_train_cell{i}, LineWidth=1, DisplayName=sprintf('u-train%d', i));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('u [N]');
title('Training and Testing Input u');
legend()