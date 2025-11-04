% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      31.10.2025
% Beschreibung:
% In diesem skript werde ich die Klasse für den Gauß-Prozess evaluieren und
% erproben, sodass sie zukünftig auch für AI_MOLE eingesetzt werden kann.
% Dazu werde ich versuchen, das Systemverhalten eines linearen dynmiaschen
% repetetiven Systems durch einen Gauß Prozess zu erlernen und anschließend
% Prädiktionen durchzuführen.
% -------------------------------------------------------------

clc
clear
close all

%% Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

% Input trajectorys
u_scale_train = [1];
N_traj = length(u_scale_train);
u_scale_test = -1;

u_vec_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    u_vec_train_cell{i} = u_scale_train(i)*sin(2*pi/T_end.*t_vec');
end

u_vec_test= u_scale_test*sin(2*pi/T_end.*t_vec');
% u_vec_test = u_scale_test * (t_vec' >= 1) .* (t_vec' <= 3);

%% Data generation
% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation train/test
x0 = [0;
    0];

y_sim_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    [~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_cell{i} = x_sim(:, 1);
end

[~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_test, t_vec), t_vec, x0, opts);
y_sim_test = x_sim(:, 1);

%% Predict System Dynamics
% Init Gaussian Process
GP_IO = GP_SISO_IO();

% Train Gaussian Process
GP_IO.train_GP_model(y_sim_train_cell, u_vec_train_cell);

% Predict new Trajectory
[y_pred_test, y_std_test] = GP_IO.predict_trajectory(u_vec_test);

%% Linearize GP at given input trajectory
% Linearisation
P = GP_IO.linearize_at_given_trajectory(u_vec_train_cell{1});

%% Plot
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test, LineWidth=1, DisplayName='y-pred-test');
for i = 1:N_traj
    plot(t_vec, y_sim_train_cell{i}, LineWidth=1, DisplayName=sprintf('y-sim-train%d', i));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and Predicted Results');
legend()

subplot(1,2,2);
plot(t_vec, u_vec_test, LineWidth=1, DisplayName='u-test'); hold on;
for i = 1:N_traj
    plot(t_vec, u_vec_train_cell{i}, LineWidth=1, DisplayName=sprintf('u-train%d', i));
end

grid on;
xlabel('Zeit [s]'); 
ylabel('u [N]');
title('Training and Testing Input u');
legend()

%% Local functions
function dx = oszillator_linear(t, x_vec, u_vec, t_vec)
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

    % Dynamics
    dx = A*x_vec + B*u;
end