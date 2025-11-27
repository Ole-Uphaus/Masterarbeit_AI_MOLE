% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      31.10.2025
% Beschreibung:
% In diesem skript werde ich die Klasse für den Gauß-Prozess evaluieren und
% erproben, sodass sie zukünftig auch für AI_MOLE eingesetzt werden kann.
% Dazu werde ich versuchen, das Systemverhalten eines nichtlinearen dynmiaschen
% repetetiven Systems durch einen Gauß Prozess zu erlernen und anschließend
% Prädiktionen durchzuführen.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', 'ILC', 'ILC_SISO');
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
u_scale_test = 2;

u_vec_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    u_vec_train_cell{i} = u_scale_train(i)*sin(2*pi/T_end.*t_vec');
end

u_vec_test= u_scale_test*sin(2*pi/T_end.*t_vec');

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
x_sim_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    [~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_cell{i} = x_sim(:, 1);
    x_sim_train_cell{i} = x_sim;
end

[~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_vec_test, t_vec), t_vec, x0, opts);
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
tic;
P = GP_IO.linearize_at_given_trajectory(u_vec_train_cell{1});
t = toc;
fprintf('Dauer der Linearisierung: %g s\n', t);

tic;
[P2, Var_P2] = GP_IO.linearize_at_given_trajectory_fast(u_vec_train_cell{1});
t = toc;
fprintf('Dauer der Linearisierung und Varianzberechnung (fast): %g s\n\n', t);

tic;
[P3, Var_P3] = GP_IO.approx_linearisation_at_given_trajectory(u_vec_train_cell{1});
t = toc;
fprintf('Dauer der Linearisierung (approximiert): %g s\n\n', t);

%% Compare Errors
% Compare Fast and slow Computation
error_P1 = P - P2;
max_error_P1 = max(abs(P(:) - P2(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von P (fast vs. slow): %.3e\n', max_error_P1);

% Compare analytic and approx Computation
error_P3 = P - P3;
max_error_P2 = max(abs(P(:) - P3(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von P (analytic vs. approx): %.3e\n', max_error_P2);

% Compare analytic and approx Computation
error_Var_P = Var_P2 - Var_P3;
max_error_Var_P = max(abs(Var_P2(:) - Var_P3(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von Var_P (analytic vs. approx): %.3e\n\n', max_error_Var_P);

% Prediction with linearized gp model
delta_u = u_vec_test - u_vec_train_cell{1};
y_pred_lin_test = GP_IO.predict_trajectory(u_vec_train_cell{1}) + P*delta_u;

%% Finite difference check
% Parameters
u0 = u_vec_train_cell{1};
K = 3;
h = 1e-6;

% First Prediction
[y0,~] = GP_IO.predict_trajectory(u0);
N = length(u0);

% Run Loop (one Direction per iteration)
for i = 1:K
    % Random Direction
    v = randn(N, 1);

    % Normalize direction vector (vector norm increases with vector size)
    v = v / norm(v);

    % Calculate change in output via Finite-Differences (from both sides)
    y_plus = GP_IO.predict_trajectory(u0 + h*v);
    y_minus = GP_IO.predict_trajectory(u0 - h*v);
    delta_y_fd = 0.5*(y_plus - y_minus);

    % Calculate Change in output via Linearized lifted Matrix
    delta_y_lin = P*(h*v);

    % Error
    max_error_y = max(abs(delta_y_fd - delta_y_lin));
    rel_error_y = norm(delta_y_fd - delta_y_lin) / max(1e-17, norm(delta_y_lin));

    % Print
    fprintf('Richtung %d: rel. Fehler = %.3e, max. abs. Fehler = %.3e\n', i, rel_error_y, max_error_y);
end
fprintf('\n');

%% Compare Lifted System representation
P_analytic = Lifted_dynamics_nonlinear_SISO(@(x) oszillator_linearized_discrete(x, Ts), N, m_delay, x_sim_train_cell{1});

% Predictions with analytic linearisation
y_lin_test = y_sim_train_cell{1} + [0; P_analytic*delta_u(1:end-1)];

% Compare lifted Matrices
P_reduced_size = P(2:end, 1:end-1);
error_P2 = P_analytic - P_reduced_size;
max_error_P2 = max(abs(P_analytic(:) - P_reduced_size(:)));
mean_error_P2 = mean(abs(P_analytic(:) - P_reduced_size(:)));

fprintf('Maximaler absoluter Unterschied P_analytic und P_lin: %.3e\n', max_error_P2);
fprintf('Mittlerer absoluter Unterschied P_analytic und P_lin: %.3e\n', mean_error_P2);

% Calculate prediction Error
error_y_pred = abs(y_sim_test - y_pred_test);
error_y_lin = abs(y_lin_test - y_pred_lin_test);

% Mean row-error P
row_err_P = zeros(size(error_P2, 1), 1);
for i = 1:size(error_P2, 1)
    row_err_P(i) = mean(abs(error_P2(i, 1:i)));
end

%% Plots
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(error_P2);
colorbar;
title('error P_{analytic} - P_{lin}');
xlabel('Input-Index');
ylabel('Output-Index');

subplot(1,2,2);
plot(t_vec, error_y_pred, LineWidth=2, DisplayName='error (y-sim - y-pred)');
hold on;
plot(t_vec, error_y_lin, LineWidth=2, DisplayName='error (y-lin - y-pred-lin)');
plot(t_vec(2:end), row_err_P, LineWidth=2, DisplayName='row-error (P-analyt - P-lin)');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Prediction Error');
legend()

figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test, LineWidth=1, DisplayName='y-pred-test');
plot(t_vec, y_pred_lin_test, LineWidth=1, DisplayName='y-pred-lin-test');
plot(t_vec, y_lin_test, LineWidth=1, DisplayName='y-lin-test');
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