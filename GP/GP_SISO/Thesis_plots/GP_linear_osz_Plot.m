% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      20.02.2026
% Beschreibung:
% In diesem skript werde ich Plots für die Masterarbeit erstellen, die
% zeigen, inwiefern der GP in der lage ist, die Systemmodelle zu
% approximieren.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
GP_Path = fullfile(base_dir, '..');
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
Plot_Path = fullfile(base_dir, '..', '..', '..', 'Plot');
addpath(GP_Path);
addpath(Model_Path);
addpath(Plot_Path);

%% Geeral
save_pdf = false;

%% Input Trajectorys
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

% Input trajectorys (Training)
u_scale_train = [1, 3];
N_traj = length(u_scale_train);

u_vec_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    u_vec_train_cell{i} = u_scale_train(i)*sin(2*pi/T_end.*t_vec');
end

% Input trajectorys (Test)
u_scale_test = 2;

u_vec_test1 = u_scale_test*sin(2*pi/T_end.*t_vec');
% u_vec_test = u_scale_test * (t_vec' >= 1) .* (t_vec' <= 3);

%% System Simulation
% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Initial condition (has to be constant)
x0 = [0;
    0];

% Simulation train

y_sim_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    [~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_cell{i} = x_sim(:, 1);
end

% Simulation test
[~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_test1, t_vec), t_vec, x0, opts);
y_sim_test1 = x_sim(:, 1);

%% Plot training Data
% Assign values (args)
args = struct();

args.x_cell = {t_vec, t_vec};
args.y_cell = {y_sim_train_cell, u_vec_train_cell};
args.x_label_cell = {'$t$', '$t$'};
args.y_label_cell = {'$y_L$', '$u_L$'};
args.title_cell = {'\textbf{a)}', '\textbf{b)}'};
args.legend_cell = {{'$y_{L,1}$', '$y_{L,2}$'}, {'$u_{L,1}$', '$u_{L,2}$'}};

args.filename = fullfile('04_Simulation_Experiment', 'Trainingsdaten_GP_IO.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 7;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_lim = {[], []};
opts.x_lim = {[], []};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
orientation = [1, 2];
plot.tiled_plot(opts, orientation);

%% Predict System Dynamics (with one training trajectory)
% Init Gaussian Process
GP_IO = GP_SISO_IO();

% Train Gaussian Process
GP_IO.train_GP_model(y_sim_train_cell, {u_vec_train_cell{1}});

% Predict new Trajectory with full covriance matrix (without measurement noise)
[y_pred_test, Cov_y_test] = GP_IO.predict_trajectory_covariance(u_vec_test1);
y_pred_test_upper = y_pred_test + 3*sqrt(diag(Cov_y_test));
y_pred_test_lower = y_pred_test - 3*sqrt(diag(Cov_y_test));

%% Plot Predictions
y_cell = {{y_pred_test, y_sim_test1}, {u_vec_test1}};
filename = fullfile('05_Ergebnisse_Diskussion', 'GP_Praediktion_linear_1u.pdf');
plot_tiled_GP_prediction(t_vec(:), y_cell, y_pred_test_upper, y_pred_test_lower, filename, save_pdf)

%% Local Functions
function plot_tiled_GP_prediction(x_plot, y_cell, y_upper, y_lower, filename, save_pdf)
    % Assign values (args)
    args = struct();
    
    args.x_cell = {x_plot, x_plot};
    args.y_cell = y_cell;
    args.x_label_cell = {'$t$', '$t$'};
    args.y_label_cell = {'$y_L$', '$u_L$'};
    args.title_cell = {'$\ell = 0.3$', '$\ell = 6$'};
    args.legend_cell = {{'3$\sigma$-Band', '$\hat{y}_{L*}$', '$y_{L*}$'}, {'$u_{L*}$'}};
    
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 7;
    opts.linewidth = 1;
    opts.y_scale = 'linear';
    opts.y_lim = {[], []};
    opts.x_lim = {[], []};
    opts.marker = 'none';
    
    % Uncertainty
    y_upper_cell = {y_upper, []};
    y_lower_cell = {y_lower, []};
    x_fill_cell = {x_plot, []};

    % Training data
    x_train_cell = cell(2, 1);
    y_train_cell = cell(2, 1);

    for i = 1:2
        x_train_cell{i} = [];
        y_train_cell{i} = [];
    end

    % Create Plot
    orientation = [1, 2];
    tiled_regression_plot = Plot_Manager(args);
    tiled_regression_plot.tiled_scatter_fill_plot(opts, orientation, x_train_cell, y_train_cell, x_fill_cell, y_upper_cell, y_lower_cell);
end