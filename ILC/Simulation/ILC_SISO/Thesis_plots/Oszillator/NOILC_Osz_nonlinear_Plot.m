% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      25.02.2026
% Beschreibung:
% In diesem skript werde ich einen NO-ILC Plot für die Masterarbeit
% erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'System_Models');
ILC_Path = fullfile(base_dir, '..', '..');
Plot_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'Plot');
addpath(Model_Path);
addpath(ILC_Path);
addpath(Plot_Path);

%% General
save_pdf = false;

%% Reference Trajectory
% Parameters
x_max = 0.5;
T_end = 5;
Ts = 0.01;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
u_init = zeros(size(r_vec, 1), 1);

%% ILC quadratic optimal design
% Lifted system dynamics
N = size(t_vec, 2);
m_delay = 1;

% Parameters
N_iter = 10;
x0 = [0;
    0]; 
W = eye(N - m_delay);
S = 0.00007*eye(N - m_delay);
R = 0*eye(N - m_delay);

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, (u_sim), t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    P = Lifted_dynamics_nonlinear_SISO(@(x) oszillator_linearized_discrete(x, Ts), N, m_delay, x_sim);
    ILC_Quadr.init_Quadr_type(W, S, R, P);
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, (u_sim), t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1);
end
% Calculate and log final error
ILC_Quadr.calculate_final_error(y_sim);

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {};
args.y_cell = {};
args.x_label_cell = {'', '', '$t$', 'Iteration'};
args.y_label_cell = {'$y$', 'RMSE', '$u$', '$\eta$'};
args.title_cell = {'', '', '', ''};
args.legend_cell = {{'$y_d$', '$y_0$', '$y_5$', '$y_{10}$'}, {}, {'$u_0$', '$u_5$', '$u_{10}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_nonlinear_Modell.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_lim = {[], [], [], []};
opts.x_lim = {[], [], [], [0, 10]};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
plot.tiled_ilc_results_plot(opts, ILC_Quadr, t_vec);