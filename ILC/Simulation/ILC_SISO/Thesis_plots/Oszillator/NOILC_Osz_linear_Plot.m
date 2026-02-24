% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      24.02.2026
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
addpath(Model_Path);

%% General
save_pdf = false;

%% System Dynamics
% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
d  = 0.5; % Ns/m
m_delay = 1;

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
S = 0.00007*eye(size(P));
R = 0*eye(size(P));

% Track Results
y_cell_quadr = cell(N_iter+1, 1);
u_cell_quadr = cell(N_iter+1, 1);

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
ILC_Quadr.init_Quadr_type(W, S, R, P)

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim), t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim), t_vec), t_vec, x0, opts);
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

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_linear_Modell.pdf');
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