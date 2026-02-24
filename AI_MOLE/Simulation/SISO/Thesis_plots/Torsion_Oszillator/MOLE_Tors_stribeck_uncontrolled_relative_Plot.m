% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      23.02.2025
% Beschreibung:
% In diesem skript werde ich einen AI-MOLE Plot für die Masterarbeit
% erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', '..', '..', '..', 'ILC', 'Simulation', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'System_Models');
MOLE_Path = fullfile(base_dir, '..', '..');
Plot_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'Plot');
GP_path = fullfile(base_dir, '..', '..', '..', '..', '..', 'GP', 'GP_SISO');
addpath(ILC_path);
addpath(Model_Path);
addpath(MOLE_Path);
addpath(Plot_Path);
addpath(GP_path);

%% General
save_pdf = false;

%% Check if .mat-File exists
% Geth Name of current script
[~, script_name, ~] = fileparts(mfilename('fullpath'));

% Name of Data file
data_name = [script_name '.mat'];

% Evaluate AI-MOLE only when no data file exists
if ~isfile(data_name)

    %% Reference Trajectory
    % Load File
    filename = 'Trajectory_02.mat';
    filepath = fullfile(base_dir, '..', '..', '..', '..', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
    load(filepath);
    
    r_vec = ref_traj.phi2;
    r_vec = r_vec(:);

    % Simulation parameters
    Ts = ref_traj.t_vec(2) - ref_traj.t_vec(1);
    T_end = ref_traj.t_vec(end);
    t_vec = ref_traj.t_vec;
    
    % Solver settings
    opts = odeset( ...
        'RelTol', 1e-6, ...         % Tolerance
        'AbsTol', 1e-8, ...  % Tolerance
        'MaxStep', Ts/5, ...        % Use smaller step size for better Results
        'InitialStep', Ts/20);
    
    %% System Model
    % Choose system model
    dynamic_model = @torsion_oszillator_linear_stribeck;
    x0 = [0; 0; 0; 0];
    
    % Initial input Trajectory
    sigma_I = 8;
    u_init = sigma_I*sin(2*pi/T_end.*t_vec');
    
    %% Initialize AI-MOLE
    
    params = struct();
    
    % Parameters
    params.m_delay = 1;
    params.N_iter = 10;
    params.H_trials = 3;
    
    % Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
    % 'Robust', 'Manual')
    params.weight_init_method = 'Stochastic';
    
    % Choose nonlinearity damping method ('none', 'relative_1', 'relative_2', 'minimize')
    params.nonlin_damping = 'relative_2';
    params.beta = 0;
    
    % Initialisation
    SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);
    
    %% Run ILC
    tic;
    % Update Loop
    u_sim = u_init;
    
    [t_sim, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1);
    for i = 1:params.N_iter
        % Update input
        u_sim = [SISO_MOLE.update_input(y_sim); 0];
    
        % Simulate the system
        [t_sim, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_sim, t_vec), t_vec, x0, opts);
        y_sim = x_sim(:, 1);
    end
    SISO_MOLE.save_final_trajectory(y_sim);
    y_sim_quadratic = y_sim;
    
    % Zeitmessung und Ausgabe
    time = toc;
    fprintf('Dauer von AI-MOLE mit %d Iterationen, %d Trials Delay (H) und jeweils %d Datenpunkten pro Trial: %g s\n', params.N_iter, params.H_trials, length(t_vec), time);

    %% Save AI-MOLE Data
    % Delete large Matrices
    % GP
    SISO_MOLE.GP_SISO.L_chol = [];
    SISO_MOLE.GP_SISO.V_transp = [];
    SISO_MOLE.GP_SISO.V = [];
    SISO_MOLE.GP_SISO.GP = [];

    % ILC
    SISO_MOLE.ILC_SISO.W = [];
    SISO_MOLE.ILC_SISO.S = [];
    SISO_MOLE.ILC_SISO.R = [];
    SISO_MOLE.ILC_SISO.L = [];
    SISO_MOLE.ILC_SISO.Q = [];
    SISO_MOLE.ILC_SISO.P = [];

    % Save Results
    save(data_name, 'SISO_MOLE', 't_vec');

else
    %% Load results
    load(data_name);

end

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {};
args.y_cell = {};
args.x_label_cell = {'', '', '$t$', 'Iteration'};
args.y_label_cell = {'$y$', 'RMSE', '$u$', '$\eta$'};
args.title_cell = {'', '', '', ''};
args.legend_cell = {{'$y_d$', '$y_0$', '$y_5$', '$y_{10}$'}, {}, {'$u_0$', '$u_5$', '$u_{10}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Torsion_Stribeck_controlled_ff_relative.pdf');
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
plot.tiled_mole_results_plot(opts, SISO_MOLE, t_vec);