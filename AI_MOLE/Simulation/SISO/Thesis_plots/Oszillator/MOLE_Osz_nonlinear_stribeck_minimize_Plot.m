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
addpath(ILC_path);
addpath(Model_Path);
addpath(MOLE_Path);
addpath(Plot_Path);

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
    % Parameters
    x_max = 0.5;
    Ts = 0.01;
    T_end = 5;
    x0 = [0;
        0];
    
    t_vec = 0:Ts:T_end;
    
    % Solver settings
    opts = odeset( ...
        'RelTol', 1e-6, ...         % Tolerance
        'AbsTol', [1e-8 1e-8], ...  % Tolerance
        'MaxStep', Ts/5, ...        % Use smaller step size for better Results
        'InitialStep', Ts/20);
    
    % Trajectory (no delay - delay is applied later)
    sigma = 1;
    [r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
    
    %% System Model
    % Use the linear system for AI-MOLE
    dynamic_model = @oszillator_nonlinear_stribeck;
    
    % Initial input Trajectory (simple sin or automatic generated)
    sigma_I = 1;  % for stibeck model = 1, otherwise = 0.1
    u_init = sigma_I*sin(2*pi/T_end.*t_vec');
    
    %% Initialize AI-MOLE
    
    params = struct();
    
    % Parameters
    params.m_delay = 1;
    params.N_iter = 20;
    params.H_trials = 3;
    
    % Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
    % 'Robust', 'Manual')
    params.weight_init_method = 'Stochastic';
    
    % Choose nonlinearity damping method ('none', 'relative_1', 'relative_2', 'minimize')
    params.nonlin_damping = 'minimize';
    params.beta = 0.5;
    
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
args.legend_cell = {{'$y_d$', '$y_0$', '$y_{10}$', '$y_{20}$'}, {}, {'$u_0$', '$u_{10}$', '$u_{20}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_nonlinear_stribeck_minimize.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_lim = {[], [], [], []};
opts.x_lim = {[], [], [], [0, 20]};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
plot.tiled_mole_results_plot(opts, SISO_MOLE, t_vec);