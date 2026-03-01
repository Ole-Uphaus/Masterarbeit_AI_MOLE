% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      28.02.2025
% Beschreibung:
% In diesem skript werde ich eine Rechenzeitanalyse von einer Iteration von
% AI-MOLE durchführen. Dabei werden verschiedene Anzahlen an Trainingsdaten
% genutzt.
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

    %% Reference Trajectorys with different sample times Ts
    % Parameters
    T_end = 5;
    x0 = [0;
        0];

    % Different sample times
    max_samples = 1000;
    samples_vec = linspace(max_samples/10, max_samples, 10);

    % Create reference Trajectory
    ref_traj_cell = cell(length(samples_vec), 1);

    for i = 1:length(samples_vec)
        % Time vector
        t_vec = linspace(0, T_end, samples_vec(i));

        % Trajectory (no delay - delay is applied later)
        sigma = 1;
        rng(43);    % Every Trajectory is the same
        [r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);

        % Save in Cell
        ref_traj_cell{i}.t_vec = t_vec;
        ref_traj_cell{i}.r_vec = r_vec;
    end
    
    %% System Model
    % Use the linear system for AI-MOLE
    dynamic_model = @oszillator_linear;
    
    %% AI-MOLE Parameters
    
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
    params.beta = 0.5;

    %% Runtime measurements
    % Empty results vector
    runtime_vec = zeros(length(samples_vec), 1);

    % Calculate runtime for first iteration each
    for i = 1:length(samples_vec)
        % Solver settings
        Ts = ref_traj_cell{i}.t_vec(2) - ref_traj_cell{i}.t_vec(1);
        opts = odeset( ...
            'RelTol', 1e-6, ...         % Tolerance
            'AbsTol', [1e-8 1e-8], ...  % Tolerance
            'MaxStep', Ts/5, ...        % Use smaller step size for better Results
            'InitialStep', Ts/20);

        % Initial input Trajectory
        sigma_I = 0.1;
        u_init = sigma_I*sin(2*pi/T_end.*ref_traj_cell{i}.t_vec');

        % System simulation
        u_sim = u_init;
        [~, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_sim, ref_traj_cell{i}.t_vec), ref_traj_cell{i}.t_vec, x0, opts);
        y_sim = x_sim(:, 1);

        % Initialize AI-MOLE
        SISO_MOLE = SISO_MOLE_IO(ref_traj_cell{i}.r_vec, u_init, params);

        % Calculate first AI-MOLE iteration
        tic;
        [~] = [SISO_MOLE.update_input(y_sim); 0];
        runtime_vec(i) = toc;

        % Delete MOLE object
        delete(SISO_MOLE);
        clear SISO_MOLE;
    end

    %% Save Runtime Data
    % Save Results
    save(data_name, 'samples_vec', 'runtime_vec');

else
    %% Load results
    load(data_name);

end

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {samples_vec};
args.y_cell = {{runtime_vec}};
args.x_label_cell = {'$N$'};
args.y_label_cell = {'$t$ in $\mathrm{s}$'};
args.title_cell = {''};
args.legend_cell = {{}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_MOLE_Rechenzeitanalyse.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_lim = {[]};
opts.x_lim = {[0, 1000]};
opts.marker = '.';

% Create Plot
plot = Plot_Manager(args);
Position = [0.27, 0.20, 0.55, 0.72];
plot.single_plot(opts, Position);