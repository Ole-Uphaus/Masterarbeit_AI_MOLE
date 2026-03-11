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
    T_end = 5;
    Ts = 0.01;
    
    t_vec = 0:Ts:T_end;
    
    % Trajectory (no delay - delay is applied later)
    sigma = 1;
    [r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
    u_init = zeros(size(r_vec, 1), 1);
    
    % Noise Parameters
    sigma_v = 0.01;      % Measurement Noise 0.01
    fc_v = 20;
    white = true;       % if white == true -> white noise is sampled - no filter
    
    %% ILC quadratic optimal design
    % Lifted system dynamics
    N = size(t_vec, 2);
    m_delay = 1;
    
    % Parameters
    N_iter = 15;
    x0 = [0;
        0]; 
    W = eye(N - m_delay);
    S = 0.002*eye(N - m_delay);
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
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
    u_sim = [ILC_Quadr.u_vec; 0];
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, (u_sim), t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
    for i = 1:N_iter
        % Update input
        P = Lifted_dynamics_nonlinear_SISO(@(x) oszillator_linearized_discrete(x, Ts), N, m_delay, x_sim);
        ILC_Quadr.init_Quadr_type(W, S, R, P);
        u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];
    
        % Simulate the system
        v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
        [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, (u_sim), t_vec), t_vec, x0, opts);
        y_sim = x_sim(:, 1) + v_vec;
    end
    % Calculate and log final error
    ILC_Quadr.calculate_final_error(y_sim);

    %% Save AI-MOLE Data
    % Delete large Matrices
    ILC_Quadr.W = [];
    ILC_Quadr.S = [];
    ILC_Quadr.R = [];
    ILC_Quadr.L = [];
    ILC_Quadr.Q = [];
    ILC_Quadr.P = [];

    % Save Results
    save(data_name, 'ILC_Quadr', 't_vec');

else
    %% Load results
    load(data_name);

end

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {};
args.y_cell = {};
args.x_label_cell = {'', '', '$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$y_L$ in $\mathrm{m}$', 'RMSE in $\mathrm{m}$', '$u_L$ in $\mathrm{N}$', '$\eta$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}', '\textbf{(c)}', '\textbf{(d)}'};
args.legend_cell = {{'$y_{L,d}$', '$y_{L,0}$', '$y_{L,8}$', '$y_{L,15}$'}, {}, {'$u_{L,0}$', '$u_{L,8}$', '$u_{L,15}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_nonlinear_NOILC_noise.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = {'linear', 'log', 'linear', 'linear'};
opts.y_lim = {[], [], [], [0, 1]};
opts.x_lim = {[], [], [], [0, 15]};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
plot.tiled_ilc_results_plot(opts, ILC_Quadr, t_vec);