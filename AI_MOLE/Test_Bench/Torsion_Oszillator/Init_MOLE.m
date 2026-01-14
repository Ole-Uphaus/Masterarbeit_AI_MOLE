% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      17.12.2025
% Beschreibung:
% Dieses Skript dient dazu, ein neuen AI-MOLE Run zu initialisieren.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
MOLE_Path = fullfile(base_dir, '..', '..', 'Simulation', 'SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(MOLE_Path);
addpath(Model_Path);

%% Controller design (discrete state feedback DR)
[sys_contin, sys_disc ,sys_disc_cont, Ts, S_gain, k_T_disc] = controller_design_torsion_oszillator();

%% Parameters
% General ILC architecture ('uncontrolled', 'serial')
architecture = 'serial';

% Choose reference trajectory
traj_name = 'Trajectory_02.mat';

%% Load reference Trajectory
% File path
traj_path = fullfile(pwd, 'Reference_Trajectories', traj_name);

% Load Trajectory Data
load(traj_path);
r_vec = ref_traj.phi2;
r_vec = r_vec(:);

%% Initialize AI-MOLE object
switch architecture
    case 'uncontrolled'
        % Use the uncontrolled system and initialize with simple input
        % trajectory

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
        
        % Initial input Trajectory (simple sin or automatic generated)
        sigma_I = 7.5;
        u_init_sin = sigma_I*sin(2*pi/ref_traj.t_vec(end).*ref_traj.t_vec');
        u_init = u_init_sin;        % u_init_sin / u_init_auto
        
        % Initialisation
        SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

    case 'serial'
        % Use the state feedback controlled (LQR) system for AI-MOLE while
        % modifying the reference trajectory if the static feedforward
        % controlled system.

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
        
        % Initial input Trajectory (simple sin or static feed forward)
        use_feedforward_control = true;

        if use_feedforward_control
            u_init = S_gain .* r_vec;
        else
            sigma_I = 20;
            u_init = sigma_I*sin(2*pi/ref_traj.t_vec(end).*ref_traj.t_vec');
        end

        % Initialisation
        SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params, 0.001);

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

%% Save MOLE object
% Runs directory
runs_dir = fullfile(pwd, 'Runs');

% Date directory (today)
today_str = datestr(datetime('today'), 'yyyy_mm_dd');
date_dir = fullfile(runs_dir, today_str);

if ~exist(date_dir, 'dir')
    mkdir(date_dir);
end

% Count runs in directory
run_files = dir(fullfile(date_dir, 'Run_*.mat'));
run_idx = numel(run_files) + 1;

% Filename
run_filename = sprintf('Run_%02d_%s.mat', run_idx, architecture);
run_filepath = fullfile(date_dir, run_filename);

% Save data
init_update_timestamp = datetime('now');    % Timestamp of initialisation
save(run_filepath, 'SISO_MOLE', 'ref_traj', 'init_update_timestamp');