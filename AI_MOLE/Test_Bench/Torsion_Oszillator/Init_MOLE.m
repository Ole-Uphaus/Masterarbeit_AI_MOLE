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
addpath(MOLE_Path);

%% Parameters
% General ILC architecture ('uncontrolled', 'serial', 'parallel')
architecture = 'uncontrolled';

% Choose reference trajectory
traj_name = 'Trajectory_01.mat';

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
        % Use the uncontrolled System and initialize with simple input
        % trajectory

        params = struct();

        % Parameters
        params.m_delay = 1;
        params.N_iter = 10;
        params.H_trials = 3;
        
        % Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
        % 'Robust', 'Manual')
        params.weight_init_method = 'Stochastic';
        
        % Choose nonlinearity damping Parameters
        params.use_nonlin_damping = true;
        params.beta = 2;
        
        % Initial input Trajectory (simple sin or automatic generated)
        sigma_I = 0.1;
        u_init_sin = sigma_I*sin(2*pi/ref_traj.t_vec(end).*ref_traj.t_vec');
        u_init = u_init_sin;        % u_init_sin / u_init_auto
        
        % Initialisation
        SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

    case 'serial'
        % Use the PI-controlled system for AI-MOLE while modifying the
        % reference trajectory of the controller. The reference trajectory
        % is used as initial input.

        params = struct();

        % Parameters
        params.m_delay = 1;
        params.N_iter = 10;
        params.H_trials = 3;
        
        % Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
        % 'Robust', 'Manual')
        params.weight_init_method = 'Stochastic';
        
        % Choose nonlinearity damping Parameters
        params.use_nonlin_damping = true;
        params.beta = 2;
        
        % Initial input Trajectory
        u_init = r_vec;
        
        % Initialisation
        SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

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