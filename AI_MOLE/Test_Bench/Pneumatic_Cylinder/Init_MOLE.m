% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      03.03.2026
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
[S_gain] = controller_design_pneumatic_cylinder();

%% Parameters
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
% Use the state feedback controlled (LQR) system for AI-MOLE while
% modifying the reference trajectory if the static feedforward
% controlled system.

params = struct();

% Parameters
params.m_delay = 1;
params.N_iter = 10;
params.H_trials = 10;

% Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
% 'Robust', 'Manual')
params.weight_init_method = 'Stochastic';

% Choose nonlinearity damping method ('none', 'relative_1', 'relative_2', 'minimize')
params.nonlin_damping = 'relative_2';
params.beta = 2;

% Initial input Trajectory (simple sin or static feed forward)
u_init = S_gain .* r_vec;

% Initialisation
SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

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
run_filename = sprintf('Run_%02d.mat', run_idx);
run_filepath = fullfile(date_dir, run_filename);

% Save data
init_update_timestamp = datetime('now');    % Timestamp of initialisation
save(run_filepath, 'SISO_MOLE', 'ref_traj', 'init_update_timestamp');