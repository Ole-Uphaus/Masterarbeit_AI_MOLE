% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      10.01.2026
% Beschreibung:
% Dieses Skript dient dazu, ein neuen ILC-Run zu initialisieren.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
MOLE_Path = fullfile(base_dir, '..', '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator');
ILC_Path = fullfile(base_dir, '..', '..', 'Simulation', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(ILC_Path);
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
traj_path = fullfile(MOLE_Path, 'Reference_Trajectories', traj_name);

% Load Trajectory Data
load(traj_path);
r_vec = ref_traj.phi2;
r_vec = r_vec(:);

% ILC parameter
N = size(ref_traj.t_vec, 2);
m_delay = 1;

%% Initialize ILC object
% ILC sample Time
Ts_ILC = ref_traj.t_vec(2) -ref_traj.t_vec(1);

switch architecture
    case 'uncontrolled'
        % Use the uncontrolled system and initialize with zero Input.

        % System dynamics (ILC sample time)
        sys_disc_ILC = c2d(sys_contin, Ts_ILC, 'zoh');
        [Ad_ILC,Bd_ILC,Cd_ILC,Dd_ILC] = ssdata(sys_disc_ILC);

        % Lifted system dynamics (uncontrolled)
        P = Lifted_dynamics_linear_SISO(Ad_ILC, Bd_ILC, Cd_ILC, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        % Simulation with s = 10 leads to expected errors
        S = 10*eye(size(P));
        R = 0*eye(size(P));

        % Initial Input Trajectory
        u_init = zeros(size(r_vec, 1), 1);
        
        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    case 'serial'
        % Use the state feedback controlled (LQR) system for AI-MOLE while
        % modifying the unput trajectory of the static feedforward
        % controlled system.

        % System dynamics (ILC sample time)
        sys_disc_ILC = c2d(sys_contin, Ts_ILC, 'zoh');
        [Ad_ILC,Bd_ILC,Cd_ILC,Dd_ILC] = ssdata(sys_disc_ILC);

        % Controlled System dynamics
        Ad_ILC_cont = Ad_ILC - Bd_ILC*k_T_disc;

        % Lifted system dynamics (controlled)
        P = Lifted_dynamics_linear_SISO(Ad_ILC_cont, Bd_ILC, Cd_ILC, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        % Simulation with s = 0.1 leads to expected errors
        S = 0.1*eye(size(P));
        R = 0*eye(size(P));
        
        % Initial input Trajectory (simple sin or static feed forward)
        use_feedforward_control = false;

        if use_feedforward_control
            u_init = S_gain .* r_vec;
        else
            u_init = zeros(size(r_vec, 1), 1);
        end

        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

%% Save ILC object
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
save(run_filepath, 'ILC_Quadr', 'ref_traj', 'init_update_timestamp');