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
addpath(ILC_Path);

%% Controller design (discrete state feedback DR)
% Sample Time
Ts = 0.001;

% Simulation parameters
J1  = 0.0299;    % kgm^2
J2  = 0.0299;    % kgm^2
c_phi = 7.309;   % Nm/rad
d_v1 = 0.055;    % Nms/rad
d_v2 = 0.0064;   % Nms/rad

% State space
A = [0, 1, 0, 0;
    -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
    0, 0, 0, 1;
    c_phi/J2, 0, -c_phi/J2, -d_v2/J2];

b = [0;
    1/J1;
    0;
    0];

c_T = [0, 0, 1, 0];

d = 0;

% Discrete System
sys_contin = ss(A, b, c_T, 0);
sys_disc = c2d(sys_contin, Ts, 'zoh');

% System matrices
Ad = sys_disc.A;
bd = sys_disc.B;
c_Td = sys_disc.C;
dd = sys_disc.D;

% LQR weighting matrices (as in DR)
Q_LQR = diag([1, 1, 10, 1]);
R_LQR = 1;

% LQR gain
k_T_disc = dlqr(Ad, bd, Q_LQR, R_LQR);

% Controlled system dynamics
Ad_cont = Ad - bd * k_T_disc;
sys_disc_cont = ss(Ad_cont, bd, c_Td, dd, Ts);

% Static gain (feedforward control signal)
S_gain = 1 / dcgain(sys_disc_cont);

%% Parameters
% General ILC architecture ('uncontrolled', 'serial')
architecture = 'uncontrolled';

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
        S = 0.1*eye(size(P));
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
        S = 0.00000001*eye(size(P));
        R = 0*eye(size(P));
        
        % Initial input Trajectory (simple sin or static feed forward)
        use_feedforward_control = true;

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