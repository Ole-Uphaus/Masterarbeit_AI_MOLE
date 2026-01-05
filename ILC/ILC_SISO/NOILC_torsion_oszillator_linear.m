% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      05.01.2025
% Beschreibung:
% In diesem skript werde ich die ILC Regelung für den linearen
% Torsionsschwinger wie am Prüfstand implementieren. ich werde dazu das
% lineare Systemmodell nutzen und auch das geregelte system betrachten.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', 'System_Models');
addpath(Model_Path);

%% Reference Trajectory
filename = 'Trajectory_01.mat';
filepath = fullfile(pwd, '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
load(filepath);

r_vec = ref_traj.phi2;
r_vec = r_vec(:);

%% System dynamics
% Simulation parameters
J1  = 0.031;    % kgm^2
J2  = 0.237;    % kgm^2
c_phi = 9;      % Nm/rad
d_v1 = 0.070;   % Nms/rad
d_v2 = 0.231;   % Nms/rad

% State space
A = [0, 1, 0, 0;
    -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
    0, 0, 0, 1;
    c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
B = [0;
    1/J1;
    0;
    0];
C = [0, 0, 1, 0];
D = 0;

%% Parameters
% ILC architcture ('uncontrolled', 'serial_LQR', 'serial_LQR_ff')
architecture = 'uncontrolled';

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

% ILC parameter
N = size(t_vec, 2);
m_delay = 1;

%% ILC Initialisation
switch architecture
    case 'uncontrolled'
        % Use the uncontrolled System for model based ILC

        % System dynamics
        sys_cont = ss(A,B,C,D);
        sys_disc = c2d(sys_cont, Ts, 'zoh');
        [Ad,Bd,Cd,Dd] = ssdata(sys_disc);

        % Choose system model
        system_dynamics = @torsion_oszillator_linear;
        x0 = [0; 0; 0; 0];

        % Lifted system dynamics
        P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        S = 0.0001*eye(size(P));
        R = 0*eye(size(P));

        % Initial Input Trajectory
        u_init = zeros(size(r_vec, 1), 1);
        
        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

%% Run ILC
% Track Results
y_cell_quadr = cell(N_iter+1, 1);
u_cell_quadr = cell(N_iter+1, 1);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 3);
y_cell_quadr{1} = y_sim;
u_cell_quadr{1} = u_sim;
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 3);
    y_cell_quadr{i+1} = y_sim;
    u_cell_quadr{i+1} = u_sim;
end
% Calculate and log final error
ILC_Quadr.calculate_final_error(y_sim);
