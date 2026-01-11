% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      10.01.2026
% Beschreibung:
% Dieses Skript dient dazu, die Aktualisierte Eingangstrajektorie zu
% Nutzen, und das System mit dieser Eingangstrajektorie zu simulieren. Am
% realen Prüfstand wird anstelle dieses Skriptes ein Versuchslauf
% durchgeführt.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(Model_Path);

%% Controller design (discrete state feedback DR)
[sys_contin, sys_disc ,sys_disc_cont, Ts, S_gain, k_T_disc] = controller_design_torsion_oszillator();

%% Choose run for Simulation
% Run
date_string = '2026_01_11';
run_filename = 'Run_01_serial.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

% Extract architecture
name = erase(run_filename, '.mat');
parts = split(name, '_');
architecture = parts{3};

% Load file
load(run_filepath)

%% System simulation
% Extract current input trajectory (find last nonempty cell)
idx = find(~cellfun('isempty', ILC_Quadr.u_cell), 1, 'last');
u_vec = ILC_Quadr.u_cell{idx};

switch architecture
    case 'uncontrolled'
        % Use the uncontrolled System

        % Choose system model
        system_dynamics = sys_disc;
        x0 = [0; 0; 0; 0];

    case 'serial'
        % Use the PI-controlled system for AI-MOLE while modifying the
        % reference trajectory of the controller. The reference trajectory
        % is used as initial input.

        % Choose system model
        system_dynamics = sys_disc_cont;
        x0 = [0; 0; 0; 0];

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

% Upsample input trajectory
t_vec_sys = 0:Ts:ref_traj.t_vec(end);
t_vec_u = ref_traj.t_vec;

u_vec_sys = interp1(t_vec_u, u_vec, t_vec_sys, 'previous', 'extrap');

% Simulate system
[~, ~, x_sim] = lsim(system_dynamics, u_vec_sys(:), t_vec_sys(:), x0);
y_vec_sys = x_sim(:, 3);

% Downsample simulation results for ILC
Ts_u = t_vec_u(2) - t_vec_u(1);
sample_factor = round(Ts_u / Ts);

idx_down = 1:sample_factor:length(t_vec_sys);
y_vec = y_vec_sys(idx_down);

%% Save Results
% Name
res_name = sprintf('Trial_%s.mat', date_string);

% Save
sim_trial_timestamp = datetime('now');  % Timestamp of simulation
save(res_name, 'y_vec', 'sim_trial_timestamp');