% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      17.12.2025
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

%% Choose run for Simulation
% Run
date_string = '2025_12_29';
run_filename = 'Run_01_uncontrolled.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

% Extract architecture
name = erase(run_filename, '.mat');
parts = split(name, '_');
architecture = parts{3};

% Load file
load(run_filepath)

%% System simulation
% Extract current input trajectory (find last nonempty cell)
idx = find(~cellfun('isempty', SISO_MOLE.u_cell), 1, 'last');
u_vec = SISO_MOLE.u_cell{idx};

switch architecture
    case 'uncontrolled'
        % Use the uncontrolled System

        % Choose system model
        system_dynamics = @torsion_oszillator_linear;
        x0 = [0; 0; 0; 0];

    case 'serial'
        % Use the PI-controlled system for AI-MOLE while modifying the
        % reference trajectory of the controller. The reference trajectory
        % is used as initial input.

        % Choose system model
        system_dynamics = @torsion_oszillator_linear_PI;
        x0 = [0; 0; 0; 0; 0];

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

% Solver setting
Ts = ref_traj.t_vec(2) - ref_traj.t_vec(1);
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', 1e-8, ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation step
[t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_vec, ref_traj.t_vec), ref_traj.t_vec, x0, opts);
y_vec = x_sim(:, 3);

%% Save Results
% Name
res_name = sprintf('Sim_Run_%s_%s.mat', parts{2}, date_string);

% Save
sim_trial_timestamp = datetime('now');  % Timestamp of simulation
save(res_name, 'y_vec', 'sim_trial_timestamp');