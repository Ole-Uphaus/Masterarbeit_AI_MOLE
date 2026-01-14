% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      17.12.2025
% Beschreibung:
% Dieses Skript dient dazu, einen einzelnen AI-MOLE Iterationsschritt
% am Prüfstand durchzuführen. Es werden dazu die aktuellen Messwerte
% ausgelesen und anschließend die optimierte Eingangstrajektorie
% zurückgegeben.
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

%% Load MOLE and simulation/trial 
% MOLE object
date_string = '2026_01_14';
run_filename = 'Run_05_serial.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

% Current Simulation/Trial
sim_trial_filename = sprintf('Trial_%s.mat', date_string);

% Extract architecture
name = erase(run_filename, '.mat');
parts = split(name, '_');
architecture = parts{3};

% Load files
load(run_filepath);
load(sim_trial_filename);

% Variable for saving results (only true if AI-MOLE does an update)
save_results = false;

%% MOLE Update
% Check if timestamps are in the right order (this prevents from performing
% an update before performing a new simulation/trial)
if init_update_timestamp < sim_trial_timestamp
    % Chech if this was the last iteration (no update will be performed - just
    % the last trajectory will be saved)
    idx_u = find(~cellfun('isempty', SISO_MOLE.u_cell), 1, 'last');
    idx_y = find(~cellfun('isempty', SISO_MOLE.y_cell), 1, 'last');
    if idx_u <= SISO_MOLE.N_iter
        % Perform MOLE update
        [~] = SISO_MOLE.update_input(y_vec);
        disp('AI-MOLE Update durchgeführt.')

        save_results = true;
    elseif (idx_u == SISO_MOLE.N_iter+1) && (idx_y == SISO_MOLE.N_iter)
        % Save last trajectory
        SISO_MOLE.save_final_trajectory(y_vec);
        disp('Letzte trajektorie gespeichert.')

        save_results = true;
    else
        disp('Es wurde kein Update durchgeführt, da die maximale Anzahl an Iterationen erreicht wurde.')
    end
else
    disp('Es wurde kein Update durchgeführt, da zuerst ein neuer Trial/Simulation benötigt wird.')
end

%% Save results
if save_results
    % Overwrite existing .mat-file
    init_update_timestamp = datetime('now');    % Timestamp of update
    save(run_filepath, 'SISO_MOLE', 'ref_traj', 'init_update_timestamp');
end

%% Calculate approximate actuator input
% Get the latest u input
idx_u = find(~cellfun('isempty', SISO_MOLE.u_cell), 1, 'last');

switch architecture
    case 'serial'
        % Controller design
        [sys_contin, sys_disc ,sys_disc_cont, Ts, S_gain, k_T_disc] = controller_design_torsion_oszillator();

        % Upsample input trajectory
        t_vec_sys = 0:Ts:ref_traj.t_vec(end);
        t_vec_u = ref_traj.t_vec;
        
        u_vec_sys = interp1(t_vec_u, SISO_MOLE.u_cell{idx_u}, t_vec_sys, 'previous', 'extrap');

        % Simulate system
        x0 = [0; 0; 0; 0];
        [~, ~, x_sim] = lsim(sys_disc_cont, u_vec_sys(:), t_vec_sys(:), x0);

        % Calculate actuator input
        u_vec_actuator = u_vec_sys(:) - x_sim*k_T_disc.';

    case 'uncontrolled'
        % Controller design
        [sys_contin, sys_disc ,sys_disc_cont, Ts, S_gain, k_T_disc] = controller_design_torsion_oszillator();

        % Upsample input trajectory
        t_vec_sys = 0:Ts:ref_traj.t_vec(end);
        t_vec_u = ref_traj.t_vec;
        
        u_vec_actuator = interp1(t_vec_u, SISO_MOLE.u_cell{idx_u}, t_vec_sys, 'previous', 'extrap');
    otherwise
        error('Architektur nicht erkannt.')
end

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(ref_traj.t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:SISO_MOLE.N_iter
    if ~isempty(SISO_MOLE.y_cell{i})
        plot(ref_traj.t_vec, SISO_MOLE.y_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
    end
end
if ~isempty(SISO_MOLE.y_cell{SISO_MOLE.N_iter+1})
    plot(ref_traj.t_vec, SISO_MOLE.y_cell{SISO_MOLE.N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', SISO_MOLE.N_iter));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Compare desired and simulated Trajectory');
legend('Location', 'best');

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
% plot(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
semilogy(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,4);   % 1 Zeile, 2 Spalten, erster Plot
hold on;
for i = 1:SISO_MOLE.N_iter
    if ~isempty(SISO_MOLE.u_cell{i})
        plot(ref_traj.t_vec, SISO_MOLE.u_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
    end
end
if ~isempty(SISO_MOLE.u_cell{SISO_MOLE.N_iter+1})
    plot(ref_traj.t_vec, SISO_MOLE.u_cell{SISO_MOLE.N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', SISO_MOLE.N_iter));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('T [Nm]');
title('Input Signal');
legend('Location', 'best');

subplot(2,2,2);
% plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Noise');
legend()

figure;
plot(t_vec_sys, u_vec_actuator, LineWidth=1); hold on;
yline(9,  'r--', 'LineWidth', 1);
yline(-9, 'r--', 'LineWidth', 1);
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Approximate Actuator Input Signal');
ylim([-10 10]);