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
addpath(MOLE_Path);

%% Load MOLE and simulation/trial 
% MOLE object
date_string = '2026_01_08';
run_filename = 'Run_01_serial.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

% Current Simulation/Trial
sim_trial_filename = sprintf('Trial_%s.mat', date_string);

% Load files
load(run_filepath);
load(sim_trial_filename);

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
    elseif (idx_u == SISO_MOLE.N_iter+1) && (idx_y == SISO_MOLE.N_iter)
        % Save last trajectory
        SISO_MOLE.save_final_trajectory(y_vec);
        disp('Letzte trajektorie gespeichert.')
    else
        disp('Es wurde kein Update durchgeführt, da die maximale Anzahl an Iterationen erreicht wurde.')
    end
else
    disp('Es wurde kein Update durchgeführt, da zuerst ein neuer Trial/Simulation benötigt wird.')
end

%% Save results
% Overwrite existing .mat-file
init_update_timestamp = datetime('now');    % Timestamp of update
save(run_filepath, 'SISO_MOLE', 'ref_traj', 'init_update_timestamp');

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
ylabel('x [m]');
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
ylabel('F [N]');
title('Input Signal');
legend('Location', 'best');

subplot(2,2,2);
% plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Noise');
legend()