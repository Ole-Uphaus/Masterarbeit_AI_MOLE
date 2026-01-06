% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      29.12.2025
% Beschreibung:
% Dieses Skript dient dazu, die Ergebnisse aus einem AI-MOLE durchlauf zu
% plotten.
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
MOLE_Path = fullfile(base_dir, '..', '..', 'Simulation', 'SISO');
ILC_Path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');
GP_Path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
addpath(MOLE_Path);
addpath(ILC_Path);
addpath(GP_Path);

%% Load MOLE object
% MOLE object
date_string = '2025_12_29';
run_filename = 'Run_01_uncontrolled.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

load(run_filepath);

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
plot(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
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