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
ILC_Path = fullfile(base_dir, '..', '..', 'Simulation', 'ILC_SISO');
MOLE_Path = fullfile(base_dir, '..', '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator');
GP_Path = fullfile(base_dir, '..', '..', '..', 'GP', 'GP_SISO');
addpath(MOLE_Path);
addpath(ILC_Path);
addpath(GP_Path);

%% Load MOLE object
% MOLE object
date_string = '2026_01_14';
run_filename = 'Run_01_serial.mat';
run_filepath = fullfile(pwd, 'Runs', date_string, run_filename);

load(run_filepath);

% Iteration to plot
plot_iter = 10;

%% Plot results
% Get number of iterations
N_iter = length(ILC_Quadr.u_cell) - 1;

figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(ref_traj.t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:N_iter
    if ~isempty(ILC_Quadr.y_cell{i})
        plot(ref_traj.t_vec, ILC_Quadr.y_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
    end
end
if ~isempty(ILC_Quadr.y_cell{N_iter+1})
    plot(ref_traj.t_vec, ILC_Quadr.y_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Compare desired and simulated Trajectory');
legend('Location', 'best');

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
% plot(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
semilogy(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,4);   % 1 Zeile, 2 Spalten, erster Plot
hold on;
for i = 1:N_iter
    if ~isempty(ILC_Quadr.u_cell{i})
        plot(ref_traj.t_vec, ILC_Quadr.u_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
    end
end
if ~isempty(ILC_Quadr.u_cell{N_iter+1})
    plot(ref_traj.t_vec, ILC_Quadr.u_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
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

% Plot 2

figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(ref_traj.t_vec, ref_traj.phi2, LineWidth=1, DisplayName='desired'); hold on;
if ~isempty(ILC_Quadr.y_cell{plot_iter+1})
    plot(ref_traj.t_vec, ILC_Quadr.y_cell{plot_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', plot_iter));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Compare desired and simulated Trajectory');
legend('Location', 'best');

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
% plot(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
semilogy(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
trial_vec = 0:(length(ILC_Quadr.RMSE_log)-1);
plot(trial_vec(plot_iter+1), ILC_Quadr.RMSE_log(plot_iter+1), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r')
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,4);   % 1 Zeile, 2 Spalten, erster Plot
hold on;
if ~isempty(ILC_Quadr.u_cell{plot_iter+1})
    plot(ref_traj.t_vec, ILC_Quadr.u_cell{plot_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', plot_iter));
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