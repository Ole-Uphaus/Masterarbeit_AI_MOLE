% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.11.2025
% Beschreibung:
% In diesem skript werde ich AI-MOLE für den nichtlinearen Oszillator
% erproben.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

%% Reference Trajectory
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);

%% Initialize AI-MOLE
% Parameters
m_delay = 1;
N_iter = 10;
x0 = [0;
    0];

% Initialisation
SISO_MOLE = SISO_MOLE_IO(r_vec, m_delay);

%% Plot Results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
% plot(t_vec, y_sim_quadratic, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
% plot(1:length(ILC_Quadr.RMSE_log), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
% plot(1:length(ILC_PD.RMSE_log), ILC_PD.RMSE_log, LineWidth=1, DisplayName='ILC PD');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

%% Local Functions