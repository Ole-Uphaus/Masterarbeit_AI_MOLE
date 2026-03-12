% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      12.03.2026
% Beschreibung:
% Dieses Skript dient dazu, den Wiederholungsfehler von zwei
% Prüfstandsiterationen zu berechnen.
% -------------------------------------------------------------

clc
clear
close all

%% Merge and load measurement data
% Load measurement data
data_1 = load(fullfile(pwd, 'Trajectory_2', 'Messwerte_2', 'Messungen_gesamt.mat'));
data_2 = load(fullfile(pwd, 'Trajectory_2', 'Messwerte_3', 'Messungen_gesamt.mat'));

%% Separate relevant data part
% ILC Time
t_start = 5;
t_end = 10;

% Cut out results vectors
idx_ilc = (data_1.time_ges >= t_start & data_1.time_ges <= t_end);

t_vec_ILC = data_1.time_ges(idx_ilc);
t_vec_ILC = t_vec_ILC - t_vec_ILC(1);   % Adjust range

data_vec_ILC_1 = data_1.Data_ges(idx_ilc, :);
data_vec_ILC_2 = data_2.Data_ges(idx_ilc, :);

%% Downsample measurement data for ILC
% Demanded sample Time (has to fit the actual sample time)
Ts_ILC = 0.01;

% Compare to actual sample time
Ts = t_vec_ILC(2) - t_vec_ILC(1);
sample_factor = round(Ts_ILC / Ts);

% Downsample
idx_ilc_down = 1:sample_factor:length(t_vec_ILC);

t_vec_ILC_down = t_vec_ILC(idx_ilc_down);
data_vec_ILC_down_1 = data_vec_ILC_1(idx_ilc_down, :);
data_vec_ILC_down_2 = data_vec_ILC_2(idx_ilc_down, :);

%% Calculate RMSE
% Output Vectors
y_vec_1 = data_vec_ILC_down_1(:, 3);
y_vec_2 = data_vec_ILC_down_2(:, 3);

% Error Vector
error_vec = y_vec_1 - y_vec_2;
RMSE = sqrt(mean(error_vec.^2));

% Display
fprintf('RMSE: %e\n', RMSE);

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(data_1.time_ges, data_1.Data_ges(:, 3), LineWidth=1); hold on;
plot(data_2.time_ges, data_2.Data_ges(:, 3), LineWidth=1); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Full Run');
legend('Location', 'best');

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec_ILC, data_vec_ILC_1(:, 3), LineWidth=1, DisplayName='Original1'); hold on;
plot(t_vec_ILC, data_vec_ILC_2(:, 3), LineWidth=1, DisplayName='Original2'); hold on;
plot(t_vec_ILC_down, y_vec_1, LineWidth=1, DisplayName='Downsampled1'); hold on;
plot(t_vec_ILC_down, y_vec_2, LineWidth=1, DisplayName='Downsampled2'); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('ILC Run');
legend('Location', 'best');