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
data_1 = load(fullfile(pwd, 'Trajectory_4', 'Messungen_67', 'Messungen_gesamt.mat'));
data_2 = load(fullfile(pwd, 'Trajectory_4', 'Messungen_68', 'Messungen_gesamt.mat'));

%% Detect Start Time
% Reference starting Point
bool_reference_1 = data_1.Data_ges(:, 16);
bool_reference_2 = data_2.Data_ges(:, 16);

% Find first index with one
idx_start_1 = find(bool_reference_1 == 1, 1);
idx_start_2 = find(bool_reference_2 == 1, 1);

t_start_1 = data_1.time_ges(idx_start_1);
t_start_2 = data_2.time_ges(idx_start_2);

%% Separate relevant data part
% ILC Time 
delta_t = 0.4;    % Trajectory lasts 1s
t_end_1 = t_start_1 + delta_t;
t_end_2 = t_start_2 + delta_t;

% Cut out results vectors
idx_ilc_1 = (data_1.time_ges >= t_start_1 & data_1.time_ges <= t_end_1);
idx_ilc_2 = (data_2.time_ges >= t_start_2 & data_2.time_ges <= t_end_2);

t_vec_ILC_1 = data_1.time_ges(idx_ilc_1);
t_vec_ILC_2 = data_2.time_ges(idx_ilc_2);

t_vec_ILC_1 = t_vec_ILC_1 - t_vec_ILC_1(1);   % Adjust range
t_vec_ILC_2 = t_vec_ILC_2 - t_vec_ILC_2(1);

data_vec_ILC_1 = data_1.Data_ges(idx_ilc_1, :);
data_vec_ILC_2 = data_2.Data_ges(idx_ilc_2, :);

%% Calculate RMSE
% Output Vectors
y_vec_1 = data_vec_ILC_1(:, 1);
y_vec_2 = data_vec_ILC_2(:, 1);

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
plot(t_vec_ILC_1, y_vec_1, LineWidth=1, DisplayName='Original1'); hold on;
plot(t_vec_ILC_2, y_vec_2, LineWidth=1, DisplayName='Original2'); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('ILC Run');
legend('Location', 'best');