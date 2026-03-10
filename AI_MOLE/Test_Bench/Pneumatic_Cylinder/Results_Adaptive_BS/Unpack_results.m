% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      10.03.2026
% Beschreibung:
% Dieses Skript dient dazu die Regelungsergebnisse vom Adaptive
% Backstepping zu entpacken.
% -------------------------------------------------------------

clc
clear
close all

%% Merge and load measurement data
% Load measurement data
load(fullfile(pwd, 'Results_BS_Trajectory_04', 'mess_0.mat'));
Data_ges = data.signals.values;
t_ges = data.time;

%% Detect Start Time
% Reference starting Point
bool_reference = Data_ges(:, 16);

% Find first index with one
idx_start = find(bool_reference == 1, 1);
t_start = t_ges(idx_start);

%% Separate relevant data part
% ILC Time 
delta_t = 0.4;    % Trajectory lasts 1s
t_end = t_start + delta_t;

% Cut out results vectors
idx_ilc = (t_ges >= t_start & t_ges <= t_end);

t_vec_ILC = t_ges(idx_ilc);
t_vec_ILC = t_vec_ILC - t_vec_ILC(1);   % Adjust range
data_vec_ILC = Data_ges(idx_ilc, :);

%% Downsample measurement data for ILC
% Demanded sample Time (has to fit the actual sample time)
Ts_ILC = 0.001;

% Compare to actual sample time
Ts = t_vec_ILC(2) - t_vec_ILC(1);
sample_factor = round(Ts_ILC / Ts);

% Downsample
idx_ilc_down = 1:sample_factor:length(t_vec_ILC);

t_vec_ILC_down = t_vec_ILC(idx_ilc_down);
data_vec_ILC_down = data_vec_ILC(idx_ilc_down, :);

%% Save data
% Output vector
y_vec = data_vec_ILC_down(:, 1);
t_vec_meas = t_vec_ILC_down;

% Save file
res_name = fullfile(pwd, 'Results_BS_Trajectory_04', 'Results_ABS.mat');
save(res_name, 'y_vec', 't_vec_meas');

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_ges, Data_ges(:, 1), LineWidth=1); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Full Run');
legend('Location', 'best');

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec_ILC, data_vec_ILC(:, 1), LineWidth=1, DisplayName='Original'); hold on;
plot(t_vec_meas, y_vec, LineWidth=1, DisplayName='Downsampled'); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('ILC Run');
legend('Location', 'best');





