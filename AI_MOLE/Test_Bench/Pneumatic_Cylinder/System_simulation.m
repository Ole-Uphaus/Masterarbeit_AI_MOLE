% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.01.2026
% Beschreibung:
% Dieses Skript dient dazu, nach einem Simulationsdurchlauf in Simulink die
% daten für AI-MOLE vorzubereiten und abzuspeichern.

clc
close all

%% Load data
time_ges = simout.time;
y_vec_ges = simout.signals.values;

%% Separate relevant data part
% ILC Time
t_start = 0;
t_end = 1;

% Cut out results vectors
idx_ilc = (time_ges >= t_start & time_ges <= t_end);

t_vec_ILC = time_ges(idx_ilc);
t_vec_ILC = t_vec_ILC - t_vec_ILC(1);   % Adjust range
y_vec_ILC = y_vec_ges(idx_ilc, :);

%% Downsample measurement data for ILC
% Demanded sample Time (has to fit the actual sample time)
Ts_ILC = 0.01;

% Compare to actual sample time
Ts = t_vec_ILC(2) - t_vec_ILC(1);
sample_factor = round(Ts_ILC / Ts);

% Downsample
idx_ilc_down = 1:sample_factor:length(t_vec_ILC);

t_vec_ILC_down = t_vec_ILC(idx_ilc_down);
y_vec = y_vec_ILC(idx_ilc_down);

%% Save data
% Output vector
t_vec_meas = t_vec_ILC_down;

% Timestamp of trial
sim_trial_timestamp = datetime('now');

% Save file
date_string = datestr(sim_trial_timestamp, 'yyyy_mm_dd');
res_name = sprintf('Trial_%s.mat', date_string);
save(res_name, 'y_vec', 'sim_trial_timestamp', 't_vec_meas');

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(time_ges, y_vec_ges, LineWidth=1); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('Full Run');
legend('Location', 'best');

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec_ILC, y_vec_ILC, LineWidth=1, DisplayName='Original'); hold on;
plot(t_vec_meas, y_vec, LineWidth=1, DisplayName='Downsampled'); hold on;
grid on;
xlabel('Zeit [s]'); 
ylabel('phi2 [rad]');
title('ILC Run');
legend('Location', 'best');
