% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.01.2026
% Beschreibung:
% Dieses Skript dient dazu, nach einem Durchlauf am Prüfstand die aktuelle
% Ausgangstrajektorie aus dem entsprechenden Messdatenordner zu Laden und
% als bearbeitetes file in diesem Verzeichnis abzuspeichern. Auf Basis von
% diesem File kann anschließend der AI-MOLE step durchgeführt werden.
% -------------------------------------------------------------

clc
clear
close all

%% Find the right directory
% Path to base dict
base_dir = fullfile(pwd, '..', '..', '..', '..');

% Find measurement folders
folder_names = dir(fullfile(base_dir, 'Messungen_*'));

% Extract numbers of folders
numbers  = nan(numel(folder_names), 1);

for i = 1:numel(folder_names)
    % Slpit name
    parts = split(folder_names(i).name, '_');

    % Save number
    numbers(i) = str2double(parts{2});
end

% Find folder with largest number
[~, idx] = max(numbers);
latest_folder_name = folder_names(idx).name;
latest_folder_path = fullfile(base_dir, latest_folder_name);

%% Merge and load measurement data
% Use script to merge data
run('C:\Users\Pneumatik\Desktop\ILC_Uphaus\Run_Zusammenfueren.m')

% Load measurement data
load(fullfile(latest_folder_path, 'Messungen_gesamt.mat'));

%% Detect Start Time
% Reference starting Point
bool_reference = Data_ges(:, 16);

% Find first index with one
idx_start = find(bool_reference == 1, 1);
t_start = t_ges(idx_start);

%% Separate relevant data part
% ILC Time 
delta_t = 1;    % Trajectory lasts 1s
t_end = t_start + delta_t;

% Cut out results vectors
idx_ilc = (time_ges >= t_start & time_ges <= t_end);

t_vec_ILC = time_ges(idx_ilc);
t_vec_ILC = t_vec_ILC - t_vec_ILC(1);   % Adjust range
data_vec_ILC = Data_ges(idx_ilc, :);

%% Downsample measurement data for ILC
% Demanded sample Time (has to fit the actual sample time)
Ts_ILC = 0.01;

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

% Timestamp of trial
sim_trial_timestamp = datetime(TimeStamp, 'InputFormat', 'dd-MM-yyyy HH:mm:ss');

% Save file
date_string = datestr(sim_trial_timestamp, 'yyyy_mm_dd');
res_name = sprintf('Trial_%s.mat', date_string);
save(res_name, 'y_vec', 'sim_trial_timestamp', 't_vec_meas');

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(time_ges, Data_ges(:, 1), LineWidth=1); hold on;
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





