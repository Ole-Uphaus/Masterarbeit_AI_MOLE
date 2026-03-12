% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      10.03.2026
% Beschreibung:
% In diesem skript werde ich einen AI-MOLE Plot für die Masterarbeit
% erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', '..', '..', 'ILC', 'Simulation', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', '..', 'System_Models');
MOLE_Path = fullfile(base_dir, '..', '..', '..', 'Simulation', 'SISO');
Plot_Path = fullfile(base_dir, '..', '..', '..', '..', 'Plot');
GP_path = fullfile(base_dir, '..', '..', '..', '..', 'GP', 'GP_SISO');
addpath(ILC_path);
addpath(Model_Path);
addpath(MOLE_Path);
addpath(Plot_Path);
addpath(GP_path);

%% General
save_pdf = false;

%% Load Results
% Compare algorithms
data_tors_meindl = load(fullfile(base_dir, '..', 'Runs', '2026_01_14', 'Run_01_serial.mat'));
data_tors_relative = load(fullfile(base_dir, '..', 'Runs', '2026_01_13', 'Run_05_serial.mat'));
data_tors_minimize = load(fullfile(base_dir, '..', 'Runs', '2026_01_13', 'Run_02_serial.mat'));

% Compare Ts
data_tors_Ts_05 = load(fullfile(base_dir, '..', 'Runs', '2026_01_14', 'Run_03_serial.mat'));
data_tors_Ts_10 = load(fullfile(base_dir, '..', 'Runs', '2026_01_14', 'Run_04_serial.mat'));

%% Plot
iter_vec = 0:data_tors_meindl.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {iter_vec, iter_vec};
args.y_cell = {{data_tors_meindl.SISO_MOLE.ILC_SISO.RMSE_log, data_tors_relative.SISO_MOLE.ILC_SISO.RMSE_log, data_tors_minimize.SISO_MOLE.ILC_SISO.RMSE_log}, 
    {data_tors_relative.SISO_MOLE.ILC_SISO.RMSE_log, data_tors_Ts_05.SISO_MOLE.ILC_SISO.RMSE_log, data_tors_Ts_10.SISO_MOLE.ILC_SISO.RMSE_log}};
args.x_label_cell = {'Iteration', 'Iteration'};
args.y_label_cell = {'RMSE in $\mathrm{rad}$', ''};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'MOLE', 'REL', 'MIN'}, {'$T_s^{mole} = 0{,}01\,\mathrm{s}$', '$T_s^{mole} = 0{,}05\,\mathrm{s}$', '$T_s^{mole} = 0{,}1\,\mathrm{s}$'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Tors_Pruef_vergleich_alg_Ts.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = {'linear', 'linear'};
opts.y_lim = {[], []};
opts.x_lim = {[], []};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
orientation = [1, 2];
plot.tiled_plot(opts, orientation);