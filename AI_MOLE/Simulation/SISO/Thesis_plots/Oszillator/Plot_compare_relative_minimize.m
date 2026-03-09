% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      09.03.2026
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
ILC_path = fullfile(base_dir, '..', '..', '..', '..', '..', 'ILC', 'Simulation', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'System_Models');
MOLE_Path = fullfile(base_dir, '..', '..');
Plot_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'Plot');
GP_path = fullfile(base_dir, '..', '..', '..', '..', '..', 'GP', 'GP_SISO');
addpath(ILC_path);
addpath(Model_Path);
addpath(MOLE_Path);
addpath(Plot_Path);
addpath(GP_path);

%% General
save_pdf = false;

%% Load Results
% Compare nonlinear
data_nonlin_minimize = load('MOLE_Osz_nonlinear_minimize_Plot.mat');
data_nonlin_relative = load('MOLE_Osz_nonlinear_relative_Plot.mat');
data_nonlin_modelbased = load(fullfile(ILC_path, 'Thesis_plots' ,'Oszillator', 'NOILC_Osz_nonlinear_Plot.mat'));

% Compare stribeck
data_stribeck_minimize = load('MOLE_Osz_nonlinear_stribeck_minimize_Plot.mat');
data_stribeck_relative = load('MOLE_Osz_nonlinear_stribeck_relative_Plot.mat');
data_stribeck_modelbased = load(fullfile(ILC_path, 'Thesis_plots' ,'Oszillator', 'NOILC_Osz_nonlinear_stribeck_Plot.mat'));

%% Plot
iter_vec_nonlin = 0:data_nonlin_minimize.SISO_MOLE.N_iter;
iter_vec_stribeck = 0:data_stribeck_minimize.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {iter_vec_nonlin, iter_vec_nonlin(2:end), iter_vec_stribeck, iter_vec_stribeck(2:end)};
args.y_cell = {{data_nonlin_relative.SISO_MOLE.ILC_SISO.RMSE_log, data_nonlin_minimize.SISO_MOLE.ILC_SISO.RMSE_log, data_nonlin_modelbased.ILC_Quadr.RMSE_log}, 
    {data_nonlin_relative.SISO_MOLE.alpha_log, data_nonlin_minimize.SISO_MOLE.alpha_log},
    {data_stribeck_relative.SISO_MOLE.ILC_SISO.RMSE_log, data_stribeck_minimize.SISO_MOLE.ILC_SISO.RMSE_log, data_stribeck_modelbased.ILC_Quadr.RMSE_log},
    {data_stribeck_relative.SISO_MOLE.alpha_log, data_stribeck_minimize.SISO_MOLE.alpha_log}};
args.x_label_cell = {'', '', 'Iteration', 'Iteration'};
args.y_label_cell = {'RMSE in $\mathrm{m}$', '$\eta$', 'RMSE in $\mathrm{m}$', '$\eta$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}', '\textbf{(c)}', '\textbf{(d)}'};
args.legend_cell = {{'rel', 'min', 'NOILC'}, {'rel', 'min'}, {'rel', 'min', 'NOILC'}, {'rel', 'min'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_nichtlin_vergleich.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = {'log', 'linear', 'log', 'linear'};
opts.y_lim = {[], [0, 1], [], [0, 1]};
opts.x_lim = {[], [0, 15], [], [0, 20]};
opts.marker = {'none', '.', 'none', '.'};

% Create Plot
plot = Plot_Manager(args);
orientation = [2, 2];
plot.tiled_plot(opts, orientation);