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
% Controlled 
data_controlled_heuristic = load('MOLE_Tors_controlled_meindl_Plot.mat');
data_controlled_minimize = load('MOLE_Tors_controlled_minimize_Plot.mat');
data_controlled_relative = load('MOLE_Tors_controlled_relative_Plot.mat');

%% Plot
iter_vec = 0:data_controlled_heuristic.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {iter_vec, iter_vec(2:end)};
args.y_cell = {{data_controlled_heuristic.SISO_MOLE.ILC_SISO.RMSE_log, data_controlled_relative.SISO_MOLE.ILC_SISO.RMSE_log, data_controlled_minimize.SISO_MOLE.ILC_SISO.RMSE_log}, 
    {data_controlled_heuristic.SISO_MOLE.alpha_log, data_controlled_relative.SISO_MOLE.alpha_log, data_controlled_minimize.SISO_MOLE.alpha_log}};
args.x_label_cell = {'Iteration', 'Iteration'};
args.y_label_cell = {'RMSE in $\mathrm{rad}$', '$\eta$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'MOLE', 'REL', 'MIN'}, {'MOLE', 'REL', 'MIN'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Tors_Sim_vergleich_min_rel.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = {'log', 'linear'};
opts.y_lim = {[], []};
opts.x_lim = {[], [0, 10]};
opts.marker = {'none', '.'};

% Create Plot
plot = Plot_Manager(args);
orientation = [1, 2];
plot.tiled_plot(opts, orientation);