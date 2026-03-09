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

%% Load Results
% Compare Heuristic, stochastic, model-based
data_heuristic = load('MOLE_Osz_linear_Meindl_Plot.mat');
data_stochastic = load('MOLE_Osz_linear_relative_Plot.mat');

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {};
args.y_cell = {};
args.x_label_cell = {'', '', '$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$y_L$ in $\mathrm{m}$', 'RMSE in $\mathrm{m}$', '$u_L$ in $\mathrm{N}$', '$\eta$'};
args.title_cell = {'', '', '', ''};
args.legend_cell = {{'$y_{L,d}$', '$y_{L,0}$', '$y_{L,5}$', '$y_{L,10}$'}, {}, {'$u_{L,0}$', '$u_{L,5}$', '$u_{L,10}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_linear_relative.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_lim = {[], [], [], []};
opts.x_lim = {[], [], [], [0, 10]};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
log_scale = true;   % Set true for log error plot
plot.tiled_mole_results_plot(opts, SISO_MOLE, t_vec, log_scale);