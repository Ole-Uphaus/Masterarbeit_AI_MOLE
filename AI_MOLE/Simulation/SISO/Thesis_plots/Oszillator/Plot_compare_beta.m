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
% Compare Heuristic, stochastic, model-based
data_heuristic = load('MOLE_Osz_linear_Meindl_Plot.mat');
data_stochastic = load('MOLE_Osz_linear_relative_Plot.mat');
data_modelbased = load(fullfile(ILC_path, 'Thesis_plots' ,'Oszillator', 'NOILC_Osz_linear_Plot.mat'));

% Compare impact of beta
data_varbeta = load('MOLE_Osz_linear_relative_beta_Plot.mat');

%% Plot
iter_vec = 0:data_heuristic.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {iter_vec, iter_vec(2:end)};
args.y_cell = {
    {data_varbeta.SISO_MOLE_Cell{1, 1}.ILC_SISO.RMSE_log, data_stochastic.SISO_MOLE.ILC_SISO.RMSE_log, data_varbeta.SISO_MOLE_Cell{2, 1}.ILC_SISO.RMSE_log, data_varbeta.SISO_MOLE_Cell{3, 1}.ILC_SISO.RMSE_log},
    {data_varbeta.SISO_MOLE_Cell{1, 1}.alpha_log, data_stochastic.SISO_MOLE.alpha_log, data_varbeta.SISO_MOLE_Cell{2, 1}.alpha_log, data_varbeta.SISO_MOLE_Cell{3, 1}.alpha_log}};
args.x_label_cell = {'Iteration', 'Iteration'};
args.y_label_cell = {'RMSE in $\mathrm{m}$', '$\eta$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'$\beta_{rel} = 0$', '$\beta_{rel} = 0{,}5$', '$\beta_{rel} = 2$', '$\beta_{rel} = 5$'}, {'$\beta_{rel} = 0$', '$\beta_{rel} = 0{,}5$', '$\beta_{rel} = 2$', '$\beta_{rel} = 5$'},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Osz_linear_vergleich_beta.pdf');
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