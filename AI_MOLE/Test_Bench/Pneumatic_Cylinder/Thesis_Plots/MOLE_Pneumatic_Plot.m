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
data_mole = load(fullfile(base_dir, '..', 'Runs', '2026_03_05', 'Run_01.mat'));

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {};
args.y_cell = {};
args.x_label_cell = {'', '', '$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$y_Z$ in $\mathrm{m}$', 'RMSE in $\mathrm{m}$', '$u_{Z,V}$ in $\mathrm{N}$', '$\eta$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}', '\textbf{(c)}', '\textbf{(d)}'};
args.legend_cell = {{'$y_{Z,d}$', '$y_{Z,0}$', '$y_{Z,10}$', '$y_{Z,20}$'}, {}, {'$u_{Z,V0}$', '$u_{Z,V10}$', '$u_{Z,V20}$'}, {},};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Pneumatk_MOLE.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 10;
opts.linewidth = 1.5;
opts.y_scale = {'linear', 'log', 'linear', 'linear'};
opts.y_lim = {[], [], [], [0, 1]};
opts.x_lim = {[], [], [], [0, 20]};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
plot.tiled_mole_results_plot(opts, data_mole.SISO_MOLE, data_mole.ref_traj.t_vec);