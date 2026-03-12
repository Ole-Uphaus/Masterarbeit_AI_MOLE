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
% Compare stribeck
data_stribeck_controlled = load('MOLE_Tors_stribeck_controlled_ff_relative_Plot.mat');
data_stribeck_uncontrolled = load('MOLE_Tors_stribeck_uncontrolled_relative_Plot.mat');

%% Plot
iter_vec = 0:data_stribeck_controlled.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {data_stribeck_controlled.t_vec, iter_vec};
args.y_cell = {{data_stribeck_controlled.SISO_MOLE.ILC_SISO.r_vec, data_stribeck_controlled.SISO_MOLE.y_cell{end}, data_stribeck_uncontrolled.SISO_MOLE.y_cell{end}}, 
    {data_stribeck_controlled.SISO_MOLE.ILC_SISO.RMSE_log, data_stribeck_uncontrolled.SISO_MOLE.ILC_SISO.RMSE_log}};
args.x_label_cell = {'$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$y_{T,10}$ in $\mathrm{rad}$', 'RMSE in $\mathrm{rad}$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'Soll', 'LQR', 'ungeregelt'}, {'LQR', 'ungeregelt'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Tors_Sim_vergleich_ungeregelt.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = {'linear', 'log'};
opts.y_lim = {[], []};
opts.x_lim = {[], []};
opts.marker = 'none';

% Create Plot
plot = Plot_Manager(args);
orientation = [1, 2];
plot.tiled_plot(opts, orientation);