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
data_ABS = load(fullfile(base_dir, '..', 'Results_Adaptive_BS', 'Results_BS_Trajectory_04', 'Results_ABS.mat'));

%% Cumulative error
% Error vectors
error_MOLE = data_mole.ref_traj.phi2(:) - data_mole.SISO_MOLE.y_cell{end}(:);
error_ABS = data_mole.ref_traj.phi2(:) - data_ABS.y_vec;

% Calculate Errors
error_abs_cum_MOLE = cumtrapz(data_mole.ref_traj.t_vec, abs(error_MOLE));
error_abs_cum_ABS = cumtrapz(data_mole.ref_traj.t_vec, abs(error_ABS));

%% Plot
iter_vec = 0:data_mole.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {data_mole.ref_traj.t_vec, data_mole.ref_traj.t_vec};
args.y_cell = {{data_mole.ref_traj.phi2, data_mole.SISO_MOLE.y_cell{end}, data_ABS.y_vec}, 
    {error_abs_cum_MOLE, error_abs_cum_ABS}};
args.x_label_cell = {'$t$ in $\mathrm{s}$', '$t$ in $\mathrm{s}$'};
args.y_label_cell = {'$y_{Z}$ in $\mathrm{Nm}$', 'IAE in $\mathrm{rad \cdot s}$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'Soll', 'MOLEs', 'ABS'}, {'MOLEs', 'ABS'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Pneumatk_vergleich_MOLE_ABS.pdf');
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