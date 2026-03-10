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
% Compare MOLEs
data_tors_mole_ff = load(fullfile(base_dir, '..', 'Runs', '2026_01_13', 'Run_05_serial.mat'));

% Compare Modelbased
data_tors_modelbased_ff = load(fullfile(base_dir, '..', '..', '..', '..', 'ILC', 'Test_Bench', 'Torsion_Oszillator', 'Runs', '2026_01_14', 'Run_01_serial.mat'));

%% Cumulative error
% Error vectors
error_mole = data_tors_modelbased_ff.ref_traj.phi2(:) - data_tors_mole_ff.SISO_MOLE.y_cell{end}(:);
error_modelbased = data_tors_modelbased_ff.ref_traj.phi2(:) - data_tors_modelbased_ff.ILC_Quadr.y_cell{end}(:);

% Calculate Errors
error_abs_cum_mole = cumtrapz(data_tors_modelbased_ff.ref_traj.t_vec, abs(error_mole));
error_abs_cum_modelbased = cumtrapz(data_tors_modelbased_ff.ref_traj.t_vec, abs(error_modelbased));

%% Plot
iter_vec = 0:data_tors_mole_ff.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {data_tors_modelbased_ff.ref_traj.t_vec, data_tors_modelbased_ff.ref_traj.t_vec};
args.y_cell = {{data_tors_mole_ff.ref_traj.phi2, data_tors_mole_ff.SISO_MOLE.y_cell{end}, data_tors_modelbased_ff.ILC_Quadr.y_cell{end}}, 
    {error_abs_cum_mole, error_abs_cum_modelbased}};
args.x_label_cell = {'$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$y_{T}$ in $\mathrm{Nm}$', 'IAE in $\mathrm{rad \cdot s}$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{'Soll', 'MOLEs', 'NOILC'}, {'MOLEs', 'NOILC'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Tors_Pruef_vergleich_kum_fehler.pdf');
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