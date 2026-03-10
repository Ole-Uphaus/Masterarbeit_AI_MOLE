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
data_tors_mole_no_ff = load(fullfile(base_dir, '..', 'Runs', '2026_01_27', 'Run_01_serial.mat'));

% Compare Modelbased
data_tors_modelbased_ff = load(fullfile(base_dir, '..', '..', '..', '..', 'ILC', 'Test_Bench', 'Torsion_Oszillator', 'Runs', '2026_01_14', 'Run_01_serial.mat'));
data_tors_modelbased_no_ff = load(fullfile(base_dir, '..', '..', '..', '..', 'ILC', 'Test_Bench', 'Torsion_Oszillator', 'Runs', '2026_01_14', 'Run_05_serial.mat'));

% Set last input value to the before value
data_tors_mole_ff.SISO_MOLE.u_cell{end}(end) = data_tors_mole_ff.SISO_MOLE.u_cell{end}(end-1);
data_tors_mole_no_ff.SISO_MOLE.u_cell{11}(end) = data_tors_mole_no_ff.SISO_MOLE.u_cell{11}(end-1);

data_tors_modelbased_ff.ILC_Quadr.u_cell{end}(end) = data_tors_modelbased_ff.ILC_Quadr.u_cell{end}(end-1);
data_tors_modelbased_no_ff.ILC_Quadr.u_cell{end}(end) = data_tors_modelbased_no_ff.ILC_Quadr.u_cell{end}(end-1);

%% Plot
iter_vec = 0:data_tors_mole_ff.SISO_MOLE.N_iter;

% Assign values (args)
args = struct();

args.x_cell = {data_tors_modelbased_ff.ref_traj.t_vec, iter_vec};
args.y_cell = {{data_tors_mole_ff.SISO_MOLE.u_cell{end}, data_tors_mole_no_ff.SISO_MOLE.u_cell{11}, data_tors_modelbased_ff.ILC_Quadr.u_cell{end}, data_tors_modelbased_no_ff.ILC_Quadr.u_cell{end}}, 
    {data_tors_mole_ff.SISO_MOLE.ILC_SISO.RMSE_log, data_tors_mole_no_ff.SISO_MOLE.ILC_SISO.RMSE_log(1:11), data_tors_modelbased_ff.ILC_Quadr.RMSE_log, data_tors_modelbased_no_ff.ILC_Quadr.RMSE_log}};
args.x_label_cell = {'$t$ in $\mathrm{s}$', 'Iteration'};
args.y_label_cell = {'$u_{T,V}$ in $\mathrm{Nm}$', 'RMSE in $\mathrm{rad}$'};
args.title_cell = {'\textbf{(a)}', '\textbf{(b)}'};
args.legend_cell = {{}, {'MOLEsff', 'MOLEs', 'NOILCff', 'NOILC'}};

args.filename = fullfile('05_Ergebnisse_Diskussion', 'Ergebnis_Tors_Pruef_vergleich_ff_no_ff.pdf');
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