% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      21.01.2025
% Beschreibung:
% In diesem skript werde ich einige Plots für den Grundlagenteil meiner
% Masterarbeit erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(38);

% Generate Dynamic file Paths
Plot_path = fullfile(pwd, '..', '..', 'Plot');
addpath(Plot_path);

%% General
save_pdf = false;

%% Ground truth data
% GT Function
function_GT = @(x) 1*sin(0.8*x);

% Plot data
x_min = -4; x_max = 4;
x_plot = linspace(x_min, x_max, 1000)';
y_plot = function_GT(x_plot);

% Sample noisy training data
n_train = 12;
sigma_n_GT = 0.12;
x_train = sort(x_min + (x_max - x_min)*rand(n_train, 1));
y_train = function_GT(x_train) + sigma_n_GT*randn(size(x_train));

%% Plot
% Assign values (args)
args = struct();

args.x = x_plot;
args.y_cell = {{y_plot}};
args.x_label = 'x';
args.y_label_cell = {'y'};
args.title_cell = {'Title'};
args.legend_cell = {{'GT Funktion', 'Trainingsdaten'}};

args.print_legend = true;
args.filename = fullfile('02_Grundlagen', 'Trainingsdaten_Funktion.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 8;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0.05;
opts.x_rel_offset = 0;
opts.marker = 'none';

% Create Plot
train_plot = Plot_Manager(args);
train_plot.single_scatter_plot(opts, x_train, y_train);