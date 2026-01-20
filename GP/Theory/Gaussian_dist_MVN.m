% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      20.01.2025
% Beschreibung:
% In diesem skript werde ich einige Plots für den Grundlagenteil meiner
% Masterarbeit erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(42);

% Generate Dynamic file Paths
Plot_path = fullfile(pwd, '..', '..', 'Plot');
addpath(Plot_path);

%% Gaussian Distribution
% Random data
N = 1000;
mu = 0;
sigma = 1;

x = mu + sigma*randn(N, 1);

% Gaussian probability density
x_max_abs = max(abs(x));
x_plot = linspace(-x_max_abs, x_max_abs, 1000);
gauss_prob = normpdf(x_plot, mu, sigma);

%% Plot
% Assign values (args)
args = struct();

args.x = x_plot;
args.y_cell = {{gauss_prob}};
args.x_label = 'x';
args.y_label_cell = {'y'};
args.title_cell = {'Title'};
args.legend_cell = {{''}};

args.print_legend = false;
args.filename = [];
args.save_pdf = [];

% Assign values (opts)
opts = struct();
opts.fig_height = 8;
opts.linewidth = 1.5;
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;

% Create Plot
gaus_dist_plot = Plot_Manager(args);
gaus_dist_plot.single_histo_plot(opts, x);
