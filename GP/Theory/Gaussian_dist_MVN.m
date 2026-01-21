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
args.filename = fullfile('02_Grundlagen', 'Gaußverteilung.pdf');
args.save_pdf = false;

% Assign values (opts)
opts = struct();
opts.fig_height = 8;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;
opts.marker = 'none';

% Create Plot
gaus_dist_plot = Plot_Manager(args);
gaus_dist_plot.single_histo_plot(opts, x);

%% Idenpendent gaussian vectors
% Vector size = 2
X_vec_2 = randn(2, 10);

% Assign values (args)
args = struct();

args.x = 1:size(X_vec_2, 1);
temp = cell(10, 1);
for i = 1:length(X_vec_2)
    temp{i} = X_vec_2(:, i);
end
args.y_cell = {temp};
args.x_label = 'x';
args.y_label_cell = {'y'};
args.title_cell = {'Title'};
args.legend_cell = {{''}};

args.print_legend = false;
args.filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_2.pdf');
args.save_pdf = false;

% Assign values (opts)
opts = struct();
opts.fig_height = 8;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;
opts.marker = '.';

% Create Plot
gaus_ind_plot_2 = Plot_Manager(args);
gaus_ind_plot_2.single_plot(opts);

% Vector size = 10
X_vec_10 = randn(10, 10);

% Assign values (args)
args = struct();

args.x = 1:size(X_vec_10, 1);
temp = cell(10, 1);
for i = 1:length(X_vec_10)
    temp{i} = X_vec_10(:, i);
end
args.y_cell = {temp};
args.x_label = 'x';
args.y_label_cell = {'y'};
args.title_cell = {'Title'};
args.legend_cell = {{''}};

args.print_legend = false;
args.filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_10.pdf');
args.save_pdf = false;

% Assign values (opts)
opts = struct();
opts.fig_height = 8;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;
opts.marker = '.';

% Create Plot
gaus_ind_plot_2 = Plot_Manager(args);
gaus_ind_plot_2.single_plot(opts);