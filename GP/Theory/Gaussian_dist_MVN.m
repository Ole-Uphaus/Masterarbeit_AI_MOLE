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

%% Geeral
save_pdf = false;

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
args.x_label_cell = {'x'};
args.y_label_cell = {'y'};
args.title_cell = {'Title'};
args.legend_cell = {{''}};

args.print_legend = false;
args.filename = fullfile('02_Grundlagen', 'Gaußverteilung.pdf');
args.save_pdf = save_pdf;

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

% Vector size = 10
X_vec_10 = randn(10, 10);

%% Plot
% Create Plot 1
filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_2.pdf');
plot_gaussian_vectors(X_vec_2, filename, save_pdf)

% Create Plot 2
filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_10.pdf');
plot_gaussian_vectors(X_vec_10, filename, save_pdf)

%% Gaussian vectors with squared exponantial covariance
% Dependent vectors (Dimension 10)
X_vec_dep_10 = dep_gaussian_vectores(10, 10, 0.5, 1);

% Dependent vectors (Dimension 10)
X_vec_dep_100 = dep_gaussian_vectores(100, 10, 0.5, 1);

%% Plot
% Create Plot 1
filename = fullfile('02_Grundlagen', 'Abhaengige_Zufallsvektoren_10.pdf');
plot_gaussian_vectors(X_vec_dep_10, filename, save_pdf)

% Create Plot 2
filename = fullfile('02_Grundlagen', 'Abhaengige_Zufallsvektoren_100.pdf');
plot_gaussian_vectors(X_vec_dep_100, filename, save_pdf)

%% Local functions
function X_vec_dep = dep_gaussian_vectores(size_x, n_x, ell, sigma)
    % Idizes that are close together get small covariances
    t = linspace(0, 1, size_x)';  
    
    % Squared exponantial covariance
    Dists2 = (t - t.').^2;
    K = sigma^2 * exp(-0.5 * Dists2 / (ell^2));
    
    % Cholesky
    jitter = 1e-10;
    L = chol(K + jitter*eye(size_x), 'lower');
    
    % Dependent samples
    X_vec_dep = L * randn(size_x, n_x);
end

function plot_gaussian_vectors(X_vec, filename, save_pdf)
    % Assign values (args)
    args = struct();
    
    args.x = linspace(0, 1, size(X_vec, 1));
    temp = cell(size(X_vec, 2), 1);
    for i = 1:size(X_vec, 2)
        temp{i} = X_vec(:, i);
    end
    args.y_cell = {temp};
    args.x_label_cell = {'y'};
    args.y_label_cell = {'x'};
    args.title_cell = {'Title'};
    args.legend_cell = {{''}};
    
    args.print_legend = false;
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 8;
    opts.linewidth = 1.5;
    opts.y_scale = 'linear';
    opts.y_rel_offset = 0;
    opts.x_rel_offset = 0;
    opts.marker = '.';
    
    % Create Plot
    plot = Plot_Manager(args);
    plot.single_plot(opts);
end