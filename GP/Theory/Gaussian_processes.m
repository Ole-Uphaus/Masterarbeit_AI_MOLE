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

args.x_cell = {x_plot};
args.y_cell = {{y_plot}};
args.x_label_cell = {'$x$'};
args.y_label_cell = {'$f(x)$'};
args.title_cell = {''};
args.legend_cell = {{'$f(x)$', '$y_t(x_t)$'}};

args.filename = fullfile('02_Grundlagen', 'Trainingsdaten_Funktion.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0.05;
opts.x_rel_offset = 0;
opts.marker = 'none';

% Create Plot
Position = [0.27, 0.20, 0.55, 0.72];
train_plot = Plot_Manager(args);
train_plot.single_scatter_plot(opts, Position, x_train, y_train);

%% Train GP with baseline hyperparameters
% Hyperparameters
hyp_0.ell = 0.8;
hyp_0.sigma_f = 1;
hyp_0.sigma_n = 0.15;

% Run GP prediction
[mu_star_0, sigma_star_0] = gp_predict_rbf(x_train, y_train, x_plot, hyp_0);

%% Plot
% Plot baseline regression
filename = fullfile('02_Grundlagen', 'GP_baseline_Plot.pdf');
plot_single_GP_prediction(x_plot, y_plot, x_train, y_train, mu_star_0, sigma_star_0, filename, save_pdf)

%% Investigate hyperparameter influence
% Hyperparameter variation
ells = [0.3, 6.0];
sigfs = [0.4, 2.0];
signs = [0.03, 0.5];

% Lengthscale
mu_ell_cell = cell(numel(ells), 1);
sigma_ell_cell = cell(numel(ells), 1);
for i = 1:numel(ells)
    % Hyperparameter
    hyp = hyp_0;
    hyp.ell = ells(i);

    % Train GP
    [mu_star, sigma_star] = gp_predict_rbf(x_train, y_train, x_plot, hyp);
    mu_ell_cell{i} = mu_star;
    sigma_ell_cell{i} = sigma_star;
end

% Sigma_f
mu_sigf_cell = cell(numel(ells), 1);
sigma_sigf_cell = cell(numel(ells), 1);
for i = 1:numel(ells)
    % Hyperparameter
    hyp = hyp_0;
    hyp.sigma_f = sigfs(i);

    % Train GP
    [mu_star, sigma_star] = gp_predict_rbf(x_train, y_train, x_plot, hyp);
    mu_sigf_cell{i} = mu_star;
    sigma_sigf_cell{i} = sigma_star;
end

% Sigma_n
mu_sign_cell = cell(numel(ells), 1);
sigma_sign_cell = cell(numel(ells), 1);
for i = 1:numel(ells)
    % Hyperparameter
    hyp = hyp_0;
    hyp.sigma_n = signs(i);

    % Train GP
    [mu_star, sigma_star] = gp_predict_rbf(x_train, y_train, x_plot, hyp);
    mu_sign_cell{i} = mu_star;
    sigma_sign_cell{i} = sigma_star;
end

%% Plot
filename = fullfile('02_Grundlagen', 'GP_einfluss_laengenskala.pdf');
plot_tiled_GP_prediction(x_plot, x_train, y_train, mu_ell_cell, sigma_ell_cell, filename, save_pdf)

filename = fullfile('02_Grundlagen', 'GP_einfluss_sigmaf.pdf');
plot_tiled_GP_prediction(x_plot, x_train, y_train, mu_sigf_cell, sigma_sigf_cell, filename, save_pdf)

filename = fullfile('02_Grundlagen', 'GP_einfluss_sigman.pdf');
plot_tiled_GP_prediction(x_plot, x_train, y_train, mu_sign_cell, sigma_sign_cell, filename, save_pdf)

%% Optimal hyperparameter
% Use fitrgp
optimal_GP = fitrgp( ...
        x_train, y_train, ...
        'BasisFunction','none', ...
        'KernelFunction','squaredexponential', ...
        'Standardize', false);

% Extract optimal Hyperparameters
hyp_opt.ell = optimal_GP.KernelInformation.KernelParameters(1);
hyp_opt.sigma_f = optimal_GP.KernelInformation.KernelParameters(2);
hyp_opt.sigma_n = optimal_GP.Sigma;

% Run GP prediction
[mu_star_opt, sigma_star_opt] = gp_predict_rbf(x_train, y_train, x_plot, hyp_opt);

%% Plot
% Plot optimal regression
filename = fullfile('02_Grundlagen', 'GP_optimal_Plot.pdf');
plot_single_GP_prediction(x_plot, y_plot, x_train, y_train, mu_star_opt, sigma_star_opt, filename, save_pdf)

%% Local functions
function [mu_star, sigma_star] = gp_predict_rbf(x_train, y_train, x_star, hyp)
    % Compute kernel value
    K = rbf_kernel(x_train, x_train, hyp.ell, hyp.sigma_f);
    n = size(x_train, 1);
    K_y = K + (hyp.sigma_n^2)*eye(n);

    % Cholesky
    L = chol(K_y, 'lower');

    % alpha
    alpha = L'\(L\y_train);

    % Cross kernel
    K_s = rbf_kernel(x_train, x_star, hyp.ell, hyp.sigma_f);

    % Posterior mean
    mu_star = K_s'*alpha;

    % Posterior covariance
    v = L\K_s;
    K_ss = rbf_kernel(x_star, x_star, hyp.ell, hyp.sigma_f);
    Cov = K_ss - v.'*v;

    % Sandard deviation
    Var = diag(Cov);
    sigma_star = sqrt(Var);
end

function K = rbf_kernel(x1, x2, ell, sigma_f)
    % Squared distances
    Dists2 = (x1 - x2.').^2;

    % Kernel matrix
    K = sigma_f^2 * exp(-0.5 * Dists2 /(ell^2));
end

function plot_single_GP_prediction(x_plot, y_plot, x_train, y_train, mu_star, sigma_star, filename, save_pdf)
    % Assign values (args)
    args = struct();
    
    args.x_cell = {x_plot};
    args.y_cell = {{y_plot, mu_star}};
    args.x_label_cell = {'$x$'};
    args.y_label_cell = {'$f(x)$'};
    args.title_cell = {''};
    args.legend_cell = {{'3$\sigma$-Band', '$f(x)$', '$\bar{f}_*(x)$', '$y_t(x_t)$'}};
    
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 6.5;
    opts.linewidth = 1.5;
    opts.y_scale = 'linear';
    opts.y_rel_offset = 0.05;
    opts.x_rel_offset = 0;
    opts.marker = 'none';
    
    % Create Plot
    y_upper = mu_star + 3*sigma_star;
    y_lower = mu_star - 3*sigma_star;

    Position = [0.27, 0.20, 0.55, 0.72];
    regression_plot = Plot_Manager(args);
    regression_plot.single_scatter_fill_plot(opts, Position, x_train, y_train, x_plot, y_upper, y_lower);
end

function plot_tiled_GP_prediction(x_plot, x_train, y_train, mu_star_cell, sigma_star_cell, filename, save_pdf)
    % Assign values (args)
    n = numel(mu_star_cell);
    args = struct();
    
    args.x_cell = {x_plot, x_plot};
    args.y_cell = {{mu_star_cell{1}}, {mu_star_cell{2}}};
    args.x_label_cell = {'$x$', '$x$'};
    args.y_label_cell = {'${f}(x)$', ''};
    args.title_cell = {'$\ell = 0.3$', '$\ell = 6$'};
    args.legend_cell = {{}, {'3$\sigma$-Band', '$\bar{f}_*(x)$', '$y_t(x_t)$'}};
    
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 7;
    opts.linewidth = 1.5;
    opts.y_scale = 'linear';
    opts.y_rel_offset = 0.05;
    opts.x_rel_offset = 0;
    opts.marker = 'none';
    
    % Uncertainty
    y_upper_cell = cell(n, 1);
    y_lower_cell = cell(n, 1);
    x_fill_cell = cell(n, 1);

    for i = 1:n
        y_upper_cell{i} = mu_star_cell{i} + 3*sigma_star_cell{i};
        y_lower_cell{i} = mu_star_cell{i} - 3*sigma_star_cell{i};
        x_fill_cell{i} = x_plot;
    end

    % Training data
    x_train_cell = cell(n, 1);
    y_train_cell = cell(n, 1);

    for i = 1:n
        x_train_cell{i} = x_train;
        y_train_cell{i} = y_train;
    end

    % Create Plot
    orientation = [1, n];
    tiled_regression_plot = Plot_Manager(args);
    tiled_regression_plot.tiled_scatter_fill_plot(opts, orientation, x_train_cell, y_train_cell, x_fill_cell, y_upper_cell, y_lower_cell);
end