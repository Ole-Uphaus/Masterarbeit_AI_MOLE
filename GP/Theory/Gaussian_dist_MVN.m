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

args.x_cell = {x_plot};
args.y_cell = {{gauss_prob}};
args.x_label_cell = {'$\psi$'};
args.y_label_cell = {'$p_{\Psi}(\psi)$'};
args.title_cell = {''};
args.legend_cell = {{}};

args.filename = fullfile('02_Grundlagen', 'Gaußverteilung.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 6.5;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;
opts.marker = 'none';

% Create Plot
Position = [0.27, 0.20, 0.55, 0.72];
gaus_dist_plot = Plot_Manager(args);
gaus_dist_plot.single_histo_plot(opts, Position, x);

%% Multivariate Gaussian distribution
% Plot vectors
x_plot_1 = linspace(-4, 4, 20);
x_plot_2 = linspace(-4, 4, 20);
[X1, X2] = meshgrid(x_plot_1, x_plot_2);    % Create Mesh

% Parameters
mu = [0, 0];
Sigma = [1, 0.6;
    0.6, 1];

% Put mesh in vectors
X = [X1(:), X2(:)];     
invS = inv(Sigma);
detS = det(Sigma);
dists = X - mu;

% Calculate density
Z = exp(-0.5*sum((dists*invS).*dists, 2)) / (2*pi*sqrt(detS));

% Put z in mesh order
Z = reshape(Z, size(X1));

%% Plot
% Assign values (args)
args = struct();

args.x_cell = {x_plot_1};
args.y_cell = {{x_plot_2}};
args.x_label_cell = {'$\psi_1$'};
args.y_label_cell = {'$\psi_2$'};
args.title_cell = {''};
args.legend_cell = {{}};

args.filename = fullfile('02_Grundlagen', 'Multivariate_Normalverteilung.pdf');
args.save_pdf = save_pdf;

% Assign values (opts)
opts = struct();
opts.fig_height = 7;
opts.linewidth = 1.5;
opts.y_scale = 'linear';
opts.y_rel_offset = 0;
opts.x_rel_offset = 0;
opts.marker = 'none';

% Create Plot
Position = [0.30, 0.20, 0.45, 0.72];
View = [-15, 30];
z_label = '$p_{\mathbf{\Psi}}(\mathbf{\psi})$';
gaus_mvn_plot = Plot_Manager(args);
gaus_mvn_plot.single_3d_plot(opts, Position, View, Z, z_label);

%% Idenpendent gaussian vectors
% Vector size = 2
X_vec_2 = randn(2, 10);

% Vector size = 10
X_vec_10 = randn(10, 10);

%% Plot
% % Create Plot 1
% filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_2.pdf');
% plot_gaussian_vectors(X_vec_2, filename, save_pdf)
% 
% % Create Plot 2
% filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_10.pdf');
% plot_gaussian_vectors(X_vec_10, filename, save_pdf)

% Create Plot 3
filename = fullfile('02_Grundlagen', 'Unabhaengige_Zufallsvektoren_2_10.pdf');
plot_2_gaussian_vectors(X_vec_2, X_vec_10, filename, save_pdf)

%% Gaussian vectors with squared exponantial covariance
% Dependent vectors (Dimension 10)
X_vec_dep_10 = dep_gaussian_vectores(10, 10, 0.5, 1);

% Dependent vectors (Dimension 10)
X_vec_dep_100 = dep_gaussian_vectores(40, 10, 0.5, 1);

%% Plot
% % Create Plot 1
% filename = fullfile('02_Grundlagen', 'Abhaengige_Zufallsvektoren_10.pdf');
% plot_gaussian_vectors(X_vec_dep_10, filename, save_pdf)
% 
% % Create Plot 2
% filename = fullfile('02_Grundlagen', 'Abhaengige_Zufallsvektoren_100.pdf');
% plot_gaussian_vectors(X_vec_dep_100, filename, save_pdf)

% Create Plot 3
filename = fullfile('02_Grundlagen', 'Abhaengige_Zufallsvektoren_10_100.pdf');
plot_2_gaussian_vectors(X_vec_dep_10, X_vec_dep_100, filename, save_pdf)

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
    
    args.x_cell = {linspace(0, 1, size(X_vec, 1))};
    temp = cell(size(X_vec, 2), 1);
    for i = 1:size(X_vec, 2)
        temp{i} = X_vec(:, i);
    end
    args.y_cell = {temp};
    args.x_label_cell = {'y'};
    args.y_label_cell = {'\psi'};
    args.title_cell = {''};
    args.legend_cell = {{}};
    
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 6.5;
    opts.linewidth = 1.5;
    opts.y_scale = 'linear';
    opts.y_rel_offset = 0;
    opts.x_rel_offset = 0.05;
    opts.marker = '.';
    
    % Create Plot
    Position = [0.27, 0.20, 0.55, 0.72];
    plot = Plot_Manager(args);
    plot.single_plot(opts, Position);
end

function plot_2_gaussian_vectors(X_vec1, X_vec2, filename, save_pdf)
    % Assign values (args)
    args = struct();
    
    args.x_cell = {linspace(0, 1, size(X_vec1, 1)), linspace(0, 1, size(X_vec2, 1))};
    temp1 = cell(size(X_vec1, 2), 1);
    for i = 1:size(X_vec1, 2)
        temp1{i} = X_vec1(:, i);
    end
    temp2 = cell(size(X_vec2, 2), 1);
    for i = 1:size(X_vec2, 2)
        temp2{i} = X_vec2(:, i);
    end
    args.y_cell = {temp1, temp2};
    args.x_label_cell = {'$x$', '$x$'};
    args.y_label_cell = {'$\psi$', ''};
    args.title_cell = {'\textbf{a)}', '\textbf{b)}'};
    args.legend_cell = {{}, {}};
    
    args.filename = filename;
    args.save_pdf = save_pdf;
    
    % Assign values (opts)
    opts = struct();
    opts.fig_height = 7;
    opts.linewidth = 1.5;
    opts.y_scale = 'linear';
    opts.y_rel_offset = 0;
    opts.x_rel_offset = 0.05;
    opts.marker = '.';
    
    % Create Plot
    plot = Plot_Manager(args);
    orientation = [1, 2];
    plot.tiled_plot(opts, orientation);
end