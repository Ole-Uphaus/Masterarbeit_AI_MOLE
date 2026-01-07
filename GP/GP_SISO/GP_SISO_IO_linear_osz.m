% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      31.10.2025
% Beschreibung:
% In diesem skript werde ich die Klasse für den Gauß-Prozess evaluieren und
% erproben, sodass sie zukünftig auch für AI_MOLE eingesetzt werden kann.
% Dazu werde ich versuchen, das Systemverhalten eines linearen dynmiaschen
% repetetiven Systems durch einen Gauß Prozess zu erlernen und anschließend
% Prädiktionen durchzuführen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', 'ILC', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', 'System_Models');
addpath(ILC_path);
addpath(Model_Path);

%% Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

% Noise Parameters
sigma_v = 0.00;      % Measurement Noise 0.5
fc_v = 20;
white = true;       % if white == true -> white noise is sampled - no filter

% Input trajectorys
u_scale_train = [1];
N_traj = length(u_scale_train);
u_scale_test = 2;

u_vec_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    u_vec_train_cell{i} = u_scale_train(i)*sin(2*pi/T_end.*t_vec');
end

u_vec_test= u_scale_test*sin(2*pi/T_end.*t_vec');
% u_vec_test = u_scale_test * (t_vec' >= 1) .* (t_vec' <= 3);

%% Data generation
% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation train/test
x0 = [0;
    0];

y_sim_train_cell = cell(N_traj, 1);
for i = 1:N_traj
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
    [~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_train_cell{i}, t_vec), t_vec, x0, opts);
    y_sim_train_cell{i} = x_sim(:, 1) + v_vec;
end

[~, x_sim] = ode45(@(t,x) oszillator_linear(t, x, u_vec_test, t_vec), t_vec, x0, opts);
y_sim_test = x_sim(:, 1);

%% Predict System Dynamics
% Init Gaussian Process
GP_IO = GP_SISO_IO();

% Train Gaussian Process
GP_IO.train_GP_model(y_sim_train_cell, u_vec_train_cell);

% Predict new Trajectory (with measurement noise)
[y_pred_test, y_std_test] = GP_IO.predict_trajectory_measurement(u_vec_test);
% y_pred_test_upper = y_pred_test + 2*y_std_test;
% y_pred_test_lower = y_pred_test - 2*y_std_test;

% Predict new Trajectory with full covriance matrix (without measurement noise)
[~, Cov_y_test] = GP_IO.predict_trajectory_covariance(u_vec_test);
y_pred_test_upper = y_pred_test + 2*sqrt(diag(Cov_y_test));
y_pred_test_lower = y_pred_test - 2*sqrt(diag(Cov_y_test));

% Error predicted variance (with measurement noise)
Var_y_test = y_std_test.^2;
error_var_test = Var_y_test - (diag(Cov_y_test) + GP_IO.GP.Sigma^2);
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von Var_y: %.3e \n\n', max(abs(error_var_test)));

%% Linearize GP at given input trajectory
% Linearisation
tic;
P = GP_IO.linearize_at_given_trajectory(u_vec_train_cell{1});
t = toc;
fprintf('Dauer der Linearisierung: %g s\n', t);

tic;
[P2, Cov_dy_dv_cell2, Cov_dy_du_cell2] = GP_IO.linearize_covariance_at_given_trajectory_fast(u_vec_train_cell{1});
t = toc;
fprintf('Dauer der Linearisierung und Varianzberechnung (fast): %g s\n\n', t);

% tic;
% [P3, Cov_dy_dv_cell3] = GP_IO.approx_linearisation_at_given_trajectory(u_vec_train_cell{1});
% t = toc;
% fprintf('Dauer der Linearisierung (approximiert): %g s\n\n', t);

%% Compare Errors
% Compare Fast and slow Computation
error_P1 = P - P2;
max_error_P1 = max(abs(P(:) - P2(:)));
max_rel_error_P1 = max(abs(P(:) - P2(:)) ./ max(abs(P(:)), eps));
fprintf('Maximaler absoluter (relativer) Unterschied bei der Bestimmung von P (fast vs. slow): %.3e (%.3e)\n\n', max_error_P1, max_rel_error_P1);

% % Compare analytic and approx Computation (Linearisation P)
% error_P3 = P - P3;
% max_error_P2 = max(abs(P(:) - P3(:)));
% max_rel_error_P2 = max(abs(P(:) - P3(:)) ./ max(abs(P(:)), eps));
% fprintf('Maximaler absoluter (relativer) Unterschied bei der Bestimmung von P (analytic vs. approx): %.3e (%.3e)\n\n', max_error_P2, max_rel_error_P2);
% 
% % Compare analytic and approx Computation (Covariance C)
% max_error_Cov_P = 0;
% max_rel_error_Cov_P = 0;
% 
% sum_abs_diff_cov = 0;
% sum_rel_diff_cov = 0;
% n_total_etries = 0;
% 
% for i = 1:GP_IO.N
%     C_analytic = Cov_dy_dv_cell2{i};
%     C_approx = Cov_dy_dv_cell3{i};
% 
%     diff_C = C_analytic - C_approx;
% 
%     % absolute error
%     max_error_Cov_P = max(max_error_Cov_P, max(abs(diff_C(:))));
% 
%     % relative error
%     rel_diff = abs(diff_C(:)) ./ max(abs(C_analytic(:)), eps);
%     max_rel_error_Cov_P = max(max_rel_error_Cov_P, max(rel_diff));
% 
%     % mean abs error
%     sum_abs_diff_cov = sum_abs_diff_cov + sum(abs(diff_C(:)));
%     sum_rel_diff_cov = sum_rel_diff_cov + sum(abs(rel_diff(:)));
%     n_total_etries = n_total_etries + numel(diff_C);
% end
% 
% mean_abs_diff_Cov_P = sum_abs_diff_cov / n_total_etries;
% mean_rel_diff_Cov_P = sum_rel_diff_cov / n_total_etries;
% 
% fprintf('Maximaler absoluter (relativer) Unterschied bei der Bestimmung von Cov_P (analytic vs. approx): %.3e (%.3e)\n', max_error_Cov_P, max_rel_error_Cov_P);
% fprintf('Mittlerer absoluter (relativer) Unterschied bei der Bestimmung von Cov_P (analytic vs. approx): %.3e (%.3e)\n\n', mean_abs_diff_Cov_P, mean_rel_diff_Cov_P);

% Prediction with linearized gp model
delta_u = u_vec_test - u_vec_train_cell{1};
y_pred_lin_test = GP_IO.predict_trajectory_measurement(u_vec_train_cell{1}) + P*delta_u;

% Calculate linearisation uncertanty as Variance of a new GP
Var_delta_y = linearisation_prediction_variance(GP_IO, delta_u, Cov_dy_dv_cell2);
Sigma_delta_y = sqrt(abs(Var_delta_y));

y_pred_lin_test_upper2 = y_pred_lin_test + 2*Sigma_delta_y;
y_pred_lin_test_lower2 = y_pred_lin_test - 2*Sigma_delta_y;

% %% Estimate Uncertainty in P (delta_P)
% % Use variance of every Element in P
% delta_P1 = element_wise_variance(GP_IO, Cov_dy_dv_cell2);
% 
% % Compare Spectral morm
% norm_P = norm(P, 2);
% norm_delta_P = norm(delta_P1, 2);
% 
% % Approx norm
% norm_delta_P2 = approx_norm(GP_IO, Cov_dy_dv_cell2);
% 
% % Sample covariance Matrix random
% [norms_delta_P, delta_P_cell] = sample_deltaP_norms(GP_IO, Cov_dy_dv_cell2, 100);
% [norm_delta_P3, idx] = max(norms_delta_P);
% delta_P3 = delta_P_cell{idx};
% 
% % Norm from bachelor Thesis
% norm_delta_P4 = norm_bachelor_thesis(GP_IO, Cov_dy_dv_cell2);
% 
% fprintf(['Spektralnorm von P: %.3e, \n' ...
%     'Spektralnorm von delta_P (elementweise Satndardabweichung): %.3e, \n' ...
%     'Spektralnorm von delta_P (vektoriesierung Lifted Matrix): %.3e, \n' ...
%     'Spektralnorm von delta_P (Sampling basiert): %.3e, \n' ...
%     'Spektralnorm von delta_P (aus Bachelorarbeit): %.3e \n\n'], norm_P, norm_delta_P, norm_delta_P2, norm_delta_P3, norm_delta_P4);
% 
% % Calculate linearisation uncertaincy based on Sampling
% y_pred_lin_test_upper3 = y_pred_lin_test + delta_P3*abs(delta_u);
% y_pred_lin_test_lower3 = y_pred_lin_test - delta_P3*abs(delta_u);

%% Compare Lifted System representation

N = length(u_vec_train_cell{1});

% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
d  = 0.5; % Ns/m

% State space representation 
A = [0, 1;
    -c1/m, -d/m];
B = [0;
    1/m];
C = [1, 0];
D = 0;

% Discrete System
sys_cont = ss(A,B,C,D);
sys_disc = c2d(sys_cont, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sys_disc);

% Calculate analytic lifted Matrix
P_analytic = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, 1);

% Compare lifted Matrices
P_reduced_size = P(2:end, 1:end-1);
error_P2 = P_analytic - P_reduced_size;
max_error_P2 = max(abs(P_analytic(:) - P_reduced_size(:)));
mean_error_P2 = mean(abs(P_analytic(:) - P_reduced_size(:)));
norm_error_P2 = norm(abs(error_P2), 2);

fprintf('Maximaler absoluter Unterschied P_analytic und P_lin: %.3e\n', max_error_P2);
fprintf('Mittlerer absoluter Unterschied P_analytic und P_lin: %.3e\n', mean_error_P2);
fprintf('Spektralnorm des absoluten Unterschieds zwischen P_analytic und P_lin: %.3e\n\n', norm_error_P2);

% Calculate prediction Error
error_y_pred = abs(y_sim_test - y_pred_test);
error_y_lin = abs(y_sim_test - y_pred_lin_test);

% %% Monotonic Convergence Condition
% % Compute bound
% % bound = monotonic_convergence_condition(P, norm_error_P2, delta_P3(2:end, 1:end-1));
% bound = monotonic_convergence_condition(P, norm_error_P2, error_P2);
% 
% fprintf('Monotone Konvergenzbedingung (abhängig von unsicherheitsschätzung): %.3e (< 1)\n', bound);

%% Plot
% 1. Plot
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);
imagesc(error_P2);
colorbar;
title('error P_{analytic} - P_{lin}');
xlabel('Input-Index');
ylabel('Output-Index');

subplot(1,2,2);
plot(t_vec, error_y_pred, LineWidth=2, DisplayName='error (y-sim - y-pred)');
hold on;
plot(t_vec, error_y_lin, LineWidth=2, DisplayName='error (y-sim - y-pred-lin)');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Prediction Error');
legend()

% 2. Plot
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test, LineWidth=1, DisplayName='y-pred-test');
fill([t_vec, fliplr(t_vec)], [y_pred_test_upper', fliplr(y_pred_test_lower')], ...
     [0.4 0.7 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
     'DisplayName','±2σ-Band');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and Predicted Results with confidence (2*sigma)');
legend()

subplot(2,2,2);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_lin_test, LineWidth=1, DisplayName='y-pred-lin-test');
% fill([t_vec, fliplr(t_vec)], [y_pred_lin_test_upper3', fliplr(y_pred_lin_test_lower3')], ...
%      [0.4 0.7 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
%      'DisplayName','±2σ-Band-geschätzt');
fill([t_vec, fliplr(t_vec)], [y_pred_lin_test_upper2', fliplr(y_pred_lin_test_lower2')], ...
     [1 0.6 0.4], 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
     'DisplayName','±2σ-Band-mathematisch');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and linearised Results with confidence (2*sigma)');
legend()

subplot(2,2,3);
plot(t_vec, y_sim_test, LineWidth=1, DisplayName='y-sim-test'); hold on;
plot(t_vec, y_pred_test, LineWidth=1, DisplayName='y-pred-test');
plot(t_vec, y_pred_lin_test, LineWidth=1, DisplayName='y-pred-lin-test');
for i = 1:N_traj
    plot(t_vec, y_sim_train_cell{i}, LineWidth=1, DisplayName=sprintf('y-sim-train%d', i));
end
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated and Predicted Results');
legend()

subplot(2,2,4);
plot(t_vec, u_vec_test, LineWidth=1, DisplayName='u-test'); hold on;
for i = 1:N_traj
    plot(t_vec, u_vec_train_cell{i}, LineWidth=1, DisplayName=sprintf('u-train%d', i));
end

grid on;
xlabel('Zeit [s]'); 
ylabel('u [N]');
title('Training and Testing Input u');
legend()

%% Local Functions
function Var_delta_y = linearisation_prediction_variance(GP_IO, delta_u, Cov_dy_dv_cell)

    % Empty variance Vector
    Var_delta_y = zeros(GP_IO.N, 1);
    
    % Construct test matrix
    V_test = GP_IO.constr_test_matrix(delta_u);

    % Calculate Variances
    for i = 1:GP_IO.N
        % Extract regression Vector
        vn = V_test(i, :);
        vn = vn(:);

        % Calculate Variance
        Var_delta_y(i) = vn.' * Cov_dy_dv_cell{i} * vn;
    end
end

function delta_P = element_wise_variance(GP_IO, Cov_dy_dv_cell)

    % Empty variance Matrix
    Var_P = zeros(GP_IO.N);
    
    for i = 1:GP_IO.N
        % Extract variances
        Var_dy_dv = diag(Cov_dy_dv_cell{i});
    
        % Update jacobi matrix
        if i > 1
            Var_P(i, 1:(i-1)) = Var_dy_dv((i-1):-1:1);
        end
    end
    
    % Calculate standard deviation of P
    Sigma_P = sqrt(Var_P);
    
    % Calculate uncertyincy by using confidence (3*sigma)
    delta_P = 3*Sigma_P;
end

function norm_delta_P = approx_norm(GP_IO, Cov_dy_dv_cell)
    
    % Initialize sum
    trace_sum = 0;

    for i = 1:GP_IO.N
        trace_sum = trace_sum + sum(diag(Cov_dy_dv_cell{i}));
    end

    % Approx norm
    norm_delta_P = sqrt(trace_sum);
end

function norm_delta_P = norm_bachelor_thesis(GP_IO, Cov_dy_dv_cell)
    
    % Empty variance Matrix
    Var_P = zeros(GP_IO.N);
    
    for i = 1:GP_IO.N
        % Extract variances
        Var_dy_dv = diag(Cov_dy_dv_cell{i});
    
        % Update jacobi matrix
        if i > 1
            Var_P(i, 1:(i-1)) = Var_dy_dv((i-1):-1:1);
        end
    end

    B = sqrt(Var_P);

    rowEnergies = sqrt(sum(B.^2, 2)); 
    sigma = max(rowEnergies);

    sigma_star = max(B, [], "all");

    norm_delta_P = sigma + sigma_star * sqrt(log(GP_IO.N));
    % norm_delta_P = sigma_star * sqrt(GP_IO.N);
end

function [norms_deltaP, delta_P_cell] = sample_deltaP_norms(GP_IO, Cov_dy_dv_cell, N_samples)

    % Empty vector and cell
    norms_deltaP = zeros(N_samples,1);
    delta_P_cell = cell(GP_IO.N, 1);

    % Empty cell for precomputed cholesky
    L_cell = cell(GP_IO.N, 1);

    % Pre calculate cholesky
    for i = 1:GP_IO.N
        % Numerically symmetric covariance
        C = Cov_dy_dv_cell{i};
        C = (C + C.')/2;
        % C = diag(Cov_dy_dv_cell{i});
        % C = diag(C);

        % Cholesky
        L_cell{i} = chol(C, 'lower');
    end

    for k = 1:N_samples
        % Empty realisation matrix
        DeltaP = zeros(GP_IO.N);
        z = randn(size(C,1),1);

        for i = 1:GP_IO.N          
            % Cholesky sample random gradient
            % z = randn(size(C,1),1);
            g_i = L_cell{i} * z;

            % Map random gradient values to lifted matrix
            if i > 1
                DeltaP(i, 1:(i-1)) = g_i((i-1):-1:1);
            end
        end

        % Spektral norm of current realisation
        norms_deltaP(k) = norm(DeltaP, 2);
        delta_P_cell{k} = DeltaP;
    end
end

function bound = monotonic_convergence_condition(P, norm_delta_P, delta_P)

    % Bring P to the right size
    P = P(2:end, 1:end-1);

    % Weighting matrix
    Q = eye(size(P));

    % A = P' Q P
    lambda = 0;      % for robust inversion
    A = P.'*Q*P + lambda*eye(size(P,2));

    % B = P' Q
    B = P.'*Q;
    % B = P.'*Q*delta_P;

    % Constant term X
    inv_A = inv(A);
    X = inv_A * B;

    % Calculate bound (using uncertainty norm)
    bound = norm(X, 2)*norm_delta_P;

    % Calculate bound (using uncertainty matrix)
    % bound = norm(X, 2);
end
