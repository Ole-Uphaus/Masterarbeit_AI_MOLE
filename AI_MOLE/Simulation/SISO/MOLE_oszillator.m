% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.11.2025
% Beschreibung:
% In diesem skript werde ich AI-MOLE für den nichtlinearen (bzw. wenn man
% den nichtlinearen Koeffizienten auf null setzt - linearen) Oszillator.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(ILC_path);
addpath(Model_Path);

%% Reference Trajectory
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;
x0 = [0;
    0];

t_vec = 0:Ts:T_end;

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Noise Parameters
sigma_v = 0.0;      % Measurement Noise 0.01
fc_v = 20;
white = true;       % if white == true -> white noise is sampled - no filter

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);

%% System model

% Choose system model ('linear', 'nonlinear', 'nonlinear_stribeck',)
system_model = 'nonlinear';

switch system_model
    case 'linear'
        % Use the linear system for AI-MOLE
        dynamic_model = @oszillator_linear;

        % Initial input Trajectory (simple sin or automatic generated)
        sigma_I = 0.1;  % for stibeck model = 1, otherwise = 0.1
        u_init_sin = sigma_I*sin(2*pi/T_end.*t_vec');
        u_init = u_init_sin;        % u_init_sin / u_init_auto

    case 'nonlinear'
        % Use the nonlinear system for AI-MOLE
        dynamic_model = @oszillator_nonlinear;

        % Initial input Trajectory (simple sin or automatic generated)
        sigma_I = 0.1;  % for stibeck model = 1, otherwise = 0.1
        u_init_sin = sigma_I*sin(2*pi/T_end.*t_vec');
        u_init = u_init_sin;        % u_init_sin / u_init_auto

    case 'nonlinear_stribeck'
        % Use the nonlinear system with stribeck friction for AI-MOLE
        dynamic_model = @oszillator_nonlinear_stribeck;

        % Initial input Trajectory (simple sin or automatic generated)
        sigma_I = 1;  % for stibeck model = 1, otherwise = 0.1
        u_init_sin = sigma_I*sin(2*pi/T_end.*t_vec');
        u_init = u_init_sin;        % u_init_sin / u_init_auto

    otherwise
        error('Unbekanntes Systemmodell ausgewählt.')

end

%% Determine initial input Trajectory
% Parameters
alpha = 0.3;   % Initial alpha (percentage of amplitude of reference signal)
growth_factor = 1.5;
max_iter = 10;

% Criterium var(y) > k * var(noise)
k_sign_var = 5;

% Run Loop
for i = 1:max_iter
    % Calculate sigmaI from alpha and reference amplitude
    sigmaI = alpha * max(abs(r_vec));

    % Initial input Trajectory
    u_init_auto = Initial_input_trajectory(r_vec, Ts, sigmaI);

    % System Simulation
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
    [t_sim, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_init_auto, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;

    % Calculate output variance
    y_var = var(y_sim);

    fprintf('Iter %d: alpha=%.3f, Var(y)=%.4e, NoiseVar=%.4e\n', ...
            i, alpha, y_var, sigma_v^2);

    % Test if variance is high enough
    if y_var >= k_sign_var * sigma_v^2
        fprintf('Ausreichende Anregung erreicht.\n\n');
        break;
    else
        % Update alpha
        alpha = alpha * growth_factor;
    end
end

%% Initialize AI-MOLE

params = struct();

% Parameters
params.m_delay = 1;
params.N_iter = 10;
params.H_trials = 3;

% Choose weight initialisation Method ('Meindl', 'Stochastic', 'Heuristic',
% 'Robust', 'Manual')
params.weight_init_method = 'Stochastic';

% Choose nonlinearity damping method ('none', 'relative_1', 'relative_2', 'minimize')
params.nonlin_damping = 'relative_2';
params.beta = 0.5;

% Initialisation
SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

%% Run ILC
tic;
% Update Loop
u_sim = u_init;
v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
[t_sim, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1) + v_vec;
for i = 1:params.N_iter
    % Update input
    u_sim = [SISO_MOLE.update_input(y_sim); 0];

    % Simulate the system
    v_vec = Gen_noise_Butter(t_vec, sigma_v, fc_v, white);
    [t_sim, x_sim] = ode45(@(t,x) dynamic_model(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1) + v_vec;
end
SISO_MOLE.save_final_trajectory(y_sim);
y_sim_quadratic = y_sim;

% Zeitmessung und Ausgabe
time = toc;
fprintf('Dauer von AI-MOLE mit %d Iterationen, %d Trials Delay (H) und jeweils %d Datenpunkten pro Trial: %g s\n', params.N_iter, params.H_trials, length(t_vec), time);

%% Plot Results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:params.N_iter
    % plot(t_vec, SISO_MOLE.y_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
end
plot(t_vec, SISO_MOLE.y_cell{params.N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', params.N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend('Location', 'best');

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
% plot(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
semilogy(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,4);   % 1 Zeile, 2 Spalten, erster Plot
hold on;
for i = 1:params.N_iter
    plot(t_vec, SISO_MOLE.u_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
end
plot(t_vec, SISO_MOLE.u_cell{params.N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', params.N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Input Signal');
legend('Location', 'best');

subplot(2,2,2);
plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Noise');
legend()