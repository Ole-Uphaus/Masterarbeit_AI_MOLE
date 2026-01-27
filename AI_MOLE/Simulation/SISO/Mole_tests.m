% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      27.01.2025
% Beschreibung:
% In diesem skript werde ich einige Tests im Rahmen des Prüfstandsversuches
% durchführen
% -------------------------------------------------------------

clc
clear
close all

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(Model_Path);

%% Reference Trajectory
filename = 'Trajectory_02.mat';
filepath = fullfile(pwd, '..', '..', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
load(filepath);

r_vec = ref_traj.phi2;
r_vec = r_vec(:);

%% Parameters
% Simulation parameters
Ts = ref_traj.t_vec(2) - ref_traj.t_vec(1);
T_end = ref_traj.t_vec(end);
t_vec = ref_traj.t_vec;

% Process Noise
sigma_w = 0.2;

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', 1e-8, ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

%% MOLE initialisation
% Choose system model
system_dynamics = @torsion_oszillator_linear_LQR_stribeck;
x0 = [0; 0; 0; 0];

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
params.beta = 0;

% Initial input Trajectory (simple sin or automatic generated)
sigma_I = 20;
u_init = sigma_I*sin(2*pi/T_end.*t_vec');

% Static gain (feedforward control signal)
% S = 3.316624790355372;
% 
% u_init = S .* r_vec;

% Initialisation
SISO_MOLE = SISO_MOLE_IO(r_vec, u_init, params);

%% Run ILC
tic;
% Update Loop
u_sim = u_init;

[t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 3);
for i = 1:params.N_iter
    % Update input
    u_sim = [SISO_MOLE.update_input(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 3);
end
SISO_MOLE.save_final_trajectory(y_sim);

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
% plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Noise');
legend()