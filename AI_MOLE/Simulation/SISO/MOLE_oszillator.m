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

%% Reference Trajectory
tic;
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Generate Dynamic file Path
base_dir = fileparts(mfilename("fullpath"));
ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');
addpath(ILC_path);

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);

%% Initialize AI-MOLE
% Parameters
m_delay = 1;
N_iter = 5;
H_trials = 3;
x0 = [0;
    0];

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Initial input Trajectory (simple sin)
sigma_I = 0.1;
u_init = sigma_I*sin(2*pi/T_end.*t_vec');

% Initialisation
SISO_MOLE = SISO_MOLE_IO(r_vec, m_delay, u_init, N_iter, H_trials);

%% Run ILC
% Update Loop
u_sim = u_init;
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    u_sim = [SISO_MOLE.update_input(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1);
end
SISO_MOLE.save_final_trajectory(y_sim);
y_sim_quadratic = y_sim;

%% Plot Results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:N_iter
    plot(t_vec, SISO_MOLE.y_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, SISO_MOLE.y_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend('Location', 'best');

subplot(2,2,3);   % 1 Zeile, 2 Spalten, erster Plot
plot(0:(length(SISO_MOLE.ILC_SISO.RMSE_log)-1), SISO_MOLE.ILC_SISO.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,4);   % 1 Zeile, 2 Spalten, erster Plot
hold on;
for i = 1:N_iter
    plot(t_vec, SISO_MOLE.u_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, SISO_MOLE.u_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Input Signal');
legend('Location', 'best');

% Zeitmessung und Ausgabe
time = toc;
fprintf('Dauer von AI-MOLE mit %d Iterationen, %d Trials Delay (H) und jeweils %d Datenpunkten pro Trial: %g s\n', N_iter, H_trials, length(t_vec), time);

%% Local Functions
function dx = oszillator_nonlinear(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 2; % N/m^3
    d  = 0.5; % Ns/m

    % States
    x = x_vec(1);
    xp = x_vec(2);

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % Dynamics
    dx = zeros(2, 1);
    dx(1) = xp;
    dx(2) = 1/m*(-c1*x - c2*x^3 - d*xp + u);
end