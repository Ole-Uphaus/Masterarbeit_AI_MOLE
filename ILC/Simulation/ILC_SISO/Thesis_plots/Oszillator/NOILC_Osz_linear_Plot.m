% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      24.02.2026
% Beschreibung:
% In diesem skript werde ich einen NO-ILC Plot für die Masterarbeit
% erstellen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', '..', '..', '..', 'System_Models');
addpath(Model_Path);

%% System Dynamics
% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
d  = 0.5; % Ns/m
m_delay = 1;

% State space representation 
A = [0, 1;
    -c1/m, -d/m];
B = [0;
    1/m];
C = [1, 0];
D = 0;

sys_cont = ss(A,B,C,D);
G_cont = tf(sys_cont);

% Discrete state space
Ts = 1.e-2;
sys_disc = c2d(sys_cont, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sys_disc);

%% Reference Trajectory
% Parameters
x_max = 0.5;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);
u_init = zeros(size(r_vec, 1), 1);

%% ILC quadratic optimal design
% Lifted system dynamics
N = size(t_vec, 2);
P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);

% Parameters
N_iter = 10;
x0 = [0;
    0]; 
W = eye(size(P));
S = 0.00007*eye(size(P));
R = 0*eye(size(P));

% Track Results
y_cell_quadr = cell(N_iter+1, 1);
u_cell_quadr = cell(N_iter+1, 1);

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
ILC_Quadr.init_Quadr_type(W, S, R, P)

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim), t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_linear(t, x, (u_sim), t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 1);
end
% Calculate and log final error
ILC_Quadr.calculate_final_error(y_sim);

%% Plot results
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
for i = 1:N_iter
    plot(t_vec, ILC_Quadr.y_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
end
plot(t_vec, ILC_Quadr.y_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(2,2,3);
% plot(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
semilogy(0:(length(ILC_Quadr.RMSE_log)-1), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

subplot(2,2,2);
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare process and measurement noise');
legend()

subplot(2,2,4);
hold on;
for i = 1:N_iter
    plot(t_vec, ILC_Quadr.u_cell{i}, LineWidth=1, Color=[0.5 0.5 0.5], HandleVisibility='off');
end
plot(t_vec, ILC_Quadr.u_cell{N_iter+1}, LineWidth=1, DisplayName=sprintf('Iteration %d', N_iter));
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Input Signal');
legend()