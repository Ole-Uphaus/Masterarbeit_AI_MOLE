% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      26.10.2025
% Beschreibung:
% In diesem skript werde ich ILC Regelung für eine einfaches nichtlineares
% System entwerfen.
% -------------------------------------------------------------

clc
clear

%% Reference Trajectory
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
r_vec = x_max*sin(2*pi/T_end.*t_vec');

%% ILC PD-Type design
% Parameters
kp = 0.1;
kd = 0.05;
m_delay = 1;
N_iter = 10;
x0 = [0;
    0];

% Initialisation
ILC_PD = ILC_SISO(r_vec, m_delay);
ILC_PD.init_PD_type(kp, kd);

% Update Loop
u_sim = [ILC_PD.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, Ts), t_vec, x0);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_PD.PD_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, Ts), t_vec, x0);
    y_sim = x_sim(:, 1);
end
y_sim_pd = y_sim;

%% ILC quadratic optimal design
% Parameters
N = size(t_vec, 2);

W = eye(N-m_delay);
S = 0*eye(N-m_delay);

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay);
ILC_Quadr.init_Quadr_type(W, S)

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, Ts), t_vec, x0, opts);
y_sim = x_sim(:, 1);
for i = 1:N_iter
    % Update input
    P = Lifted_dynamics_nonlinear_SISO(@(x) linear_discrete_system(x, Ts), N, m_delay, x_sim);
    u_sim = [ILC_Quadr.Quadr_update(y_sim, P); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, Ts), t_vec, x0, opts);
    y_sim = x_sim(:, 1);
end
y_sim_quadratic = y_sim;

% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, y_sim_pd, LineWidth=1, DisplayName='ILC PD');
plot(t_vec, y_sim_quadratic, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(1:length(ILC_Quadr.RMSE_log), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr'); hold on;
plot(1:length(ILC_PD.RMSE_log), ILC_PD.RMSE_log, LineWidth=1, DisplayName='ILC PD');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

%% Local Functions
function dx = oszillator_nonlinear(t, x_vec, u_vec, Ts)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 2; % N/m^3
    d  = 0.5; % Ns/m

    % Get current time index
    k = floor(t/Ts) + 1;
    k = max(1, min(k, numel(u_vec)));   % Check if k ist valid

    % States
    x = x_vec(1);
    xp = x_vec(2);

    % Input
    u = u_vec(k);

    % Dynamics
    dx = zeros(2, 1);
    dx(1) = xp;
    dx(2) = 1/m*(-c1*x - c2*x^3 - d*xp + u);
end

function [Ad, Bd, Cd, Dd] = linear_discrete_system(x_star, Ts)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 2; % N/m^3
    d  = 0.5; % Ns/m

    % States
    x = x_star(1);
    xp = x_star(2);

    % Linearisation
    A_lin = [0, 1;
        (-c1/m - 3*c2/m*x^2), -d/m];
    B_lin = [0;
        1/m];
    C_lin = [1, 0];
    D_lin = 0;
    
    sys_cont = ss(A_lin, B_lin, C_lin, D_lin);
    
    % Discrete
    sys_disc = c2d(sys_cont, Ts, 'zoh');
    [Ad, Bd, Cd, Dd] = ssdata(sys_disc);
end