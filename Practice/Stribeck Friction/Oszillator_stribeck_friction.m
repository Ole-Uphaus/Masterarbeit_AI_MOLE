% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      18.11.2025
% Beschreibung:
% In diesem Skript werde ich das Systemmodell des nichtlinearen Oszillators
% um eine Stribeck reibung erweitern und erste untersuchungen zum
% Systemverhalten durchführen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

%% Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;
N = length(t_vec);

% Input Trajectory
u_scale = 2;

u_sim= u_scale*sin(2*pi/T_end.*t_vec');

%% System Simulation
% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', [1e-8 1e-8], ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% Simulation without friction
x0 = [0;
    0];

[~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 1);

% Simulation with friction
[~, x_sim] = ode45(@(t,x) oszillator_nonlinear_stribeck(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim_fric = x_sim(:, 1);

%% Calculate friction forces
% Friction Parameters (Stribeck)
Fc   = 0.4;     % Coulomb-Friction [N]
Fs   = 0.5;     % Static Friction [N] (Fs > Fc)
vs   = 0.1;     % Stribeck-Velocity [m/s]
dv   = 0.1;     % viskos damping constant [Ns/m]
v_eps = 0.001;   % Regularisation tanh()

% General friction graph
v_plot = linspace(0, 1, N);
F_fric_general = (Fc + (Fs - Fc)*exp(-(v_plot/vs).^2)) .* tanh(v_plot/v_eps) + dv*v_plot;

% Friction force
v_sim_fric = x_sim(:, 2);
F_fric = (Fc + (Fs - Fc)*exp(-(v_sim_fric/vs).^2)) .* tanh(v_sim_fric/v_eps) + dv*v_sim_fric;

%% Plots
figure;
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1);
plot(v_plot, F_fric_general, LineWidth=1, DisplayName='F-stribeck'); hold on;
grid on;
xlabel('v [m/s]'); 
ylabel('F [N]');
title('General Friction Force');
legend()

subplot(2,2,3);
plot(t_vec, y_sim, LineWidth=1, DisplayName='y-sim'); hold on;
plot(t_vec, y_sim_fric, LineWidth=1, DisplayName='y-sim-frict');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Simulated Results');
legend()

subplot(2,2,4);
plot(t_vec, u_sim, LineWidth=1, DisplayName='u-sim'); hold on;
plot(t_vec, F_fric, LineWidth=1, DisplayName='F-frict');
grid on;
xlabel('Zeit [s]'); 
ylabel('F [N]');
title('Forces on oszillator');
legend()

%% Local Functions
function dx = oszillator_nonlinear_stribeck(t, x_vec, u_vec, t_vec)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 2; % N/m^3
    d  = 0.5; % Ns/m

    % Friction Parameters (Stribeck)
    Fc   = 0.4;     % Coulomb-Friction [N]
    Fs   = 0.5;     % Static Friction [N] (Fs > Fc)
    vs   = 0.1;     % Stribeck-Velocity [m/s]
    dv   = 0.1;     % viskos damping constant [Ns/m]
    v_eps = 0.001;   % Regularisation tanh()

    % States
    x = x_vec(1);
    xp = x_vec(2);

    % Input
    u = interp1(t_vec, u_vec, t, 'previous', 'extrap');

    % Stribeck-friction
    F_fric = (Fc + (Fs - Fc)*exp(-(xp/vs)^2)) * tanh(xp/v_eps) + dv*xp;
    % F_fric = 0;

    % Dynamics
    dx = zeros(2, 1);
    dx(1) = xp;
    dx(2) = 1/m*(-c1*x - c2*x^3 - d*xp - F_fric + u);
end

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