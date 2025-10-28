% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      28.10.2025
% Beschreibung:
% In diesem skript werde ich ILC Regelung für eine den linearen Oszillator
% nochmal zum Vergelich erstellen.
% -------------------------------------------------------------

clc
clear

%% System Dynamics
% Simulation parameters
m  = 2; % kg
c1 = 2; % N/m
c2 = 1; % N/m^3
d  = 0.5; % Ns/m
m_delay = 1;

% State space representation x = [i; omega]
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
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
r_vec = x_max*sin(2*pi/T_end.*t_vec');

%% ILC quadratic optimal design
% Lifted system dynamics
N = size(t_vec, 2);
P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);

% Parameters
N_iter = 10;
x0 = [0;
    0]; 
W = eye(size(P));
S = 0*eye(size(P));

% Initialisation
ILC_Quadr = ILC_SISO(r_vec, m_delay);
ILC_Quadr.init_Quadr_type(W, S, P)

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[y_vec_Quadr, ~, ~] = lsim(sys_disc, u_sim, t_vec, x0);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_vec_Quadr); 0];

    % Simulate the system
    [y_vec_Quadr, ~, ~] = lsim(sys_disc, u_sim, t_vec, x0);
end

% Plot results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);   % 1 Zeile, 2 Spalten, erster Plot
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired'); hold on;
plot(t_vec, y_vec_Quadr, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('Compare desired and simulated Trajectory');
legend()

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(1:length(ILC_Quadr.RMSE_log), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()

%% Compare lifted represetation
% Parameters
x_sim = zeros(N, 2);

% Lifted Dynamics
P_nonlin = Lifted_dynamics_nonlinear_SISO(@(x) linear_discrete_system(x, Ts), N, m_delay, x_sim);

% Compare Results
delta_P = P - P_nonlin;
max_error_P = max(abs(P(:) - P_nonlin(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von P: %.3e\n', max_error_P);

%% Local Functions
function [Ad, Bd, Cd, Dd] = linear_discrete_system(x_star, Ts)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 0; % N/m^3
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