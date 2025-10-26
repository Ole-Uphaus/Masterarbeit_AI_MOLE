% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      23.10.2025
% Beschreibung:
% In diesem skript werde ich eine einfache ILC Regelung für einen DC-Motor
% implementieren (linear, SISO).
% -------------------------------------------------------------

clc
clear

%% System Dynamics
% Simulation parameters
J  = 1.5e-3;     % kg m^2
d  = 1.2e-3;     % N m s/rad
R  = 2.0;        % Ohm
L  = 5.0e-3;     % H
Kt = 0.08;       % N m / A
Ke = 0.08;       % V s / rad
m_delay = 1;     % Time-step delay of the System

% State space representation x = [i; omega]
A = [-R/L, -Ke/L;
    Kt/J, -d/J];
B = [1/L;
    0];
C = [0, 1];
D = 0;

sys_cont = ss(A,B,C,D);
G_cont = tf(sys_cont);

% Discrete state space
Ts = 1.e-2;
sys_disc = c2d(sys_cont, Ts, 'zoh');
[Ad,Bd,Cd,Dd] = ssdata(sys_disc);

%% Reference Trajectory
% Parameters
omega_max = 100;
T_end = 2;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
r_vec = omega_max*sin(2*pi/T_end.*t_vec');

%% ILC PD-Type design
% Parameters
kp = 0.1;
kd = 0.05;
N_iter = 10;
x0 = [0;
    0]; 

% Initialisation
ILC_PD = ILC_linear_SISO(r_vec, m_delay);
ILC_PD.init_PD_type(kp, kd);

% Update Loop
u_sim = [ILC_PD.u_vec; 0];
[y_vec_pd, ~, ~] = lsim(sys_disc, u_sim, t_vec, x0);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_PD.PD_update(y_vec_pd); 0];

    % Simulate the system
    [y_vec_pd, ~, ~] = lsim(sys_disc, u_sim, t_vec, x0);
end

%% ILC quadratic optimal design
% Lifted system dynamics
N = size(t_vec, 2);
P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);

% Parameters
W = eye(size(P));
S = 1.5*eye(size(P));

% Initialisation
ILC_Quadr = ILC_linear_SISO(r_vec, m_delay);
ILC_Quadr.init_Quadr_type(P, W, S)

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
plot(t_vec, y_vec_pd, LineWidth=1, DisplayName='ILC PD');
plot(t_vec, y_vec_Quadr, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Zeit [s]'); 
ylabel('\Omega [rad/s]');
title('Compare desired and simulated Trajectory');
legend()

subplot(1,2,2);   % 1 Zeile, 2 Spalten, erster Plot
plot(1:length(ILC_PD.RMSE_log), ILC_PD.RMSE_log, LineWidth=1, DisplayName='ILC PD'); hold on;
plot(1:length(ILC_Quadr.RMSE_log), ILC_Quadr.RMSE_log, LineWidth=1, DisplayName='ILC Quadr');
grid on;
xlabel('Iteration'); 
ylabel('RMSE');
title('Compare error development');
legend()