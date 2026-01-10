% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      05.01.2025
% Beschreibung:
% In diesem skript werde ich die ILC Regelung für den linearen
% Torsionsschwinger wie am Prüfstand implementieren. ich werde dazu das
% lineare Systemmodell nutzen und auch das geregelte system betrachten.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(Model_Path);

%% Reference Trajectory
filename = 'Trajectory_02.mat';
filepath = fullfile(pwd, '..', '..', '..', 'AI_MOLE', 'Test_Bench', 'Torsion_Oszillator', 'Reference_Trajectories', filename);
load(filepath);

r_vec = ref_traj.phi2;
r_vec = r_vec(:);

%% System dynamics
% Simulation parameters
J1  = 0.0299;    % kgm^2
J2  = 0.0299;    % kgm^2
c_phi = 7.309;   % Nm/rad
d_v1 = 0.055;    % Nms/rad
d_v2 = 0.0064;   % Nms/rad

% State space
A = [0, 1, 0, 0;
    -c_phi/J1, -d_v1/J1, c_phi/J1, 0;
    0, 0, 0, 1;
    c_phi/J2, 0, -c_phi/J2, -d_v2/J2];
B = [0;
    1/J1;
    0;
    0];
C = [0, 0, 1, 0];
D = 0;

%% Parameters
% ILC architcture ('uncontrolled', 'serial_LQR', 'serial_LQR_ff')
architecture = 'uncontrolled';

% Simulation parameters
Ts = ref_traj.t_vec(2) - ref_traj.t_vec(1);
T_end = ref_traj.t_vec(end);
t_vec = ref_traj.t_vec;

% Solver settings
opts = odeset( ...
    'RelTol', 1e-6, ...         % Tolerance
    'AbsTol', 1e-8, ...  % Tolerance
    'MaxStep', Ts/5, ...        % Use smaller step size for better Results
    'InitialStep', Ts/20);

% ILC parameter
N = size(t_vec, 2);
m_delay = 1;

%% ILC Initialisation
switch architecture
    case 'uncontrolled'
        % Use the uncontrolled System for model based ILC

        % System dynamics
        sys_cont = ss(A,B,C,D);
        sys_disc = c2d(sys_cont, Ts, 'zoh');
        [Ad,Bd,Cd,Dd] = ssdata(sys_disc);

        % Choose system model
        system_dynamics = @torsion_oszillator_linear;
        x0 = [0; 0; 0; 0];

        % Lifted system dynamics
        P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        S = 0.1*eye(size(P));
        R = 0*eye(size(P));

        % Initial Input Trajectory
        u_init = zeros(size(r_vec, 1), 1);
        
        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    case 'serial_LQR'
        % Use the controlled System for model based ILC

        % Controlled System Dynamics
        k_T = [12.524133472585133, 1.268619349231718, -9.207508682229756, 0.314246813584626];
        A_cont = A - B*k_T;
        B_cont = B;
        C_cont = C;
        D_cont = D;

        % System dynamics
        sys_cont = ss(A_cont,B_cont,C_cont,D_cont);
        sys_disc = c2d(sys_cont, Ts, 'zoh');
        [Ad,Bd,Cd,Dd] = ssdata(sys_disc);

        % Choose system model
        system_dynamics = @torsion_oszillator_linear_LQR;
        x0 = [0; 0; 0; 0];

        % Lifted system dynamics
        P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        S = 0.00000001*eye(size(P));
        R = 0*eye(size(P));

        % Initial Input Trajectory
        u_init = zeros(size(r_vec, 1), 1);
        
        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    case 'serial_LQR_ff'
        % Use the controlled System for model based ILC

        % Controlled System Dynamics (with static feedforward gain)
        k_T = [12.524133472585133, 1.268619349231718, -9.207508682229756, 0.314246813584626];
        S_gain = 3.316624790355372;
        A_cont = A - B*k_T;
        B_cont = B * S_gain;
        C_cont = C;
        D_cont = D;

        % System dynamics
        sys_cont = ss(A_cont,B_cont,C_cont,D_cont);
        sys_disc = c2d(sys_cont, Ts, 'zoh');
        [Ad,Bd,Cd,Dd] = ssdata(sys_disc);

        % Choose system model
        system_dynamics = @torsion_oszillator_linear_LQR_ff;
        x0 = [0; 0; 0; 0];

        % Lifted system dynamics
        P = Lifted_dynamics_linear_SISO(Ad, Bd, Cd, N, m_delay);
        
        % Parameters
        N_iter = 10;
        W = eye(size(P));
        S = 0.001*eye(size(P));
        R = 0*eye(size(P));

        % Initial Input Trajectory
        u_init = r_vec;
        
        % Initialisation
        ILC_Quadr = ILC_SISO(r_vec, m_delay, u_init, N_iter);
        ILC_Quadr.init_Quadr_type(W, S, R, P)
        % ILC_Quadr.init_Q_lowpass(Q_fc, Q_order, Ts);

    otherwise
        error('Unbekannte Architektur ausgewählt.')
end

%% Run ILC

% Update Loop
u_sim = [ILC_Quadr.u_vec; 0];
[t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
y_sim = x_sim(:, 3);
for i = 1:N_iter
    % Update input
    u_sim = [ILC_Quadr.Quadr_update(y_sim); 0];

    % Simulate the system
    [t_sim, x_sim] = ode45(@(t,x) system_dynamics(t, x, u_sim, t_vec), t_vec, x0, opts);
    y_sim = x_sim(:, 3);
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
% plot(t_vec, w_vec, LineWidth=1, DisplayName='proc'); hold on;
% plot(t_vec, w_rep_vec, LineWidth=1, DisplayName='proc rep');
% plot(t_vec, v_vec, LineWidth=1, DisplayName='meas');
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