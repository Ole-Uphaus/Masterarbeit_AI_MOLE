% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      17.12.2025
% Beschreibung:
% Dieses Skript dient dazu, Solltrajektorien für AI-MOLE zu erstellen und
% diese im entsprechenden ordner abzuspeichern.
% -------------------------------------------------------------

clc
clear
close all

%% Trajectory Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;

% Time Points
T_acc = 1.5 * T_end/5;      % Acceleration time
T_hold = 2 * T_end/5;       % Constant velocity time
T_dec = 1.5 * T_end/5;      % Deceleration time

% Constant velocityS
omega = deg2rad(90);    % deg/s --> rad/s

% Save trajectory
save_traj = false;

%% Reference trajectory
% Time points
time_points = [ ...
    0, ...
    T_acc, ...
    T_acc + T_hold, ...
    T_acc + T_hold + T_dec];

% Waypoints (for acceleration and deceleration use mean velocity)
phi0 = 0;
phi1 = phi0 + 0.5*omega*T_acc;
phi2 = phi1 + omega*T_hold;
phi3 = phi2 + 0.5*omega*T_dec;

waypoints = [phi0, phi1, phi2, phi3];

% Velocity and acceleration boundaries
vel_bounds = [0, omega, omega, 0];
acc_bounds = [0, 0, 0, 0];

[phi2, phi2_p, phi2_pp] = quinticpolytraj( ...
    waypoints, time_points, t_vec, ...
    'VelocityBoundaryCondition', vel_bounds, ...
    'AccelerationBoundaryCondition', acc_bounds);

%% Save reference trajectory
if save_traj
    % Number of .mat-files in dir
    files = dir(fullfile(pwd, '*.mat'));
    N_traj = numel(files);
    
    % Create new file path
    filename = sprintf('Trajectory_%02d.mat', N_traj+1);
    filepath = fullfile(pwd, filename);

    % Create struct
    ref_traj = struct();

    ref_traj.t_vec = t_vec;
    ref_traj.phi2 = phi2;
    ref_traj.phi2_p = phi2_p;
    ref_traj.phi2_pp = phi2_pp;
    
    % Save data
    save(filepath, 'ref_traj');
end

%% Plots
figure();
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1)
plot(t_vec, phi2, 'LineWidth', 1.5); grid on;
ylabel('phi (rad/s)')
title('Winkel')

subplot(2,2,2)
plot(t_vec, phi2_p, 'LineWidth', 1.5); grid on;
ylabel('phi-p (rad/s)')
xlabel('t (s)')
title('Winkelgeschwindigkeit')

subplot(2,2,3)
plot(t_vec, phi2_pp, 'LineWidth', 1.5); grid on;
ylabel('phi-pp (rad/s^2)')
xlabel('t (s)')
title('Winkelbeschleunigung')
