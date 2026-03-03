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
T_end = 1;
t_vec = 0:Ts:T_end;

% Endposition
x_end = 0.4;

% Save trajectory
save_traj = false;

%% Reference trajectory
% Time points
time_points = [0, T_end];

% Waypoints
waypoints = [0.1, x_end];

% Velocity and acceleration boundaries
vel_bounds = [0, 0];
acc_bounds = [0, 0];

[phi2, phi2_p, phi2_pp, pp] = quinticpolytraj( ...
    waypoints, time_points, t_vec, ...
    'VelocityBoundaryCondition', vel_bounds, ...
    'AccelerationBoundaryCondition', acc_bounds);

% 3rd Derivative
c1 = pp.coefs(2, 1);
c2 = pp.coefs(2, 2);
c3 = pp.coefs(2, 3);
phi2_ppp = 60*c1*t_vec.^2 + 24*c2*t_vec + 6*c3;

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
    ref_traj.phi2_ppp = phi2_ppp;
    
    % Save data
    save(filepath, 'ref_traj');
end

%% Plots
figure();
set(gcf, 'Position', [100 100 1200 800]);

subplot(2,2,1)
plot(t_vec, phi2, 'LineWidth', 1.5); grid on;
ylabel('phi (rad/s)')
title('Position')

subplot(2,2,2)
plot(t_vec, phi2_p, 'LineWidth', 1.5); grid on;
ylabel('phi-p (rad/s)')
xlabel('t (s)')
title('Geschwindigkeit')

subplot(2,2,3)
plot(t_vec, phi2_pp, 'LineWidth', 1.5); grid on;
ylabel('phi-pp (rad/s^2)')
xlabel('t (s)')
title('Beschleunigung')

subplot(2,2,4)
plot(t_vec, phi2_ppp, 'LineWidth', 1.5); grid on;
ylabel('phi-ppp (rad/s^3)')
xlabel('t (s)')
title('Ruck')