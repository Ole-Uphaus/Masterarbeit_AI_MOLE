% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      28.10.2025
% Beschreibung:
% In diesem skript werde ich die AI-MOLE IMplementierung von Michael Meindl
% ausprobieren und sein Modell für meine Systeme nutzen.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Path to Meindls implementation
base_dir = fileparts(mfilename("fullpath"));
Meindl_path = fullfile(base_dir, '..', '..', '..', '..', 'Repo_Meindl_AI_MOLE', '1_matlab_ws', '1_SISO_MOLE', '0_lib');
ILC_path = fullfile(base_dir, '..', '..', '..', 'ILC', 'ILC_SISO');
addpath(Meindl_path);
addpath(ILC_path);

%% Reference Trajectory
tic;
% Parameters
x_max = 0.5;
Ts = 0.01;
T_end = 5;

t_vec = 0:Ts:T_end;

% Trajectory (no delay - delay is applied later)
sigma = 1;
[r_vec, ~, ~] = Random_C2_trajectory_1D(2, t_vec, sigma);

% Initial input Trajectory (simple sin)
sigma_I = 0.1;
u_init = sigma_I*sin(2*pi/T_end.*t_vec');

%% AI-MOLE
% Dynamic Function
dyn_func = @(u) solve_ODE(u, t_vec);

% AI-MOLE
SISO_MOLE = CIOMOLE(dyn_func);

% Parameters (make u_init, r_vec the right dimensions)
N_trials = 11;  % Train one trial more, because there is a lag in the loop)
u_init_MOLE = u_init(1:end-1);
r_vec_MOLE = r_vec(2:end);

% Run AI-MOLE
tic;
[ev, ec, yc, uc] = SISO_MOLE.run_iomole(r_vec_MOLE, u_init_MOLE, N_trials);
toc

%% Plot Results
figure;
set(gcf, 'Position', [100 100 1200 500]);

subplot(1,2,1);
plot(t_vec, r_vec, LineWidth=1, DisplayName='desired');
hold on;
for i = 1:N_trials-1
    plot(t_vec, [0; yc{i}], LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, [0; yc{N_trials}], LineWidth=1, DisplayName=sprintf('Iteration %d', N_trials-1));
grid on;
xlabel('Zeit [s]'); 
ylabel('x [m]');
title('System Response');
legend()

subplot(1,2,2);
hold on;
for i = 1:N_trials-1
    plot(t_vec, [uc{i}; 0], LineWidth=1, Color=[0.5 0.5 0.5], DisplayName=sprintf('Iteration %d', i-1));
end
plot(t_vec, [uc{N_trials}; 0], LineWidth=1, DisplayName=sprintf('Iteration %d', N_trials-1));
grid on;
xlabel('F [N]'); 
ylabel('x [m]');
title('System Input');
legend()

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

function y = solve_ODE(u, t_vec)
    % Fill u Vector for Simulation to match Dimensions
    u_sim = [u; 0];

    % Initial Condition
    x0 = [0;
        0];

    % Step size
    Ts = t_vec(2) - t_vec(1);

    % Solver settings
    opts = odeset( ...
        'RelTol', 1e-6, ...         % Tolerance
        'AbsTol', [1e-8 1e-8], ...  % Tolerance
        'MaxStep', Ts/5, ...        % Use smaller step size for better Results
        'InitialStep', Ts/20);

    % Simulation (leave out the first output, becuase u doesnt influence
    % it)
    [~, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_sim, t_vec), t_vec, x0, opts);
    y = x_sim(2:end, 1);
end