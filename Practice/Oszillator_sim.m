% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      26.10.2025
% Beschreibung:
% In diesem skript werde ich den nichtlinearen Oszillator untersuchen
% (Linearisierung, Diskretisierung, lifted Dynamics).
% -------------------------------------------------------------

clc
clear

%% System Simulation
% Simulation parameters
x0 = [0; 0];
Ts = 0.1;
t_vec = 0:Ts:50;
t_step = 10;

u_inp = zeros(size(t_vec));
u_inp(t_vec >= t_step) = 1;


% Simulation
[t_sim, x_sim] = ode45(@(t,x) oszillator_nonlinear(t, x, u_inp, Ts), t_vec, x0);
y_sim = x_sim(:, 1);

%% Lifted Dynamics
N = size(t_vec, 2);
nx = 2;

% Precompute time varying matrices A_k, B_k, C_k, D_k
Ad_seq = cell(N, 1);
Bd_seq = cell(N, 1);
Cd_seq = cell(N, 1);
Dd_seq = cell(N, 1);
for k = 1:N
    [Ad_seq{k}, Bd_seq{k}, Cd_seq{k}, Dd_seq{k}] = linear_discrete_system(x_sim(k, :), Ts);
end

% Compute lifted Matrix
P = zeros(N, N);
for k = 1:N
    % Feedthrough term (Diagonal)
    P(k, k) = Dd_seq{k};
    % Terms from past inputs
    for i = 1:(k-1)
        A_pow = eye(nx);
        % A_pow = A(k-1)*...*A(i+2)*A(i+1)
        for j = (i+1):(k-1)
            A_pow = Ad_seq{j} * A_pow;
        end
        P(k, i) = Cd_seq{k} * A_pow * Bd_seq{i};
    end
end

P1 = zeros(N, N);
for k = 1:N
    % Feedthrough term (Diagonal)
    P1(k, k) = Dd_seq{k};
    % Terms from past inputs
    A_pow = eye(nx);
    for i = (k-1):-1:1
        % Backwards
        P1(k, i) = Cd_seq{k} * A_pow * Bd_seq{i};
        A_pow = A_pow * Ad_seq{i};
    end
end

% Calculate Error
max_abs_error = max(abs(P(:) - P1(:)));
fprintf('Maximaler absoluter Unterschied bei der Bestimmung von P: %.3e\n', max_abs_error);

%% Calculate System output for changeing input
% New input
u_inp_test = zeros(size(t_vec));
u_inp_test(t_vec >= t_step) = 1.3;
delta_u_inp = u_inp_test - u_inp;

% Simulate System
[t_sim, x_sim_test] = ode45(@(t,x) oszillator_nonlinear(t, x, u_inp_test, Ts), t_vec, x0);
y_sim_test = x_sim_test(:, 1);

% Calculate system response using the lifted dynamics
y_lifted_test = y_sim + P*delta_u_inp';

% Plot results
figure;
plot(t_sim, y_sim, LineWidth=1, DisplayName='sim-star'); hold on;
plot(t_sim, y_sim_test, LineWidth=1, DisplayName='sim');
plot(t_sim, y_lifted_test, LineWidth=1, DisplayName='lifted');
grid on;
xlabel('Zeit [s]'); ylabel('x [m]');
title('Response of Oszilator to Step');
legend()

%% Local Functions
function dx = oszillator_nonlinear(t, x_vec, u_vec, Ts)
    % Simulation parameters
    m  = 2; % kg
    c1 = 2; % N/m
    c2 = 1; % N/m^3
    d  = 1; % Ns/m

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
    c2 = 1; % N/m^3
    d  = 1; % Ns/m

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