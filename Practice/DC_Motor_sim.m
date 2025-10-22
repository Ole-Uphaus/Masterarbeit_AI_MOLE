% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      22.10.2025
% Beschreibung:
% In diesem Skript werde ich einen kleinen DC-Motor simulieren. Dazu werde
% ich das sowohl das zeitkontinuierliche System als auch das zeitdiskrete
% System nutzen. Außerdem werde ich versuchen, aus dem zeitdiskreten System
% die lifted Systemdarstellung abzuleiten (wie sie auch für ILC genutzt
% wird)
% -------------------------------------------------------------

clc
clear

%% Continuous System
% Simulation parameters
J  = 1.5e-3;     % kg m^2
d  = 1.2e-3;     % N m s/rad
R  = 2.0;        % Ohm
L  = 5.0e-3;     % H
Kt = 0.08;       % N m / A
Ke = 0.08;       % V s / rad   (oft = Kt in SI)

% Transfer function (continuously)
num_cont = [Kt];
den_cont = [L*J, L*d + R*J, R*d + Ke*Kt];
G_cont = tf(num_cont, den_cont);

% Simulation
Ts = 1e-2;
t_vec = 0:Ts:2;
u_vec = 3*ones(size(t_vec));
x_0 = zeros(1, 2);

% Simulation with ode4 (Runge-Kutta)
[y_cont, t_cont] = lsim(G_cont, u_vec, t_vec, x_0);

%% Discrete System
% Exact discretisation unsing fundamental Matrix (zero order hold)
G_disc = c2d(G_cont, Ts, 'zoh');

% Discrete "simulation"
[y_disc, t_disc] = lsim(G_disc, u_vec, t_vec);

%% Lifted System representation
% State space representation of the discrete System
[Ad, Bd, Cd, Dd] = ssdata(ss(G_disc));

N = size(t_vec, 2);
nx = size(Ad, 1);

% Markov-Parameter p
p = zeros(N, 1);
p(1) = Dd;  % Feedthrough Matrix
A_pow = eye(nx);

for k = 2:N
    p(k) = Cd * A_pow * Bd;
    A_pow = A_pow * Ad;
end

% Toeplitz-Matrix P
P = toeplitz(p, [p(1), zeros(1, N-1)]); % define first column and first row

% Observability-Matrix (system reaction to nonzero initial conditions)
O = zeros(N, nx);
A_pow = eye(nx);

for i = 1:N
    O(i,:) = Cd * A_pow;
    A_pow = A_pow * Ad;
end

% Calculation of the system response
y_disc_lifted = P*u_vec' + O*x_0';

% Compare results
abs_err = abs(y_disc - y_disc_lifted);
[max_err, kmax] = max(abs_err);

% Plot results
figure;
plot(t_cont, y_cont, LineWidth=1, DisplayName='cont'); hold on;
plot(t_disc, y_disc, LineWidth=1, DisplayName='disc');
plot(t_disc, y_disc_lifted, LineWidth=1, DisplayName='disc lifted');
grid on;
xlabel('Zeit [s]'); ylabel('\Omega [rad/s]');
title('Response of DC Motor to 3V step (lsim)');
legend()