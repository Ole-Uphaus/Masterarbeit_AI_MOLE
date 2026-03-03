%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Praktikum zur Vorlesung Nichtlineare Regelungssysteme
% 
% Versuch 2: Backstepping-Regelung eines Pneumatikzylinders
%
% Initialisierungsskript
%
% Stand: 25.11.2024
% Betreuer: Sven Weishaupt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear all; rng(0)

% Parameter des Prüfstandes (aus Praktikumsanleitung)
m  = 1.7;      % Masse des Zylinderkolbens [kg]
d  = 0.04;      % Durchmesser des Zylinders [m]
L  = 0.5;      % Länge des Zylinders [m]
V0 = 1e-5;      % Totvolumen der linken und rechten Zylinderkammer [m^3]
R  = 287.058;      % spezifische Gaskonstante für trockene Luft [J/kgK]
T0 = 293.15;      % Umgebungstemperatur [K]
n  = 1;      % Polytropen-Exponent
p0 = 1e5;      % Umgebungsdruck [Pa]
pm = 3.5e5;      % pneumatischer Mitteldruck [Pa]

% Parameter des LuGre-Reibmodells (aus der Praktikumsanleitung)
F_C = 10; % [N]
F_S = 5; % [N]
v_S = 1e-3; % [m/s]
b_0 = 50; % [kg/s]
b_1 = 20; % [kg/m] 

% Berechnung weiterer Parameter
A  = pi*d^2/4;  % Kolbenfläche

% Reglerparameter für kaskadierte Regelung des Druckes und der Position
gamma_p = 500;  % Regelung Druck           [100 .. 1000]
gamma_1 = 300;  % Regelung Position        [ 50 ..  500] 
gamma_2 = 200;  % Regelung Geschwindigkeit [ 20 .. 1000]
gamma_3 = 2000;  % Schätzung Störfehler     [500 .. 3000]

%% Solltrajektorien generieren

syms a0 a1 a2 a3 a4 a5 t

% Polynomansatz
y_d = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5;
y_dp = diff(y_d, t);
y_dpp = diff(y_dp, t);

% Gleichungen definieren
eq1 = subs(y_d, t, 0) == 0.1;
eq2 = subs(y_dp, t, 0) == 0;
eq3 = subs(y_dpp, t, 0) == 0;
eq4 = subs(y_d, t, 1) == 0.4;
eq5 = subs(y_dp, t, 1) == 0;
eq6 = subs(y_dpp, t, 1) == 0;

% Gleichungen Auflösen
solution = solve([eq1, eq2, eq3, eq4, eq5, eq6], [a0, a1, a2, a3, a4, a5]);

a0s = double(solution.a0);
a1s = double(solution.a1);
a2s = double(solution.a2);
a3s = double(solution.a3);
a4s = double(solution.a4);
a5s = double(solution.a5);

% Plotten
t_vec = linspace(0, 1, 200);
y_d_plot = a0s + a1s*t_vec + a2s*t_vec.^2 + a3s*t_vec.^3 + a4s*t_vec.^4 + a5s*t_vec.^5;

plot(t_vec, y_d_plot)

%% Entwurf Optimalregelung (LQR) des mechanischen Teilsystems
% Sample Time
Ts = 0.001;

% Zustandsraummodell
A_contin = [0 1;
    0 0];
b_contin = [0;
    1/m];
c_T = [1, 0];

% Diskretisierung
sys_contin = ss(A_contin, b_contin, c_T, 0);
sys_disc = c2d(sys_contin, Ts, 'zoh');

% Diskretes Zustandsraummodell
Ad = sys_disc.A;
bd = sys_disc.B;
c_Td = sys_disc.C;
dd = sys_disc.D;

% LQR Gewichtungsmatrizen
Q_lqr = diag([1e5, 1]);
r_lqr = 1;

% Reglerverstärkung
k_T = dlqr(Ad, bd, Q_lqr, r_lqr);

% Dynamik des geregelten Systems
Ad_cont = Ad - bd * k_T;
sys_disc_cont = ss(Ad_cont, bd, c_Td, dd, Ts);

% Statische Verstärkung berechnen
S_gain = 1 / dcgain(sys_disc_cont);

%% Choose run for Simulation
% Run
date_string = '2026_03_03';
run_filename = 'Run_01.mat';
run_filepath = fullfile(pwd, '..', '..', '..', 'AI_MOLE', 'Test_Bench', 'Pneumatic_Cylinder', 'Runs', date_string, run_filename);

% Load file
load(run_filepath)

% Get the latest u input
idx_u = find(~cellfun('isempty', SISO_MOLE.u_cell), 1, 'last');

% Extract current input signal
u_sim = SISO_MOLE.u_cell{idx_u};