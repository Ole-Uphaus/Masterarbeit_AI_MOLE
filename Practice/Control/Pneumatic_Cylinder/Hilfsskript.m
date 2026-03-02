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
b_1 = 20; % [kg/s] 

% Berechnung weiterer Parameter
A  = pi*d^2/4;  % Kolbenfläche


