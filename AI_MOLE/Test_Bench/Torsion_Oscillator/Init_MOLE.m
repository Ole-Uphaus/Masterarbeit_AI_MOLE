% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      17.12.2025
% Beschreibung:
% Dieses Skript dient dazu, ein neuen AI-MOLE Run zu initialisieren.
% -------------------------------------------------------------

clc
clear
close all

%% Parameters
% General ILC Type ('free', 'serial', 'parallel')
ilc_type = 'free';

% Choose reference trajectory
traj_name = 'Trajectory_01.mat';

%% Load reference Trajectory
% File path
traj_path = fullfile(pwd, 'Reference_Trajectories', traj_name);

% Load Trajectory Data
load(traj_path);