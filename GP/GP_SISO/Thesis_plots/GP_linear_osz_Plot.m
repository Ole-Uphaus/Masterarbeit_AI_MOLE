% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      20.02.2026
% Beschreibung:
% In diesem skript werde ich Plots für die Masterarbeit erstellen, die
% zeigen, inwiefern der GP in der lage ist, die Systemmodelle zu
% approximieren.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

% Generate Dynamic file Paths
base_dir = fileparts(mfilename("fullpath"));
GP_Path = fullfile(base_dir, '..');
Model_Path = fullfile(base_dir, '..', '..', '..', 'System_Models');
addpath(GP_Path);
addpath(Model_Path);