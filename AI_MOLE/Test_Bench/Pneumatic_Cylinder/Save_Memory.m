% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.01.2026
% Beschreibung:
% Dieses Skript dient dazu, die gespeicherten AI-MOLE objekte durchzugehen
% und die großen Matrizen zu löschen, die nichtmehr benötigt werden.
% -------------------------------------------------------------

clc
clear
close all

%% Directory Name
root_dir = fullfile(pwd, 'Runs');

%% Delete large matrices
files = dir(fullfile(root_dir, "**", "*.mat"));

% Search in files
for i = 1:numel(files)
    % Get the file paths
    file_path = fullfile(files(i).folder, files(i).name);

    % Load full file
    load(file_path);

    % Delete large Matrices
    % GP
    SISO_MOLE.GP_SISO.L_chol = [];
    SISO_MOLE.GP_SISO.V_transp = [];
    SISO_MOLE.GP_SISO.V = [];
    SISO_MOLE.GP_SISO.GP = [];

    % ILC
    SISO_MOLE.ILC_SISO.W = [];
    SISO_MOLE.ILC_SISO.S = [];
    SISO_MOLE.ILC_SISO.R = [];
    SISO_MOLE.ILC_SISO.L = [];
    SISO_MOLE.ILC_SISO.Q = [];
    SISO_MOLE.ILC_SISO.P = [];

    % Resave variables
    save(file_path, 'SISO_MOLE', 'ref_traj', 'init_update_timestamp');
end