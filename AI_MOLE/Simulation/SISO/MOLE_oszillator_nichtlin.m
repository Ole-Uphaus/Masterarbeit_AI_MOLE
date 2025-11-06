% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      06.11.2025
% Beschreibung:
% In diesem skript werde ich AI-MOLE für den nichtlinearen Oszillator
% erproben.
% -------------------------------------------------------------

clc
clear
close all
rng(43);

%% Initialize AI-MOLE
% Initialisation
SISO_MOLE = SISO_MOLE_IO();