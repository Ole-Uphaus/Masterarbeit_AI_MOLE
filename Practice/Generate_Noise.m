% -------------------------------------------------------------
% Autor:      Ole Uphaus
% Datum:      29.10.2025
% Beschreibung:
% In diesem skript werde ich Rauschen erzeugen und filtern.
% -------------------------------------------------------------

clc
clear

%% Noise Parameters
% Time
Ts = 0.01;
T_end = 5;
t_vec = 0:Ts:T_end;
fs = 1/Ts;

% Variance
sigma_w = 0.2;

%% Generate Noise
% White Noise
w_raw = sigma_w*randn(size(t_vec));

% Filter (IIR - Lowpass)
fc = 5;
alpha = exp(-2*pi*fc*Ts);
b1 = 1 - alpha;
a1 = [1, -alpha];
w_A = filter(b1, a1, w_raw);
w_A = w_A * (sigma_w / std(w_A));   % Bring back initial variance

% Filter (Botterworth - 1. Order Lowpass)
[bB, aB] = butter(1, fc/(fs/2), 'low');
w_B = filter(bB, aB, w_raw);
w_B = w_B * (sigma_w / std(w_B));

% Plot
figure;
hold on;
% plot(t_vec, w_raw, LineWidth=1, DisplayName='white');
plot(t_vec, w_A, LineWidth=1, DisplayName='IIR');
plot(t_vec, w_B, LineWidth=1, DisplayName='Butter');
grid on;
xlabel('Zeit [s]'); 
ylabel('Noise');
title('Compare generated Noise');
legend()