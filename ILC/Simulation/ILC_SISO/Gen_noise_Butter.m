function w_filt = Gen_noise_Butter(t_vec, sigma, fc, white)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Parameters
Ts = t_vec(2) - t_vec(1);
fs = 1 / Ts;

% 1. order Butterworth
[b, a] = butter(1, fc/(fs/2));

% Calculate Noise
w_raw = randn(numel(t_vec),1);

% Apply filter if desired
if ~white
    w_filt = filter(b, a, w_raw);
else
    w_filt = w_raw;
end
w_filt = sigma * w_filt / std(w_filt);
end

