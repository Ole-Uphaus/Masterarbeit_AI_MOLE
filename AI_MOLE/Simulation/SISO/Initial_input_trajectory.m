function u_init = Initial_input_trajectory(r_vec, Ts, sigmaI)
%INITIAL_INPUT_TRAJECTORY Summary of this function goes here
%   Detailed explanation goes here

N = length(r_vec);
r_vec = r_vec(:);
Fs = 1/Ts;

% Get maximum significant Frequency (5% threshold)
f0 = max_significant_freq(r_vec, Ts);

% Design low-pass filter
Wn = f0 / (Fs/2);
order = 4;
[b, a] = butter(order, Wn, "low");

% Create white noise
u_white = sigmaI * randn(N, 1);

% Apply zero-phase Low-Pass
u_init = filtfilt(b, a, u_white);
end

