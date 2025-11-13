function f0 = max_significant_freq(r_vec, Ts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = length(r_vec);
r_vec = r_vec(:);
Fs = 1/Ts;
plot_fft = false;

% Remove DC-part (zero mean) --> no peak at 0 Hz
r_vec_ZM = r_vec - mean(r_vec);

% FFT
R = fft(r_vec_ZM);

% Single-sided (SS) spectrum
R_abs_SS = abs(R(1:floor(N/2)+1));      % Frequencies reach from 0, ..., Fn (Nyquist)
delta_f = Fs/N;
f_vec = (0:floor(N/2))' * delta_f;

% Find largest "significant" frequency use a threshold to decide, which
% frequency might be significant
threshold = 0.05 * max(R_abs_SS);   % 5% of maximum amplidude
idx = find(R_abs_SS > threshold);
if isempty(idx)
    % This case will only apper if we have constant signals
    [~, idx_max] = max(R_abs_SS);
    f0 = f(idx_max);
else
    f0 = f_vec(max(idx));
end

% Plot FFT
if plot_fft
    figure;
    plot(f_vec, R_abs_SS, LineWidth=1, DisplayName='FFT'); hold on;
    plot(f_vec, threshold*ones(size(f_vec)), LineWidth=1, Color='r', DisplayName='threshold');
    plot([f0, f0], [0, max(R_abs_SS)], LineWidth=1, Color='g', DisplayName='f0');
    grid on;
    xlabel('Frequency [Hz]'); 
    ylabel('x [m]');
    title('FFT of reference Trajectory');
    legend('Location', 'best');
end
end

