clc; clear; close all;
% Example signal
fs = 1000; % Sample rate (Hz)
t = 0:1/fs:1; % Time vector
f_noise = 50; % Harmonic noise frequency (Hz)
signal = sin(2*pi*100*t) + sin(2*pi*f_noise*t); % Example signal with noise
n = length(signal);
freq_bins = (-n/2:n/2-1) * (fs / n); % two-sided fft

fft_signal = fft(signal);

% Plot original signal and its PSD
figure;
subplot(3,1,1);
plot(t, signal);
title('Original Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 2)
plot(freq_bins, fftshift(fft_signal))
title('FFT of Original Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,3);
pwelch(signal, [], [], [], fs);
title('Original Signal PSD');

% Design notch filter
wo = f_noise / (fs / 2); % Normalize the frequency
bw = wo / 1; % Bandwidth of the notch filter (adjust as needed)
[b, a] = iirnotch(wo, bw);

% Apply notch filter 20 times
filtered_signal = filtfilt(b, a, signal);
% for i=1:20
%     filtered_signal = filtfilt(b, a, filtered_signal);
% end

fft_filtered_signal = fft(filtered_signal);


% Plot filtered signal and its PSD
figure;
subplot(3,1,1);
plot(t, filtered_signal);
title('Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(freq_bins, fftshift(fft_filtered_signal))
title('FFT of Filtered Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,3);
pwelch(filtered_signal, [], [], [], fs);
title('Filtered Signal PSD');
