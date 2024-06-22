%% test resample
fs = 10;
t1 = 0:1/fs:1;
x = t1;
y = resample(x,3,2);
t2 = (0:(length(y)-1))*2/(3*fs);

figure()
plot(t1,x,'*',t2,y,'o')
xlabel('Time (s)')
ylabel('Signal')
legend('Original','Resampled', ...
    'Location','NorthWest')

%% merge two signals
% Original signals and sample rates
signal1 = randn(1, 1000);  % Example signal 1
signal2 = randn(1, 2000);  % Example signal 2
fs1 = 1000;  % Sample rate of signal1
fs2 = 2000;  % Sample rate of signal2

% Target sample rate (e.g., the higher of the two)
fs_target = max(fs1, fs2);

% Resample signal1 to the target sample rate
[resampled_signal1, fs1_resampled] = resample(signal1, fs_target, fs1);

% Resample signal2 to the target sample rate
[resampled_signal2, fs2_resampled] = resample(signal2, fs_target, fs2);

% Compute the FFT of both resampled signals
fft_signal1 = fft(resampled_signal1);
fft_signal2 = fft(resampled_signal2);

% Compute the frequency bins
n = length(resampled_signal1);  % Assuming both resampled signals have the same length
freq_bins = (0:n-1) * (fs_target / n);

% Compute the magnitudes of the FFT results
magnitude1 = abs(fft_signal1);
magnitude2 = abs(fft_signal2);

% Average the magnitudes
average_magnitude = (magnitude1 + magnitude2) / 2;

% Convert the averaged magnitude back to the complex form
merged_fft = average_magnitude .* exp(1j * angle(fft_signal1));

% Inverse FFT to get the merged signal in time domain
merged_signal = ifft(merged_fft);

% Plot the results
figure;
subplot(3, 1, 1);
plot(freq_bins, abs(fft_signal1));
title('FFT of Resampled Signal 1');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3, 1, 2);
plot(freq_bins, abs(fft_signal2));
title('FFT of Resampled Signal 2');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3, 1, 3);
plot(freq_bins, abs(merged_fft));
title('Merged FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
