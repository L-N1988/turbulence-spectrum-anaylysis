clc; clear; close all;

%% read data
% before running, add mat file path here!!!
matPath_C = '.\01h-C\';
cutedge = 2;
matPath_L = strrep(matPath_C, '-C', '-L');

data_C = load(fullfile(matPath_C, 'pxx_f-C.mat'), 'u_fluc');
data_L = load(fullfile(matPath_L, 'pxx_f-L.mat'), 'u_fluc');
uc = data_C.u_fluc;
ul = data_L.u_fluc;

% 1.5:0 转速比采样频率150Hz
Fc = 150; Lc = length(uc);
Fl = 24; Ll = length(ul);

% debug figures
% PSD(ul, Fl, "PSD filtered of original burst mode");
% PSD(uc, Fc, "PSD filtered of original continuous mode");

%% filter noise components in burst mode, using notch filter
% notch filter only distorts PSD slightly!!!
noise_f = 0.025:0.025:0.325; % 1.5:0转速比
ul_filtered = ul;
for i=1:length(noise_f)
    f_noise = noise_f(i);
        
    % Design notch filter
    wo = f_noise / (Fl / 2); % Normalize the frequency
    % the larger bandwith, the more filtered signals
    % bw = wo / 50; % Bandwidth of the notch filter (adjust as needed)
    % TODO: find proper band width for each noise component
    bw = 8e-6;
    [b, a] = iirnotch(wo, bw);

    % Apply notch filter many times, because narrow band width
    for j=1:60
        ul_filtered = filtfilt(b, a, ul_filtered);
    end
end

figure;
subplot(1, 2, 1);
[pxx, f] = pwelch(ul, [], [], [], Fl);
plot(f, pxx);
grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
title("Original Signal PSD");

subplot(1, 2, 2);
[pxx, f] = pwelch(ul_filtered, [], [], [], Fl);
plot(f, pxx);
grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
title("Filtered Signal PSD");

figure;
subplot(1, 2, 1);
[pxx, f] = pwelch(ul, [], [], [], Fl);
plot(f, f.*pxx);
grid on; 
set(gca, 'XScale', 'log');
title("Original Signal pre-multiplied PSD");

subplot(1, 2, 2);
[pxx, f] = pwelch(ul_filtered, [], [], [], Fl);
plot(f, f.*pxx);
grid on; 
set(gca, 'XScale', 'log');
title("Filtered Signal pre-multiplied PSD");

%% filter noise components in continuous mode, using notch filter
% notch filter only distorts PSD slightly!!!
noise_f = 0.025:0.025:0.325; % 1.5:0转速比
uc_filtered = uc;
for i=1:length(noise_f)
    f_noise = noise_f(i);
        
    % Design notch filter
    wo = f_noise / (Fc / 2); % Normalize the frequency
    % the larger bandwith, the more filtered signals
    % bw = wo / 50; % Bandwidth of the notch filter (adjust as needed)
    % TODO: find proper band width for each noise component
    bw = 8e-6;
    [b, a] = iirnotch(wo, bw);

    % Apply notch filter many times, because narrow band width
    for j=1:30
        uc_filtered = filtfilt(b, a, uc_filtered);
    end
end

figure;
subplot(1, 2, 1);
[pxx, f] = pwelch(uc, [], [], [], Fc);
plot(f, pxx);
grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
title("Original Signal PSD");

subplot(1, 2, 2);
[pxx, f] = pwelch(uc_filtered, [], [], [], Fc);
plot(f, pxx);
grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');
title("Filtered Signal PSD");

figure;
subplot(1, 2, 1);
[pxx, f] = pwelch(uc, [], [], [], Fc);
plot(f, f.*pxx);
grid on; 
set(gca, 'XScale', 'log');
title("Original Signal pre-multiplied PSD");

subplot(1, 2, 2);
[pxx, f] = pwelch(uc_filtered, [], [], [], Fc);
plot(f, f.*pxx);
grid on; 
set(gca, 'XScale', 'log');
title("Filtered Signal pre-multiplied PSD");

%% merge two signal series
% Target sample rate (e.g., the higher of the two)
fs_target = max(Fc, Fl);

% Resample ul to the target sample rate
[resampled_ul, Fl_resampled] = resample(ul_filtered, fs_target, Fl);

% Resample uc to the target sample rate
[resampled_uc, Fc_resampled] = resample(uc_filtered, fs_target, Fc);

% Ensure both resampled signals have the same length 
% by zero-padding the shorter one
len1 = length(resampled_uc);
len2 = length(resampled_ul);
if len1 > len2
    resampled_ul = [resampled_ul; zeros(len1 - len2, 1)];
else
    resampled_uc = [resampled_uc; zeros(len2 - len1, 1)];
end

% Compute the FFT of both resampled signals
fft_uc = fft(resampled_uc);
fft_ul = fft(resampled_ul);

% Compute the frequency bins
n = length(resampled_uc);
freq_bins = (-n/2:n/2-1) * (fs_target / n); % two-sided fft

% Compute the magnitudes of the FFT results
magnitude_l = abs(fft_ul);
magnitude_c = abs(fft_uc);

% Average the magnitudes
% % merge by interval
% average_magnitude = zeros(size(magnitude_c));
% 
% shift_index = fftshift(1:length(magnitude_c));
% interval_l = (freq_bins > -cutedge) & (freq_bins < cutedge);
% average_magnitude(shift_index(interval_l)) = magnitude_l(shift_index(interval_l));
% 
% interval_c = (freq_bins <= -cutedge) | (freq_bins >= cutedge);
% average_magnitude(shift_index(interval_c)) = magnitude_c(shift_index(interval_c));

% % merge by magnitude
% average_magnitude = (magnitude_c + magnitude_l) / 2;

% Averaging the powers instead of the magnitudes to keep PSD unchanged
average_magnitude = sqrt((magnitude_c.^2 + magnitude_l.^2) / 2);

% Convert the averaged magnitude back to the complex form
merged_fft = average_magnitude .* exp(1j * angle(fft_uc));

% Inverse FFT to get the merged signal in time domain
u_merge = ifft(merged_fft);

% Plot the results
figure;
subplot(3, 1, 1);
plot(freq_bins, abs(fftshift(fft_uc)));
title('FFT of Resampled Continuous Mode');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3, 1, 2);
plot(freq_bins, abs(fftshift(fft_ul)));
title('FFT of Resampled Burst Mode');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3, 1, 3);
plot(freq_bins, abs(fftshift(merged_fft)));
title('Merged FFT');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% check PSD result
% FIXME: high uncertainty at high frequency interval
% u_merge image part is near zero!!!
% SEE: http://matlab.izmiran.ru/help/toolbox/signal/pwelch.html#702782
% pwelch complex value has [0, fs) frequency range 
[pxx_merge, f_merge] = pwelch(u_merge, 2^16, [], [], fs_target);
% pxx_merge = pxx_merge(f_merge < 110);
% f_merge = f_merge(f_merge < 110);

% Extract the PSD for positive frequencies [0, fs/2)
index_merge = fftshift(1:length(pxx_merge));
zero_freq_index = find(index_merge == 1);
pxx_positive = pxx_merge(index_merge(zero_freq_index:end));
positive_frequencies = f_merge(1:length(pxx_positive));

% FIXME: no support theorem
% Fold the negative frequencies into positive
pxx_negative = pxx_merge(index_merge(1:zero_freq_index-1));
pxx_negative = flip(pxx_negative);
if mod(length(pxx_merge), 2) == 0
    pxx_positive = pxx_positive + pxx_negative;
else
    pxx_positive(1:end) = pxx_positive(1:end) + pxx_negative;
end

figure('Position', [10 10 1000 618]);
plot(positive_frequencies, pxx_positive);

grid on; 
% Set axes properties in a single set call
ax = gca; % Get current axes handle
set(ax, 'XScale', 'log', ... % Set X-axis to logarithmic scale
        'YScale', 'log', ... % Set Y-axis to logarithmic scale
        'FontSize', 16);     % Set font size for axes
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
set(title("PSD from Merged Signal"), 'Interpreter', 'latex');
xlim([1e-3 1e3]);

figure('Position', [10 10 1000 618]);
plot(positive_frequencies, positive_frequencies .* pxx_positive);

grid on;
% Set axes properties in a single set call
ax = gca; % Get current axes handle
set(ax, 'XScale', 'log', ... % Set X-axis to logarithmic scale
        'FontSize', 16);     % Set font size for axes
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD from Merged Signal"), 'Interpreter', 'latex');

%% split frequency domain
% very large scale motion frequency in pre-multiplied spectrum
split_f = 0.3;
shift_index_merge = fftshift(1:length(merged_fft));
shift_F_merge = (-n/2:n/2-1) * (fs_target / n);
vlsm_interval = (shift_F_merge > -split_f) & (shift_F_merge < split_f);
vlsm_fft = zeros(size(merged_fft));
vlsm_fft(shift_index_merge(vlsm_interval)) = merged_fft(shift_index_merge(vlsm_interval));
u_vlsm = ifft(vlsm_fft);

% large scale motion velocity
u_lsm = u_merge - u_vlsm;
% FIXME: what is hilbert transform about complex value?
amp_lsm = hilbert(u_lsm); % amplitude wrt LSM
amp_lsm_fft = fft(abs(amp_lsm)); % Hilbert transform uses magnitudes as envelope
shift_index_amp_lsm = fftshift(1:length(amp_lsm_fft));
% VLSM part in LSM amplitude
vlsm_in_amp_lsm_fft = zeros(size(amp_lsm_fft));
vlsm_in_amp_lsm_fft(shift_index_amp_lsm(vlsm_interval)) = amp_lsm_fft(shift_index_amp_lsm(vlsm_interval));
u_lsm_vlsm = ifft(vlsm_in_amp_lsm_fft);

% amplitude corresponding analysis between LSM and VLSM
% u_lsm_vlsm: velocity about  LSM amplitude part wrt VLSM
% u_vlsm: velcoity about VLSM
C_am = sum(u_lsm_vlsm .* u_vlsm) / sqrt(sum(u_lsm_vlsm.^2) * sum(u_vlsm.^2));

save("am_analysis.mat","split_f", "u_merge", "u_vlsm", "u_lsm", "u_lsm_vlsm", "C_am");