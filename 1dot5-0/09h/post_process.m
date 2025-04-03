clc; clear; close all;

blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250 0.8];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
sky = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840 0.8];
gray = [128, 128, 128] / 255;

colors.blue      = [0.3, 0.6, 1];    % Light Blue
colors.green     = [0.5, 0.9, 0.5];  % Light Green
colors.yellow    = [1, 1, 0.5];      % Light Yellow
colors.red       = [1, 0.4, 0.4];    % Optional: Add more colors
colors.purple    = [0.7, 0.5, 1];    % Optional: Extended palette
%------------------------------------------------------------------------%
%% read data
%------------------------------------------------------------------------%

% before running, add mat file path here!!!
matPath = '.\concate_pxx_f.mat';
data = load(matPath);
pxx = data.concate_pxx;
f = data.concate_f;
urms_C = data.urms_C; U_C = data.U_C;
urms_L = data.urms_L; U_L = data.U_L;
assert(length(urms_L) == length(urms_C));
ref_pt = floor(length(urms_L) / 2) + 1;
urms = (urms_L + urms_C) / 2;
urms = urms(ref_pt);
U = (U_C + U_L) / 2;
U = U(ref_pt);
H = 0.15; % water depth 15cm

flim = 1;
% basic arguments
neps = 3;
nmode = 8;
smooth_window = {"gaussian", 200};
noise_f = 0.025:0.025:flim; % 1.5:0转速比

% eps = max(abs(diff(f(f<flim)))) * neps; % threshold is 10
valid = zeros(size(f));
for i = 2:length(f)
    eps = abs(f(i) - f(i-1)) * neps;
    valid(i) = min(abs(noise_f - f(i))) > eps; % tolerance
end
valid = logical(valid);
f_denoise = f(valid);
pxx_denoise = pxx(valid);

% EMD smooth data
[imf,residual] = emd(pxx_denoise);
reconstrct = sum(imf(:, end-nmode:end), 2) + residual;
smooth_pxx = smoothdata(reconstrct, smooth_window{:});

% Metadata structure
metadata = struct();
metadata.fileDate = datetime('now');  % Current date and time
metadata.help = table(...
    {'f_denoise'; 'pxx_denoise'; 'smooth_pxx'; 'U'; 'H'}, ...  % Field names
    {'concatenate frequency'; 'concatenate power spectral density'; 'smooth concatenate power spectral density'; 'local convection velocity'; 'water depth'}, ...  % Descriptions
    'VariableNames', {'Field', 'Description'});  % Column headers;

save('smooth_concate_pxx_f', 'f_denoise', 'pxx_denoise', 'smooth_pxx', 'U', 'H', 'metadata');


%% PSD
psd_fig = figure('Position', [10 10 1000 618]);
p1 = plot(f, pxx, 'Color', gray);
hold on
% plot part of spectrum
p2 = plot(f_denoise, ...
    pxx_denoise, ...
    'Color', blue);
p3 = plot(f_denoise, ...
    smooth_pxx, ...
    'Color', colors.green, LineWidth=3);
% reference lines
p4 = xline(noise_f(noise_f <= 1), '-.');
xlim([1e-3 1e3]); ylim([1e-10 1e-2]);

grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
% set(title("PSD"), 'Interpreter', 'latex');
legend([p1, p2, p3], {'raw $S_{uu}$', 'denoised $S_{uu}$', 'fitted line'}, "FontSize", 22, 'Interpreter', 'latex');
% legend([p1, p2], {'raw $S_{uu}$', 'denoised $S_{uu}$'}, "FontSize", 22, 'Interpreter', 'latex');
saveas(psd_fig, 'PSD-concate.eps', 'epsc');
saveas(psd_fig, 'PSD-concate.svg', 'svg');
savefig(psd_fig, 'PSD-concate.fig');

%% plot wave length spectrum
pre_psd_lamb = figure('Position', [10 10 1000 618]);
hold on
pre2_lamb = plot((U ./ f_denoise) / H, ...
    f_denoise .* pxx_denoise, ...
    'Color', blue);
pre3_lamb = plot((U ./ f_denoise) / H, ...
    f_denoise .* smooth_pxx, ...
    'Color', colors.green, LineWidth=3);

xlim([1e-2 1e3]);
grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex'); %set(gca, 'YScale', 'log'); 
set(xlabel("${\lambda}/{H}$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
% set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
% legend([pre2, pre3], {'denoised', 'EMD'}, "FontSize", 22);
legend(pre3_lamb, {'fitted line'}, "FontSize", 22, 'Interpreter', 'latex');
saveas(pre_psd_lamb, 'pre-PSD-concate-lamb.eps', 'epsc');
saveas(pre_psd_lamb, 'pre-PSD-concate-lamb.svg', 'svg');
savefig(pre_psd_lamb, 'pre-PSD-concate-lamb.fig');

%% pre-multiplied spectrum
pre_psd = figure('Position', [10 10 1000 618]);
hold on
pre2 = plot(f_denoise, ...
    f_denoise .* pxx_denoise, ...
    'Color', blue);
pre3 = plot(f_denoise, ...
    f_denoise .* smooth_pxx, ...
    'Color', colors.green, LineWidth=3);
xlim([1e-3 1e3]);

grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22); %set(gca, 'YScale', 'log'); 
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
legend([pre2, pre3], {'denoised', 'fitted line'}, "FontSize", 22);
saveas(pre_psd, 'pre-PSD-concate.eps', 'epsc');
saveas(pre_psd, 'pre-PSD-concate.svg', 'svg');
% saveas(pre_psd, 'pre-PSD-concate.emf', 'meta');
% Save the figure as an EMF file with greater than 300 DPI resolution
% print(pre_psd, 'pre-PSD-concate.emf', '-dmeta', '-r400');
savefig(pre_psd, 'pre-PSD-concate.fig');

%% plot wave number spectrum
H = 0.15; % water depth 15cm
pre_psd_k = figure('Position', [10 10 1000 618]);
hold on
pre2 = plot(f_denoise / U * H, ...
    f_denoise .* pxx_denoise / urms^2, ...
    'Color', blue);
pre3 = plot(f_denoise / U * H, ...
    f_denoise .* smooth_pxx / urms^2, ...
    'Color', colors.green, LineWidth=3);

grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22); %set(gca, 'YScale', 'log'); 
set(xlabel("$kH$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k)/u_{\mathrm{rms}}^2$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
legend([pre2, pre3], {'denoised', 'EMD'}, "FontSize", 22);
saveas(pre_psd_k, 'pre-PSD-concate-K.eps', 'epsc');
saveas(pre_psd_k, 'pre-PSD-concate-K.svg', 'svg');
% saveas(pre_psd_k, 'pre-PSD-concate-K.emf', 'meta');
% print(pre_psd_k, 'pre-PSD-concate-K.emf', '-dmeta', '-r400');
savefig(pre_psd_k, 'pre-PSD-concate-K.fig');