clc; clear; close all;

blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
sky = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840];
gray = [128, 128, 128] / 255;
%------------------------------------------------------------------------%
% read data
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

psd_fig = figure('Position', [10 10 1000 618]);
p1 = plot(f, pxx, 'Color', gray);
hold on
% plot part of spectrum
p2 = plot(f_denoise, ...
    pxx_denoise, ...
    'Color', blue);
p3 = plot(f_denoise, ...
    smoothdata(reconstrct, smooth_window{:}), ...
    'Color', yellow, LineWidth=2);
% reference lines
p4 = xline(noise_f(noise_f <= 1), '-.');
xlim([1e-3 1e3]);

grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
set(gca, 'FontSize', 16);
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
set(title("PSD"), 'Interpreter', 'latex');
legend([p1, p2, p3], {'raw', 'denoised', 'EMD'}, "FontSize", 12);
saveas(psd_fig, 'PSD-concate.eps', 'epsc');
savefig(psd_fig, 'PSD-concate.fig');

pre_psd = figure('Position', [10 10 1000 618]);
hold on
pre2 = plot(f_denoise, ...
    f_denoise .* pxx_denoise, ...
    'Color', blue);
pre3 = plot(f_denoise, ...
    f_denoise .* smoothdata(reconstrct, smooth_window{:}), ...
    'Color', yellow, LineWidth=2);
xlim([1e-3 1e3]);

grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
legend([pre2, pre3], {'denoised', 'EMD'}, "FontSize", 12);
saveas(pre_psd, 'pre-PSD-concate.eps', 'epsc');
savefig(pre_psd, 'pre-PSD-concate.fig');

% plot wave number spectrum
H = 0.15; % water depth 15cm
pre_psd_k = figure('Position', [10 10 1000 618]);
hold on
pre2 = plot(f_denoise / U * H, ...
    f_denoise .* pxx_denoise / urms^2, ...
    'Color', blue);
pre3 = plot(f_denoise / U * H, ...
    f_denoise .* smoothdata(reconstrct, smooth_window{:}) / urms^2, ...
    'Color', yellow, LineWidth=2);

grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$kH$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k)/u_{\mathrm{rms}}^2$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
legend([pre2, pre3], {'denoised', 'EMD'}, "FontSize", 12);
saveas(pre_psd_k, 'pre-PSD-concate-K.eps', 'epsc');
savefig(pre_psd_k, 'pre-PSD-concate-K.fig');