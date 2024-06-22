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
matPath_C = '.\01h-C\';
cutedge = 2;
matPath_L = strrep(matPath_C, '-C', '-L');

data_C = load(fullfile(matPath_C, 'pxx_f-C.mat'));
data_L = load(fullfile(matPath_L, 'pxx_f-L.mat'));
pxx_C = data_C.pxx; f_C = data_C.f; urms_C = data_C.u_rms; U_C = data_C.U_xt;
pxx_L = data_L.pxx; f_L = data_L.f; urms_L = data_L.u_rms; U_L = data_L.U_xt;
noise_f = 0.025:0.025:1; % 1.5:0转速比

concate_pxx = [pxx_L(f_L < cutedge); pxx_C(f_C >= cutedge)];
concate_f = [f_L(f_L < cutedge); f_C(f_C >= cutedge)];

% plot PSD demo
psd_fig_CL = figure('Position', [10 10 1000 618]);
p_C = plot(f_C, pxx_C, 'Color', blue);
hold on
p_L = plot(f_L, pxx_L, 'Color', orange);
p_CL = plot(concate_f, concate_pxx, 'Color', green);

% reference lines
p4 = xline(noise_f, '-.');

grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
set(gca, 'FontSize', 16);
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
set(title("PSD"), 'Interpreter', 'latex');
legend([p_C, p_L, p_CL], {'continure', 'burst-mode', 'concated pxx'});
saveas(psd_fig_CL, 'PSD-CL.eps', 'epsc');
savefig(psd_fig_CL, 'PSD-CL.fig');

% plot pre-PSD demo
pre_psd_CL = figure('Position', [10 10 1000 618]);
pre_C = plot(f_C, f_C.* pxx_C, 'Color', blue);
hold on
pre_L = plot(f_L, f_L.* pxx_L, 'Color', orange);
pre_CL = plot(concate_f, concate_f .* concate_pxx, 'Color', green);

grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
set(title("pre-multiplied PSD"), 'Interpreter', 'latex');
legend([pre_C, pre_L, pre_CL], {'continure', 'burst-mode', 'concated pxx'});
saveas(pre_psd_CL, 'pre-PSD-CL.eps', 'epsc');
savefig(pre_psd_CL, 'pre-PSD-CL.fig');

save("concate_pxx_f.mat", "concate_f", "concate_pxx", "urms_C", "U_C", "urms_L", "U_L");