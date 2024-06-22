clc; clear; close all;

blue   = [0.0000 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];
green  = [0.4660 0.6740 0.1880];
sky    = [0.3010 0.7450 0.9330];
maroon = [0.6350 0.0780 0.1840];
gray   = [0.5020 0.5020	0.5020];

%------------------------------------------------------------------------%
% read data
%------------------------------------------------------------------------%

% before running, add mat file path here!!!
matPath_C = '../01h-C/';
cutedge = 5;

matPath_L = strrep(matPath_C, '-C', '-L');
data_C = load(fullfile(matPath_C, 'pxx_f-C-spanwise.mat'));
data_L = load(fullfile(matPath_L, 'pxx_f-L-spanwise.mat'));
Y_C = data_C.local_Ys; Y_L = data_L.local_Ys;
U_C = data_C.local_Us; U_L = data_L.local_Us;
urms_C = data_C.u_rms; urms_L = data_L.u_rms;
pxxs_C = data_C.pxxs; pxxs_L = data_L.pxxs;
fs_C = data_C.fs; fs_L = data_L.fs;
noise_f = 0.075:0.075:1; % 4.5:0转速比

concate_pxxs = cell(length(data_C.pxxs), 1);
concate_fs = cell(length(data_C.pxxs), 1);

concate_urms = (urms_L + urms_C) / 2;
concate_Ys = (Y_C + Y_L) / 2;

if ~exist("figures-k", 'dir')
   mkdir("figures-k")
end

parfor ipxx = 1:length(Y_C)
    pxx_C = pxxs_C{ipxx}; f_C = fs_C{ipxx};
    pxx_L = pxxs_L{ipxx}; f_L = fs_L{ipxx};

    concate_pxx = [pxx_L(f_L < cutedge); pxx_C(f_C >= cutedge)];
    concate_f = [f_L(f_L < cutedge); f_C(f_C >= cutedge)];
    concate_pxxs{ipxx} = concate_pxx;
    concate_fs{ipxx} = concate_f;

    % plot PSD demo
    psd_fig_CL = figure('Position', [10 10 1000 618]);
    p_C = plot(f_C, pxx_C, 'Color', blue);
    hold on
    p_L = plot(f_L, pxx_L, 'Color', orange);
    p_CL = plot(concate_f, concate_pxx, 'Color', green);

    % reference lines, noise psd
    p4 = xline(noise_f, '-.');

    grid on; 
    set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
    set(gca, 'FontSize', 16);
    set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
    set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
    figtitle = sprintf("PSD at y=%.4f m", (Y_C(ipxx) + Y_L(ipxx)) / 2);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([p_C, p_L, p_CL], {'continure', 'burst-mode', 'concated pxx'});

    epsname = sprintf("./figures-f/PSD-CL-%d.eps", ipxx);
    figname = sprintf("./figures-f/PSD-CL-%d.fig", ipxx);
    saveas(psd_fig_CL, epsname, 'epsc');
    savefig(psd_fig_CL, figname);
    close(psd_fig_CL);

    % plot pre-PSD demo
    pre_psd_CL = figure('Position', [10 10 1000 618]);
    pre_C = plot(f_C, f_C.* pxx_C, 'Color', blue);
    hold on
    pre_L = plot(f_L, f_L.* pxx_L, 'Color', orange);
    pre_CL = plot(concate_f, concate_f .* concate_pxx, 'Color', green);

    grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
    set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
    set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
    figtitle = sprintf("pre-multiplied PSD at y=%.4f m", (Y_C(ipxx) + Y_L(ipxx)) / 2);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([pre_C, pre_L, pre_CL], {'continure', 'burst-mode', 'concated pxx'});

    epsname = sprintf("./figures-f/pre-PSD-CL-%d.eps", ipxx);
    figname = sprintf("./figures-f/pre-PSD-CL-%d.fig", ipxx);
    saveas(pre_psd_CL, epsname, 'epsc');
    savefig(pre_psd_CL, figname);
    close(pre_psd_CL);
end

save("concate_pxx_f.mat", "concate_fs", "concate_pxxs", "concate_Ys", "concate_urms", "U_C", "U_L");