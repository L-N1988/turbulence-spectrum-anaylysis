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
matPath_C = '.\1dot5-0-C\';
cutedge = 2;

matPath_L = strrep(matPath_C, '-C', '-L');
data_C = load(fullfile(matPath_C, 'pxx_f-C.mat'));
data_L = load(fullfile(matPath_L, 'pxx_f-L.mat'));
Y_C = data_C.Y; Y_L = data_L.Y;
U_C = data_C.U_half; U_L = data_L.U_half;
noise_f = 0.025:0.025:1; % 1.5:0转速比

concate_pxxs = cell(length(data_C.pxxs), 1);
concate_fs = cell(length(data_C.pxxs), 1);
concate_ks = cell(length(data_C.pxxs), 1);
concate_Ys = zeros(length(data_C.pxxs), 1);
% y in C: 68 items, y in L:69 items
% y_C[1:end] zips with y_L[2:end]
% pxx_C[1:end] zips with y_C[2:end-1]
% pxx_L[1:end] zips with y_L[2:end-1]
% ==> pxx_C[1:end] zips with pxx_L[2:end]
% ==> f_C[1:end] zips with f_L[2:end]
for ipxx = 1:length(data_C.pxxs)
    pxx_C = data_C.pxxs{ipxx}; f_C = data_C.fs{ipxx};
    pxx_L = data_L.pxxs{ipxx+1}; f_L = data_L.fs{ipxx+1};

    concate_pxx = [pxx_L(f_L < cutedge); pxx_C(f_C >= cutedge)];
    concate_f = [f_L(f_L < cutedge); f_C(f_C >= cutedge)];
    concate_k = [f_L(f_L < cutedge) ./ U_L; f_C(f_C >= cutedge) ./ U_C];
    concate_pxxs{ipxx} = concate_pxx;
    concate_fs{ipxx} = concate_f;
    concate_ks{ipxx} = concate_k;
    concate_Ys(ipxx) = (Y_C(ipxx) + Y_L(ipxx+1)) / 2;

    % plot PSD demo
    psd_fig_CL = figure('Position', [10 10 1000 618]);
    p_C = plot(f_C ./ U_C, pxx_C, 'Color', blue);
    hold on
    p_L = plot(f_L ./ U_L, pxx_L, 'Color', orange);
    p_CL = plot(concate_k, concate_pxx, 'Color', green);

    % reference lines, noise psd
    p4 = xline(noise_f ./ U_L, '-.');

    grid on; 
    set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
    set(gca, 'FontSize', 16);
    set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
    set(ylabel("$S_{uu}(k) (\rm m^3/s^2)$"), 'Interpreter', 'latex');
    figtitle = sprintf("PSD at z=%.4f m", (Y_C(ipxx) + Y_L(ipxx+1)) / 2);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([p_C, p_L, p_CL], {'continure', 'burst-mode', 'concated pxx'});

    epsname = sprintf("./figures-k/PSD-CL-%d.eps", ipxx);
    figname = sprintf("./figures-k/PSD-CL-%d.fig", ipxx);
    saveas(psd_fig_CL, epsname, 'epsc');
    savefig(psd_fig_CL, figname);
    close(psd_fig_CL);

    % plot pre-PSD demo
    pre_psd_CL = figure('Position', [10 10 1000 618]);
    pre_C = plot(f_C ./ U_C, f_C.* pxx_C, 'Color', blue);
    hold on
    pre_L = plot(f_L ./ U_L, f_L.* pxx_L, 'Color', orange);
    pre_CL = plot(concate_k, concate_f .* concate_pxx, 'Color', green);

    grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
    set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
    set(ylabel("$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
    figtitle = sprintf("pre-multiplied PSD at z=%.4f m", (Y_C(ipxx) + Y_L(ipxx+1)) / 2);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([pre_C, pre_L, pre_CL], {'continure', 'burst-mode', 'concated pxx'});

    epsname = sprintf("./figures-k/pre-PSD-CL-%d.eps", ipxx);
    figname = sprintf("./figures-k/pre-PSD-CL-%d.fig", ipxx);
    saveas(pre_psd_CL, epsname, 'epsc');
    savefig(pre_psd_CL, figname);
    close(pre_psd_CL);
end

save("concate_pxx_k.mat", "concate_ks", "concate_fs", "concate_pxxs", "concate_Ys", "U_C", "U_L");