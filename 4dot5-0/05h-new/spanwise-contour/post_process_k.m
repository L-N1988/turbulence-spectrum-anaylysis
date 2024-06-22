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
matPath = '.\concate_pxx_k.mat';
data = load(matPath);
pxxs_denoise = cell(length(data.concate_pxxs), 1);
fs_denoise = cell(length(data.concate_pxxs), 1);
ks_denoise = cell(length(data.concate_pxxs), 1);
emds = cell(length(data.concate_pxxs), 1);

Ys = data.concate_Ys;
U_C = data.U_C; U_L = data.U_L;
urms = data.concate_urms;
concate_pxxs = data.concate_pxxs;
concate_fs = data.concate_fs;
concate_ks = data.concate_ks;
if ~exist("concate_figures-k", 'dir')
   mkdir("concate_figures-k")
end
parfor ipxx = 1:length(data.concate_pxxs)
    pxx = concate_pxxs{ipxx};
    f = concate_fs{ipxx};
    k = concate_ks{ipxx};
    Y = Ys(ipxx);

    neps = 10; nmode = 8;
    smooth_window = {"gaussian", 200};
    flim = 3; noise_f = 0.075:0.075:flim; % 4.5:0转速比

    valid = zeros(size(f));
    for i = 2:length(f)
        eps = abs(f(i) - f(i-1)) * neps;
        valid(i) = min(abs(noise_f - f(i))) > eps; % tolerance
    end
    valid = logical(valid);
    f_denoise = f(valid);
    k_denoise = k(valid);
    pxx_denoise = pxx(valid);
    % store denoised values
    fs_denoise{ipxx} = f_denoise;
    ks_denoise{ipxx} = k_denoise;
    pxxs_denoise{ipxx} = pxx_denoise;

    % % EMD smooth data
    [imf,residual] = emd(pxx_denoise);
    if nmode < size(imf, 2)
        reconstrct = sum(imf(:, end-nmode:end), 2) + residual;
    else
        reconstrct = sum(imf(:, 2:end), 2) + residual;
    end
    

    psd_fig = figure('Position', [10 10 1000 618]);
    p1 = plot(k, pxx, 'Color', gray);
    hold on
    % plot part of spectrum
    p2 = plot(k_denoise, pxx_denoise, 'Color', blue);
    p3 = plot(k_denoise, ...
        smoothdata(reconstrct, smooth_window{:}), ...
        'Color', yellow, LineWidth=2);
    % reference lines
    p4 = xline(noise_f ./ U_L(ipxx), '-.');
    % ylim([0.7E-8, 1E-2]);

    grid on; 
    set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
    set(gca, 'FontSize', 16);
    set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
    set(ylabel("$S_{uu}(k) (\rm m^3/s^2)$"), 'Interpreter', 'latex');
    figtitle = sprintf("PSD at y=%.4f m", Y);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([p1, p2, p3], {'raw', 'denoised', 'EMD'}, "FontSize", 12);
    
    epsname = sprintf("./concate_figures-k/PSD-concate-%d.eps", ipxx);
    figname = sprintf("./concate_figures-k/PSD-concate-%d.fig", ipxx);
    saveas(psd_fig, epsname, 'epsc');
    savefig(psd_fig, figname);
    close(psd_fig);

    pre_psd = figure('Position', [10 10 1000 618]);
    hold on
    pre2 = plot(k_denoise, f_denoise .* pxx_denoise, 'Color', blue);
    pre3 = plot(k_denoise, ...
        f_denoise .* smoothdata(reconstrct, smooth_window{:}), ...
        'Color', yellow, LineWidth=2);
    emds{ipxx} = smoothdata(reconstrct, smooth_window{:});
    % ylim([0, 5.5E-5]);

    grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
    set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
    set(ylabel("$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
    figtitle = sprintf("pre-multiplied PSD at y=%.4f m", Y);
    set(title(figtitle), 'Interpreter', 'latex');
    legend([pre2, pre3], {'denoised', 'EMD'}, "FontSize", 12);

    epsname = sprintf("./concate_figures-k/pre-PSD-concate-%d.eps", ipxx);
    figname = sprintf("./concate_figures-k/pre-PSD-concate-%d.fig", ipxx);
    saveas(pre_psd, epsname, 'epsc');
    savefig(pre_psd, figname);
    close(pre_psd);
end

save("denoise_concate_pxx_k.mat", "ks_denoise", "fs_denoise", "pxxs_denoise", "Ys", "emds", "U_C", "U_L", "urms");