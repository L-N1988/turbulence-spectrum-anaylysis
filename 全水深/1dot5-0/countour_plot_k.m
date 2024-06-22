clc; clear; close all;

%------------------------------------------------------------------------%
% read data
%------------------------------------------------------------------------%
data = load("./denoise_concate_pxx_k.mat");
% large deviation in psd near wall
Ys = data.Ys(3:end-1);
% pxxs = data.pxxs_denoise(3:end-1);
pxxs = data.emds(3:end-1);
fs = data.fs_denoise(3:end-1);
ks = data.ks_denoise(3:end-1);

y_data = repmat(Ys, 1, length(pxxs{1}));
y_data = y_data.'; y_data = y_data(:); % convert to 1d vector
x_data = double(cell2mat(ks)); f = cell2mat(fs);
z_data = cell2mat(pxxs); z_data = z_data .* f;

assert(length(y_data) == length(x_data));
assert(length(x_data) == length(z_data));

%Create regular grid across data space
xscale = 100; yscale = 500; % increase point density
[Xmesh,Ymesh] = meshgrid(linspace(min(x_data),max(x_data),length(pxxs{1})/xscale), ...
    linspace(min(y_data),max(y_data),length(Ys)*yscale));

%create contour plot
con = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
contourf(Xmesh,Ymesh,griddata(x_data,y_data,z_data,Xmesh,Ymesh), 50, 'LineStyle','none');
colormap(jet);
col = colorbar();

set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
set(ylabel("$z(\rm m)$"), 'Interpreter', 'latex');
set(ylabel(col, "$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
saveas(con, 'PSD_contour-k.eps', 'epsc');
savefig(con, 'PSD_contour-k.fig');