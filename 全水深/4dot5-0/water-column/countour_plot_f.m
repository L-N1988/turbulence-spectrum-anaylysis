clc; clear; close all;

%------------------------------------------------------------------------%
%% read data
%------------------------------------------------------------------------%
data = load("./denoise_concate_pxx_f.mat");
% large deviation in psd near wall
Ys = data.Ys(3:end);
pxxs = data.emds(3:end);
fs = data.fs_denoise(3:end);

y_data = repmat(Ys, 1, length(pxxs{1}));
y_data = y_data.'; y_data = y_data(:); % convert to 1d vector
x_data = cell2mat(fs);
z_data = cell2mat(pxxs); z_data = z_data .* x_data;

assert(length(y_data) == length(x_data));
assert(length(x_data) == length(z_data));

%Create regular grid across data space
xscale = 0.1; yscale = 100; % increase point density
xgrid = logspace(log10(min(x_data)), log10(max(x_data)), length(pxxs{1})*xscale);
ygrid = linspace(min(y_data), max(y_data), length(Ys)*yscale);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
Zmesh = griddata(x_data,y_data,z_data,Xmesh,Ymesh);

%% create contour plot
con = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
contourf(Xmesh,Ymesh, Zmesh, 20, 'LineStyle','none');
colormap(jet);
col = colorbar();

set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$z(\rm m)$"), 'Interpreter', 'latex');
set(ylabel(col,"$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex')

%% save figure
saveas(con, 'PSD_contour-f.eps', 'epsc');
savefig(con, 'PSD_contour-f.fig');
print(con,'PSD_contour-f.jpg','-djpeg','-r300');