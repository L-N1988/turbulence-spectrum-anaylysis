clc; clear; close all;

%------------------------------------------------------------------------%
%% read data
%------------------------------------------------------------------------%
data = load("./denoise_concate_pxx_k.mat");
Ys = data.Ys(3:end);
% first 2 PSDs not convergent
pxxs = data.emds(3:end);
fs = data.fs_denoise(3:end);
ks = data.ks_denoise(3:end);
urms = data.urms; urms = urms(2:end-1);
urms = urms(3:end);
H = 0.15; % water depth 15 cm

y_data = repmat(Ys, 1, length(pxxs{1}));
y_data = y_data.'; % hate weird symbols
y_data = y_data(:); % convert to 1d vector
x_data = double(cell2mat(ks));
f = cell2mat(fs);
z_data = cell2mat(pxxs); z_data = z_data .* f;

assert(length(y_data) == length(x_data));
assert(length(x_data) == length(z_data));

%% Create regular grid across data space
xscale = 0.1; yscale = 100; % increase point density
xgrid = logspace(log10(min(x_data)), log10(max(x_data)), length(pxxs{1})*xscale);
ygrid = linspace(min(y_data), max(y_data), length(Ys)*yscale);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
Zmesh = griddata(x_data,y_data,z_data,Xmesh,Ymesh);

%% create contour plot
% dimension wave number spectrum
con = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
contourf(Xmesh, Ymesh, Zmesh, 20, 'LineStyle', 'none');
colormap(jet);
col = colorbar();

set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$k (\rm m^{-1})$"), 'Interpreter', 'latex'); 
set(ylabel("$z(\rm m)$"), 'Interpreter', 'latex');
set(ylabel(col, "$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
saveas(con, 'PSD_contour-k.eps', 'epsc');
savefig(con, 'PSD_contour-k.fig');
print(con,'PSD_contour-k.jpg','-djpeg','-r300');

%% dimensionless wave number spectrum
x_data = double(cell2mat(ks)) * H;
f = cell2mat(fs);
zzz = cellfun(@(pxx, u) pxx ./ u.^2, pxxs, num2cell(urms), 'UniformOutput', false);
z_data = cell2mat(zzz); z_data = double(z_data .* f);

%Create regular grid across data space
xgrid = xgrid*H;
% ygrid = linspace(min(y_data), max(y_data), length(Ys)*yscale);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
Zmesh = griddata(x_data,y_data,z_data,Xmesh,Ymesh);

%create contour plot
con2 = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
contourf(Xmesh, Ymesh, Zmesh, 20, 'LineStyle', 'none');
colormap(jet);
col = colorbar();

set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); %set(gca, 'YScale', 'log'); 
set(xlabel("$kH$"), 'Interpreter', 'latex'); 
set(ylabel("$z(\rm m)$"), 'Interpreter', 'latex');
set(ylabel(col, "$kS_{uu}(k)/u_{\mathrm{rms}}^2$"), 'Interpreter', 'latex');
% xlim([1e-2, 1e2]);
saveas(con2, 'PSD_contour-k-dimensionless.eps', 'epsc');
savefig(con2, 'PSD_contour-k-dimensionless.fig');
print(con2,'PSD_contour-k-dimensionless.jpg','-djpeg','-r300');
