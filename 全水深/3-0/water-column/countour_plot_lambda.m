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
omega_t = 3;
u_lid = 2*pi/60*3.3/2 * omega_t;

y_data = repmat(Ys, 1, length(pxxs{1}));
y_data = y_data.'; % hate weird symbols
y_data = y_data(:); % convert to 1d vector
x_data = double(cell2mat(ks)); % wave number
x_data = 2*pi ./ x_data; % wave length
f = cell2mat(fs);
z_data = cell2mat(pxxs); z_data = z_data .* f;

assert(length(y_data) == length(x_data));
assert(length(x_data) == length(z_data));

%% Create regular grid across data space
xscale = 0.008; yscale = 100; % increase point density
xgrid = logspace(log10(min(x_data)), log10(max(x_data)), length(pxxs{1})*xscale);
ygrid = linspace(min(y_data), max(y_data), length(Ys)*yscale);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
Zmesh = griddata(x_data,y_data,z_data,Xmesh,Ymesh);

%% create contour plot
% dimension wave number spectrum
con = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
contourf(Xmesh ./ H, Ymesh ./ H, Zmesh ./ u_lid^2, 8, 'LineStyle', '--');
colormap("sky");
col = colorbar();

xlim([1e-2 1e3]);
set(gca, 'XScale', 'log');  %set(gca, 'YScale', 'log'); 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex'; % latex for y-axis

% Set custom x-tick positions and labels
xtick_positions = 10.^(-2:1:3); % Tick positions
xticks(xtick_positions);
xticklabels(arrayfun(@(x) sprintf('$10^{%d}$', log10(x)), xtick_positions, 'UniformOutput', false));

% Set custom y-tick positions and labels
ytick_positions = 0.1:0.2:0.9; % Tick positions
yticks(ytick_positions);
yticklabels(arrayfun(@(y) sprintf('$%.1f$', y), ytick_positions, 'UniformOutput', false));

set(gca, 'FontSize', 22);
set(col, 'TickLabelInterpreter', 'latex');
set(xlabel("$\lambda/H$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
% set(ylabel(col, "$kS_{uu}(k) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
set(ylabel(col, "$kS_{uu}(k)/U_{lid}^2$"), 'Interpreter', 'latex');

savefig(con, 'PSD_contour-lambda.fig');
print('PSD_contour-lambda', '-dsvg', '-vector');
print('PSD_contour-lambda', '-depsc', '-vector');
print(con,'PSD_contour-lambda.jpg','-djpeg','-r500');

yline([0.1, 0.5, 0.9], '--', 'LineWidth', 3, 'Color', 'w');
savefig(con, 'PSD_contour-lambda-yline.fig');
print('PSD_contour-lambda-yline', '-dsvg', '-vector');
print('PSD_contour-lambda-yline', '-depsc', '-vector');
print(con,'PSD_contour-lambda-yline.jpg','-djpeg','-r500');
