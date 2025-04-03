clc; clear; close all;
%% read data
% Define the file path
filePath = 'summary.xlsx';

% Define the sheet name and range
sheetName = 'water column turbulence';
range = 'A3:X71';

% Read the data from the specified range and sheet
data = readtable(filePath, 'Sheet', sheetName, 'Range', range);

% Define colors
colors.blue      = [0.3, 0.6, 1];    % Light Blue
colors.green     = [0.5, 0.9, 0.5];  % Light Green
colors.yellow    = [1, 1, 0.5];      % Light Yellow
colors.red       = [1, 0.4, 0.4];    % Optional: Add more colors
colors.purple    = [0.7, 0.5, 1];    % Optional: Extended palette

% Display the table to check the import
% disp(data);

% Get the number of blocks
numBlocks = size(data, 2) / 4;

% Initialize cell arrays to hold each block of data
Z = cell(1, numBlocks);
urms = cell(1, numBlocks);
vrms = cell(1, numBlocks);
rss = cell(1, numBlocks);

% Loop through each block and extract the data
for i = 1:numBlocks
    startCol = (i - 1) * 4 + 1;
    Z{i} = data{:, startCol};
    urms{i} = data{:, startCol + 1};
    vrms{i} = data{:, startCol + 2};
    rss{i} = data{:, startCol + 3};
end

% % Example: Display the first block of data
% disp('First block of data:');
% disp(table(Z{1}, urms{1}, vrms{1}, rss{1}, 'VariableNames', {'Z', 'urms', 'vrms', 'rss'}));

%% average data
% 1.5:0
urms1 = [urms{2}(1); (urms{1}(1:end-1) + urms{2}(2:end)) / 2];
vrms1 = [vrms{2}(1); (vrms{1}(1:end-1) + vrms{2}(2:end)) / 2];
rss1 = -1 * [rss{2}(1); (rss{1}(1:end-1) + rss{2}(2:end)) / 2];
z1 = [Z{2}(1); (Z{1}(1:end-1) + Z{2}(2:end)) / 2];
% 3.0:0
urms2 = [urms{4}(1:3); (urms{3}(1:end-3) + urms{4}(4:end)) / 2];
vrms2 = [vrms{4}(1:3); (vrms{3}(1:end-3) + vrms{4}(4:end)) / 2];
rss2 = -1 * [rss{4}(1:3); (rss{3}(1:end-3) + rss{4}(4:end)) / 2];
z2 = [Z{4}(1:3); (Z{3}(1:end-3) + Z{4}(4:end)) / 2];
% 4.5:0
urms3 = [urms{6}(1:2); (urms{5}(1:end-3) + urms{6}(3:end-1)) / 2];
vrms3 = [vrms{6}(1:2); (vrms{5}(1:end-3) + vrms{6}(3:end-1)) / 2];
rss3 = -1 * [rss{6}(1:2); (rss{5}(1:end-3) + rss{6}(3:end-1)) / 2];
z3 = [Z{6}(1:2); (Z{5}(1:end-3) + Z{6}(3:end-1)) / 2];

%% dimension variables
% 转速比1.5:0, 3.0:0, 4.5:0
u_half = [0.0907878395000000, 0.204672000000000, 0.306719615000000];
H = 0.15;
% velocity profile
omega_t = [1.5, 3, 4.5];
u_t = 3.3/2 * omega_t * 2 * pi / 60; % m/s
u_b = zeros(1, 3);
%% subplot
fig0 = figure('Position', [10 10 1000 618]);
h = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'loose');
% subplot(2, 1, 1)
nexttile
grid on; hold on;
scatter(urms1 ./ u_t(1), z1 / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(urms2 ./ u_t(2), z2 / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(urms3 ./ u_t(3), z3 / H, 80, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\left \langle \sqrt{\overline{u'^2}}\right \rangle/U_{lid} $"), 'Interpreter', 'latex'); 

% subplot(2, 1, 2)
nexttile
grid on; hold on;
scatter(vrms1 ./ u_t(1), z1 / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(vrms2 ./ u_t(2), z2 / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(vrms3 ./ u_t(3), z3 / H, 80, "square", 'MarkerEdgeColor', colors.red);
xlim([0.01, 0.03]);
ylim([0, 1]);

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(gca, 'YTickLabel', []);  % Remove y-axis tick labels
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
set(xlabel("$\left \langle \sqrt{\overline{w'^2}} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 

%% plot
fig1 = figure('Position', [10 10 1000 618]);
grid on; hold on;
scatter(urms1 ./ u_t(1), z1 / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(urms2 ./ u_t(2), z2 / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(urms3 ./ u_t(3), z3 / H, 80, "square", 'MarkerEdgeColor', colors.red);

ylim([0, 1]);

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\sqrt{\overline{u'^2}}/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter','latex');

%%
fig2 = figure('Position', [10 10 1000 618]);
grid on; hold on;
scatter(vrms1 ./ u_t(1), z1 / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(vrms2 ./ u_t(2), z2 / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(vrms3 ./ u_t(3), z3 / H, 80, "square", 'MarkerEdgeColor', colors.red);

ylim([0, 1]);

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\sqrt{\overline{w'^2}}/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter','latex');

%%
fig3 = figure('Position', [10 10 1000 618]);
grid on; hold on;
p7 = scatter(rss1 ./ u_t(1)^2, z1 / H, 80, "o", 'MarkerEdgeColor', colors.green);
p8 = scatter(rss2 ./ u_t(2)^2, z2 / H, 80, "^", 'MarkerEdgeColor', colors.blue);
p9 = scatter(rss3 ./ u_t(3)^2, z3 / H, 80, "square", 'MarkerEdgeColor', colors.red);
% xlim([0, 0.1]);
ylim([0, 1]);

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$\left \langle -\overline{u'w'}\right \rangle/U^2_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
% title(hLg, '$\omega_t:\omega_b$', 'Interpreter','latex');

%%
% hold off
saveas(fig0, 'uvrms.eps', 'epsc');
saveas(fig0, 'uvrms.svg', 'svg');
savefig(fig0, 'uvrms.fig');
print(fig0,'uvrms.jpg','-djpeg','-r500');

saveas(fig1, 'urms.eps', 'epsc');
saveas(fig1, 'urms.svg', 'svg');
savefig(fig1, 'urms.fig');
print(fig1,'urms.jpg','-djpeg','-r500');

saveas(fig2, 'vrms.eps', 'epsc');
saveas(fig2, 'vrms.svg', 'svg');
savefig(fig2, 'vrms.fig');
print(fig2,'vrms.jpg','-djpeg','-r500');

saveas(fig3, 'rss.eps', 'epsc');
saveas(fig3, 'rss.svg', 'svg');
savefig(fig3, 'rss.fig');
print(fig3,'rss.jpg','-djpeg','-r500');
