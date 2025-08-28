clc; clear; close all;

%% Define file path and parameters
filePath = 'summary.xlsx';
H = 0.15; % Height in meters
omega_t = [1.5, 3, 4.5]; % Rotational speeds in rpm
u_t = 3.3/2 * omega_t * 2 * pi / 60; % Top velocities in m/s
u_b = zeros(1, 3); % Bottom velocities
u_half = [0.0907878395, 0.204672, 0.306719615]; % Velocities at mid-height
u_star_v = [0.00345, 0.00650, 0.00946]; % Shear velocities for normalization
u_star_h = u_star_v; % Same shear velocities for horizontal data

% Define colors
colors.blue = [0.3, 0.6, 1];    % Light Blue
colors.green = [0.5, 0.9, 0.5];  % Light Green
colors.yellow = [1, 1, 0.5];     % Light Yellow
colors.red = [1, 0.4, 0.4];      % Light Red
colors.purple = [0.7, 0.5, 1];   % Light Purple

% Velocity data at specific heights
u_h = [
    0.0152  0.099980439  0.210982265  0.31879817; % Heights and velocities at 1.5, 3.0, 4.5 rpm
    0.0751  0.092197232  0.19693763   0.29683335;
    0.1351  0.090358094  0.20222399   0.302666185;
];

%% Read velocity data (First script)
% Define range for velocity data
range_vel = 'A3:R71';
opts_vel = detectImportOptions(filePath, 'Sheet', 1, 'Range', range_vel);
opts_vel.VariableNames = {'Z_1', 'U_1', 'urms_1', 'Z_2', 'U_2', 'urms_2', ...
                          'Z_3', 'U_3', 'urms_3', 'Z_4', 'U_4', 'urms_4', ...
                          'Z_5', 'U_5', 'urms_5', 'Z_6', 'U_6', 'urms_6'};
data_vel = readtable(filePath, opts_vel);

% Compute average velocity profiles
% 1.5:0
u1 = [data_vel.U_2(1); (data_vel.U_1(1:end-1) + data_vel.U_2(2:end)) / 2];
z1 = [data_vel.Z_2(1); (data_vel.Z_1(1:end-1) + data_vel.Z_2(2:end)) / 2];
% 3.0:0
u2 = [data_vel.U_4(1:3); (data_vel.U_3(1:end-3) + data_vel.U_4(4:end)) / 2];
z2 = [data_vel.Z_4(1:3); (data_vel.Z_3(1:end-3) + data_vel.Z_4(4:end)) / 2];
% 4.5:0
u3 = [data_vel.U_6(2:3); (data_vel.U_5(1:end-4) + data_vel.U_6(4:end-1)) / 2];
z3 = [data_vel.Z_6(2:3); (data_vel.Z_5(1:end-4) + data_vel.Z_6(4:end-1)) / 2];

% Laminar flow profiles
u1_lam = z1 * (u_t(1) - u_b(1)) / H;
u2_lam = z2 * (u_t(2) - u_b(2)) / H;
u3_lam = z3 * (u_t(3) - u_b(3)) / H;

% Turbulent profiles
u1_pri = u1 - u1_lam;
u2_pri = u2 - u2_lam;
u3_pri = u3 - u3_lam;

%% Read turbulence data (Second script)
% Define sheet and range for turbulence data
sheetName = 'water column turbulence';
range_turb = 'A3:X71';
data_turb = readtable(filePath, 'Sheet', sheetName, 'Range', range_turb);

% Initialize cell arrays for turbulence data
numBlocks = size(data_turb, 2) / 4;
Z = cell(1, numBlocks);
urms = cell(1, numBlocks);
vrms = cell(1, numBlocks);
rss = cell(1, numBlocks);

% Extract turbulence data
for i = 1:numBlocks
    startCol = (i - 1) * 4 + 1;
    Z{i} = data_turb{:, startCol};
    urms{i} = data_turb{:, startCol + 1};
    vrms{i} = data_turb{:, startCol + 2};
    rss{i} = data_turb{:, startCol + 3};
end

% Compute average turbulence profiles
% 1.5:0
urms1 = [urms{2}(1); (urms{1}(1:end-1) + urms{2}(2:end)) / 2];
vrms1 = [vrms{2}(1); (vrms{1}(1:end-1) + vrms{2}(2:end)) / 2];
rss1 = -1 * [rss{2}(1); (rss{1}(1:end-1) + rss{2}(2:end)) / 2];
z1_turb = [Z{2}(1); (Z{1}(1:end-1) + Z{2}(2:end)) / 2];
% 3.0:0
urms2 = [urms{4}(1:3); (urms{3}(1:end-3) + urms{4}(4:end)) / 2];
vrms2 = [vrms{4}(1:3); (vrms{3}(1:end-3) + vrms{4}(4:end)) / 2];
rss2 = -1 * [rss{4}(1:3); (rss{3}(1:end-3) + rss{4}(4:end)) / 2];
z2_turb = [Z{4}(1:3); (Z{3}(1:end-3) + Z{4}(4:end)) / 2];
% 4.5:0
urms3 = [urms{6}(1:2); (urms{5}(1:end-3) + urms{6}(3:end-1)) / 2];
vrms3 = [vrms{6}(1:2); (vrms{5}(1:end-3) + vrms{6}(3:end-1)) / 2];
rss3 = -1 * [rss{6}(1:2); (rss{5}(1:end-3) + rss{6}(3:end-1)) / 2];
z3_turb = [Z{6}(1:2); (Z{5}(1:end-3) + Z{6}(3:end-1)) / 2];

%% Plot 1: Velocity decomposition (averageU-decompose)
subf = figure('Position', [10 10 1000 800]);
subplot(3, 2, 1);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on; hold on;
scatter(u1_lam / u_t(1), z1 / H, 80, colors.green, "o");
scatter(u2_lam / u_t(2), z2 / H, 80, colors.blue, "^");
scatter(u3_lam / u_t(3), z3 / H, 80, colors.red, "square");
title("(a) Laminar profile", 'Interpreter', 'latex');
set(xlabel("$U_L/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hold off;

subplot(3, 2, 2);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on; hold on;
scatter(u1_pri / u_t(1), z1 / H, 80, colors.green, "o");
scatter(u2_pri / u_t(2), z2 / H, 80, colors.blue, "^");
scatter(u3_pri / u_t(3), z3 / H, 80, colors.red, "square");
title("(b) Deviation", 'Interpreter', 'latex');
set(xlabel("$U_T/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hold off;

subplot(3, 2, [3, 4, 5, 6]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on; hold on;
p1 = scatter(u1 / u_t(1), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u_t(2), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u_t(3), z3 / H, 80, colors.red, "square");
scatter(u_h(:, 2) / u_t(1), u_h(:, 1) / H, 100, colors.green, 'filled', "o");
scatter(u_h(:, 3) / u_t(2), u_h(:, 1) / H, 100, colors.blue, 'filled', "^");
scatter(u_h(:, 4) / u_t(3), u_h(:, 1) / H, 100, colors.red, 'filled', "square");
xlim([0, 1]); ylim([0, 1]);
text(0.05, 1.15, '(c)', 'FontSize', 22, 'Interpreter', 'latex');
set(xlabel("$U/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend([p1, p2, p3], {'1.5:0', '3.0:0', '4.5:0'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter', 'latex');
hold off;

saveas(subf, 'averageU-decompose.eps', 'epsc');
saveas(subf, 'averageU-decompose.svg', 'svg');
savefig(subf, 'averageU-decompose.fig');
print(subf, 'averageU-decompose.jpg', '-djpeg', '-r500');

%% Plot 2: Average velocity with U_lid normalization
fig = figure('Position', [10 10 1000 618]);
grid on; hold on;
p1 = scatter(u1 / u_t(1), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u_t(2), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u_t(3), z3 / H, 80, colors.red, "square");
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$\left \langle \overline{u} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
ylim([0, 1]);
saveas(fig, 'averageU.eps', 'epsc');
saveas(fig, 'averageU.svg', 'svg');
savefig(fig, 'averageU.fig');
print(fig, 'averageU.jpg', '-djpeg', '-r500');
close(fig);

%% Plot 3: Average velocity with mid-height normalization
fig = figure('Position', [10 10 1000 618]);
grid on; hold on;
p1 = scatter(u1 / u1(round(length(u1) / 2)), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u2(round(length(u2) / 2)), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u3(round(length(u3) / 2)), z3 / H, 80, colors.red, "square");
scatter(u_h(:, 2) / u1(round(length(u1) / 2)), u_h(:, 1) / H, 100, colors.green, 'filled', "o");
scatter(u_h(:, 3) / u2(round(length(u2) / 2)), u_h(:, 1) / H, 100, colors.blue, 'filled', "^");
scatter(u_h(:, 4) / u3(round(length(u3) / 2)), u_h(:, 1) / H, 100, colors.red, 'filled', "square");
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$U/U_{0.5}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5:0', '3.0:0', '4.5:0'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter', 'latex');
xlim([0.8, 1.2]); ylim([0, 1]);
saveas(fig, 'averageU-half.eps', 'epsc');
saveas(fig, 'averageU-half.svg', 'svg');
savefig(fig, 'averageU-half.fig');
print(fig, 'averageU-half.jpg', '-djpeg', '-r500');
close(fig);

%% Plot 4: Average velocity with u* normalization
fig = figure('Position', [10 10 1000 618]);
grid on; hold on;
p1 = scatter(u1 / u_star_v(1), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u_star_v(2), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u_star_v(3), z3 / H, 80, colors.red, "square");
scatter(u_h(:, 2) / u_star_h(1), u_h(:, 1) / H, 100, colors.green, 'filled', "o");
scatter(u_h(:, 3) / u_star_h(2), u_h(:, 1) / H, 100, colors.blue, 'filled', "^");
scatter(u_h(:, 4) / u_star_h(3), u_h(:, 1) / H, 100, colors.red, 'filled', "square");
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$U/u^*$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5', '3.0', '4.5'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t$', 'Interpreter', 'latex');
ylim([0, 1]);
saveas(fig, 'averageU-ustar.eps', 'epsc');
saveas(fig, 'averageU-ustar.svg', 'svg');
savefig(fig, 'averageU-ustar.fig');
print(fig, 'averageU-ustar.jpg', '-djpeg', '-r500');
close(fig);

%% Plot 5: Turbulence profiles (urms and vrms)
fig0 = figure('Position', [10 10 1000 618]);
h = tiledlayout(1, 2, 'TileSpacing', 'tight', 'Padding', 'loose');
nexttile;
grid on; hold on;
scatter(urms1 ./ u_t(1), z1_turb / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(urms2 ./ u_t(2), z2_turb / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(urms3 ./ u_t(3), z3_turb / H, 80, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\left \langle \sqrt{\overline{u'^2}}\right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
hold off;

nexttile;
grid on; hold on;
scatter(vrms1 ./ u_t(1), z1_turb / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(vrms2 ./ u_t(2), z2_turb / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(vrms3 ./ u_t(3), z3_turb / H, 80, "square", 'MarkerEdgeColor', colors.red);
xlim([0.01, 0.03]); ylim([0, 1]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(gca, 'YTickLabel', []);
set(xlabel("$\left \langle \sqrt{\overline{w'^2}} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
hold off;

saveas(fig0, 'uvrms.eps', 'epsc');
saveas(fig0, 'uvrms.svg', 'svg');
savefig(fig0, 'uvrms.fig');
print(fig0, 'uvrms.jpg', '-djpeg', '-r500');
close(fig0);

%% Plot 6: urms profile
fig1 = figure('Position', [10 10 1000 618]);
grid on; hold on;
scatter(urms1 ./ u_t(1), z1_turb / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(urms2 ./ u_t(2), z2_turb / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(urms3 ./ u_t(3), z3_turb / H, 80, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\sqrt{\overline{u'^2}}/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter', 'latex');
saveas(fig1, 'urms.eps', 'epsc');
saveas(fig1, 'urms.svg', 'svg');
savefig(fig1, 'urms.fig');
print(fig1, 'urms.jpg', '-djpeg', '-r500');
close(fig1);

%% Plot 7: vrms profile
fig2 = figure('Position', [10 10 1000 618]);
grid on; hold on;
scatter(vrms1 ./ u_t(1), z1_turb / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(vrms2 ./ u_t(2), z2_turb / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(vrms3 ./ u_t(3), z3_turb / H, 80, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\sqrt{\overline{w'^2}}/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter', 'latex');
saveas(fig2, 'vrms.eps', 'epsc');
saveas(fig2, 'vrms.svg', 'svg');
savefig(fig2, 'vrms.fig');
print(fig2, 'vrms.jpg', '-djpeg', '-r500');
close(fig2);

%% Plot 8: rss profile
fig3 = figure('Position', [10 10 1000 618]);
grid on; hold on;
scatter(rss1 ./ u_t(1)^2, z1_turb / H, 80, "o", 'MarkerEdgeColor', colors.green);
scatter(rss2 ./ u_t(2)^2, z2_turb / H, 80, "^", 'MarkerEdgeColor', colors.blue);
scatter(rss3 ./ u_t(3)^2, z3_turb / H, 80, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$\left \langle -\overline{u'w'}\right \rangle/U^2_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
saveas(fig3, 'rss.eps', 'epsc');
saveas(fig3, 'rss.svg', 'svg');
savefig(fig3, 'rss.fig');
print(fig3, 'rss.jpg', '-djpeg', '-r500');
close(fig3);

%% Plot 9: avg U, rss, uvrms profile
fig4 = figure('Position', [10 10 3000 618*2]);
h = tiledlayout(1, 4, 'TileSpacing', 'tight', 'Padding', 'loose');

nexttile;
grid on; hold on;
scatter(u1 ./ u_t(1), z1 / H, 280, "o", 'MarkerEdgeColor', colors.green);
scatter(u2 ./ u_t(2), z2 / H, 280, "^", 'MarkerEdgeColor', colors.blue);
scatter(u3 ./ u_t(3), z3 / H, 280, "square", 'MarkerEdgeColor', colors.red);
set(gca, 'FontSize', 32, 'TickLabelInterpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
set(xlabel("$\left \langle \overline{u} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
ylim([0, 1]);
hold off;

nexttile;
grid on; hold on;
scatter(rss1 ./ u_t(1)^2, z1_turb / H, 280, "o", 'MarkerEdgeColor', colors.green);
scatter(rss2 ./ u_t(2)^2, z2_turb / H, 280, "^", 'MarkerEdgeColor', colors.blue);
scatter(rss3 ./ u_t(3)^2, z3_turb / H, 280, "square", 'MarkerEdgeColor', colors.red);
set(gca, 'FontSize', 32, 'TickLabelInterpreter', 'latex');
set(gca, 'YTickLabel', []);
set(xlabel("$\left \langle -\overline{u'w'}\right \rangle/U^2_{lid}$"), 'Interpreter', 'latex');
xlim([-1e-4, 4.5e-4]);
ylim([0, 1]);
hold off;

nexttile;
grid on; hold on;
scatter(urms1 ./ u_t(1), z1_turb / H, 280, "o", 'MarkerEdgeColor', colors.green);
scatter(urms2 ./ u_t(2), z2_turb / H, 280, "^", 'MarkerEdgeColor', colors.blue);
scatter(urms3 ./ u_t(3), z3_turb / H, 280, "square", 'MarkerEdgeColor', colors.red);
ylim([0, 1]);
set(gca, 'FontSize', 32, 'TickLabelInterpreter', 'latex');
set(gca, 'YTickLabel', []);
set(xlabel("$\left \langle \sqrt{\overline{u'^2}}\right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
ylim([0, 1]);
hold off;

nexttile;
grid on; hold on;
scatter(vrms1 ./ u_t(1), z1_turb / H, 280, "o", 'MarkerEdgeColor', colors.green);
scatter(vrms2 ./ u_t(2), z2_turb / H, 280, "^", 'MarkerEdgeColor', colors.blue);
scatter(vrms3 ./ u_t(3), z3_turb / H, 280, "square", 'MarkerEdgeColor', colors.red);
xlim([0.01, 0.03]); ylim([0, 1]);
set(gca, 'FontSize', 32, 'TickLabelInterpreter', 'latex');
set(gca, 'YTickLabel', []);
set(xlabel("$\left \langle \sqrt{\overline{w'^2}} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 32, 'Interpreter', 'latex');
% Increase legend marker size
lgd_objs = findobj(hLg, 'Type', 'Line'); % Find line/marker objects in legend
for i = 1:length(lgd_objs)
    lgd_objs(i).MarkerSize = 280; % Set larger marker size
end
ylim([0, 1]);
hold off;
saveas(fig4, 'avg_u_rss_uvrms.svg', 'svg');
saveas(fig4, 'avg_u_rss_uvrms.png', 'png');
saveas(fig4, 'avg_u_rss_uvrms.eps', 'epsc');