clc; clear; close all;
% Specify the file path
filePath = 'summary.xlsx';

% Read the entire table first to identify column names (optional)
fullData = readtable(filePath, 'Sheet', 1);

% Define the range to read (e.g., 'B2:D10' reads columns B to D and rows 2 to 10)
range = 'A3:R71';

% Read the import options for the first sheet
opts = detectImportOptions(filePath, 'Sheet', 1, 'Range', range);

% Set the variable names manually as per your headers
opts.VariableNames = {'Z_1', 'U_1', 'urms_1', 'Z_2', 'U_2', 'urms_2', ...
                      'Z_3', 'U_3', 'urms_3', 'Z_4', 'U_4', 'urms_4', ...
                      'Z_5', 'U_5', 'urms_5', 'Z_6', 'U_6', 'urms_6'};

% Read the specified range from the first sheet
data = readtable(filePath, opts);

% Define colors
colors.blue      = [0.3, 0.6, 1];    % Light Blue
colors.green     = [0.5, 0.9, 0.5];  % Light Green
colors.yellow    = [1, 1, 0.5];      % Light Yellow
colors.red       = [1, 0.4, 0.4];    % Optional: Add more colors
colors.purple    = [0.7, 0.5, 1];    % Optional: Extended palette

% velocity U at each horizontal plane
% measure height, U at 1.5rpm, 3.0rpm, 4.5rpm
u_h = [
    0.0152	0.099980439	0.210982265	0.31879817; % under 3 rotating cases
    0.0751	0.092197232	0.19693763	0.29683335;
    0.1351	0.090358094	0.20222399	0.302666185;
];

%% velocity profile
omega_t = [1.5, 3, 4.5]; % rpm
u_t = 3.3/2 * omega_t * 2 * pi / 60; % m/s
u_b = zeros(1, 3);
H = 0.15;

%% average profile
% 1.5:0
u1 = [data.U_2(1); (data.U_1(1:end-1) + data.U_2(2:end)) / 2];
z1 = [data.Z_2(1); (data.Z_1(1:end-1) + data.Z_2(2:end)) / 2];
% 3.0:0
u2 = [data.U_4(1:3); (data.U_3(1:end-3) + data.U_4(4:end)) / 2];
z2 = [data.Z_4(1:3); (data.Z_3(1:end-3) + data.Z_4(4:end)) / 2];
% 4.5:0
u3 = [data.U_6(2:3); (data.U_5(1:end-4) + data.U_6(4:end-1)) / 2];
z3 = [data.Z_6(2:3); (data.Z_5(1:end-4) + data.Z_6(4:end-1)) / 2];

% laminar flow
u1_lam = z1 * (u_t(1) - u_b(1)) / H;
u2_lam = z2 * (u_t(2) - u_b(2)) / H;
u3_lam = z3 * (u_t(3) - u_b(3)) / H;

% turbulent profile
u1_pri = u1 - u1_lam;
u2_pri = u2 - u2_lam;
u3_pri = u3 - u3_lam;

%% 3 profiles in one figure
subf = figure('Position', [10 10 1000 800]);

subplot(3, 2, 1);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on;
hold on
% scatter(u1_lam / u1(round(length(u1) / 2)), z1/H, 80, "o");
% scatter(u2_lam / u2(round(length(u2) / 2)), z2/H, 80, "^");
% scatter(u3_lam / u3(round(length(u3) / 2)), z3/H, 80, "square");
scatter(u1_lam / u_t(1), z1/H, 80, colors.green, "o");
scatter(u2_lam / u_t(2), z2/H, 80, colors.blue, "^");
scatter(u3_lam / u_t(3), z3/H, 80, colors.red, "square");
title("(a) Laminar profile", 'Interpreter', 'latex');
set(xlabel("$U_L/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hold off

subplot(3, 2, 2);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on;
hold on
% scatter(u1_pri / u1(round(length(u1) / 2)), z1/H, 80, "o");
% scatter(u2_pri / u2(round(length(u2) / 2)), z2/H, 80, "^");
% scatter(u3_pri / u3(round(length(u3) / 2)), z3/H, 80, "square");
scatter(u1_pri / u_t(1), z1/H, 80, colors.green, "o");
scatter(u2_pri / u_t(2), z2/H, 80, colors.blue, "^");
scatter(u3_pri / u_t(3), z3/H, 80, colors.red, "square");
title("(b) Deviation", 'Interpreter', 'latex');
set(xlabel("$U_T/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hold off

subplot(3, 2, [3, 4, 5, 6]);
set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
grid on;
hold on
% p1 = scatter(u1 / u1(round(length(u1) / 2)), z1 / H, 80, "o");
% p2 = scatter(u2 / u2(round(length(u2) / 2)), z2 / H, 80, "^");
% p3 = scatter(u3 / u3(round(length(u3) / 2)), z3 / H, 80, "square");
% scatter(u_h(:, 2) / u1(round(length(u1) / 2)), u_h(:, 1) / H, 100, 'blue', 'filled', "o");
% scatter(u_h(:, 3) / u2(round(length(u2) / 2)), u_h(:, 1) / H, 100, 'red', 'filled', "^");
% scatter(u_h(:, 4) / u3(round(length(u3) / 2)), u_h(:, 1) / H, 100, 'yellow', 'filled', "square");
p1 = scatter(u1 / u_t(1), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u_t(2), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u_t(3), z3 / H, 80, colors.red, "square");
scatter(u_h(:, 2) / u_t(1), u_h(:, 1) / H, 100, colors.green, 'filled', "o");
scatter(u_h(:, 3) / u_t(2), u_h(:, 1) / H, 100, colors.blue, 'filled', "^");
scatter(u_h(:, 4) / u_t(3), u_h(:, 1) / H, 100, colors.red, 'filled', "square");
xlim([0, 1]);
% Specify the position where you want to add the text
% For example, at the point (5, 0)
x_position = 0.05;
y_position = 1.15;

% Add label '(c)' at the specified position
text(x_position, y_position, '(c)', 'FontSize', 22, 'Interpreter', 'latex');
ylim([0, 1]);
set(xlabel("$U/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend([p1, p2, p3], {'1.5:0', '3.0:0', '4.5:0'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
title(hLg, '$\omega_t:\omega_b$', 'Interpreter','latex');
hold off

saveas(subf, 'averageU-decompose.eps', 'epsc');
saveas(subf, 'averageU-decompose.svg', 'svg');
savefig(subf, 'averageU-decompose.fig');
print(subf,'averageU-decompose.jpg','-djpeg','-r500');

%% plot
close all;
fig = figure('Position', [10 10 1000 618]);
grid on; hold on;
p1 = scatter(u1 / u_t(1), z1 / H, 80, colors.green, "o");
p2 = scatter(u2 / u_t(2), z2 / H, 80, colors.blue, "^");
p3 = scatter(u3 / u_t(3), z3 / H, 80, colors.red, "square");

% scatter(u_h(:, 2) / u_t(1), u_h(:, 1) / H, 100, colors.green, 'filled', "o");
% scatter(u_h(:, 3) / u_t(2), u_h(:, 1) / H, 100, colors.blue, 'filled', "^");
% scatter(u_h(:, 4) / u_t(3), u_h(:, 1) / H, 100, colors.red, 'filled', "square");

set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex');
set(xlabel("$\left \langle \overline{u} \right \rangle/U_{lid}$"), 'Interpreter', 'latex'); 
set(ylabel("$z/H$"), 'Interpreter', 'latex');
hLg = legend({'1.5 rpm', '3.0 rpm', '4.5 rpm'}, 'Location', 'best', "FontSize", 22, 'Interpreter', 'latex');
% title(hLg, '$\omega_t$', 'Interpreter','latex');
% set(hLg, 'Box', 'off');
ylim([0, 1]); % xlim([0 1]);
% set(gca, 'YScale', 'log');

% hold off
saveas(fig, 'averageU.eps', 'epsc');
saveas(fig, 'averageU.svg', 'svg');
savefig(fig, 'averageU.fig');
print(fig,'averageU.jpg','-djpeg','-r500');

%% plot
close all;
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
title(hLg, '$\omega_t:\omega_b$', 'Interpreter','latex');
xlim([0.8, 1.2]);
ylim([0, 1]);
% set(gca, 'YScale', 'log');

% hold off
saveas(fig, 'averageU-half.eps', 'epsc');
saveas(fig, 'averageU-half.svg', 'svg');
savefig(fig, 'averageU-half.fig');
print(fig,'averageU-half.jpg','-djpeg','-r500');

%% plot ustar normalization
u_star_v = [0.00345, 0.00650, 0.00946];
u_star_h = [0.00345, 0.00650, 0.00946];
close all;
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
title(hLg, '$\omega_t$', 'Interpreter','latex');
% xlim([0.8, 1.2]);
ylim([0, 1]);
% set(gca, 'YScale', 'log');
% hold off
saveas(fig, 'averageU-ustar.eps', 'epsc');
saveas(fig, 'averageU-ustar.svg', 'svg');
savefig(fig, 'averageU-ustar.fig');
print(fig,'averageU-ustar.jpg','-djpeg','-r500');