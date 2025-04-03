clc; clear; close all;
%-------------------------------------------------------------------------%
% % Define colors
colors.blue      = [0.3, 0.6, 1];    % Light Blue
colors.green     = [0.5, 0.9, 0.5];  % Light Green
colors.yellow    = [1, 1, 0.5];      % Light Yellow
colors.red       = [1, 0.4, 0.4];    % Optional: Add more colors
colors.purple    = [0.7, 0.5, 1];    % Optional: Extended palette
%-------------------------------------------------------------------------%
% % Define the main folders
main_folders = {
    '1dot5-0', ...
    '3-0', ...
    '4dot5-0'
};
% % Define the subfolders
subfolders = {
    '01h', ...
    '05h', ...
    '09h'
};
% % Initialize variables to store data
data_01h = {};
data_05h = {};
data_09h = {};

% Loop through each subfolder
for i = 1:length(subfolders)
    subfolder = subfolders{i};
    temp_data = struct('pxx', {}, 'f', {}, 'U', {}, 'H', {});
    
    % Loop through each main folder
    for j = 1:length(main_folders)
        main_folder = main_folders{j};
        file_path = fullfile(main_folder, subfolder, 'smooth_concate_pxx_f.mat'); % Adjust file name if needed
        
        % Load the data
        if isfile(file_path)
            loaded_data = load(file_path);
            temp_data(j).pxx = loaded_data.smooth_pxx; % Load pxx data
            temp_data(j).f = loaded_data.f_denoise; % Load frequency data
            temp_data(j).U = loaded_data.U; % Load U data
            temp_data(j).H = loaded_data.H; % Load H data
        else
            warning('File not found: %s', file_path);
        end
    end
    
    % Assign the data to the corresponding variable
    switch subfolder
        case '01h'
            data_01h = temp_data; % Store as a cell array
        case '05h'
            data_05h = temp_data; % Store as a cell array
        case '09h'
            data_09h = temp_data; % Store as a cell array
    end
end

% Now data_01h, data_05h, and data_09h contain the organized data
urms_01h = [0.007455196	0.013489619	0.01907937];
urms_05h = [0.004664566	0.006939675	0.009064086];
urms_09h = [0.008226628	0.014336447	0.018437306];

%% Plotting
% Create a figure for the PSD plot at 01h
pre_psd_01h = figure('Position', [10 10 1000 618]);
p1 = plot((data_01h(1).U ./ data_01h(1).f) / data_01h(1).H, data_01h(1).f .* data_01h(1).pxx ./ urms_01h(1)^2, 'Color', colors.green, 'LineWidth', 3);
hold on;
p2 = plot((data_01h(2).U ./ data_01h(2).f) / data_01h(2).H, data_01h(2).f .* data_01h(2).pxx ./ urms_01h(2)^2, 'Color', colors.blue, 'LineWidth', 3);
p3 = plot((data_01h(3).U ./ data_01h(3).f) / data_01h(3).H, data_01h(3).f .* data_01h(3).pxx ./ urms_01h(3)^2, 'Color', colors.red, 'LineWidth', 3);

xlim([1e-2 1e2]); % ylim([0 4e-5]);
grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex'); %set(gca, 'YScale', 'log'); 
set(xlabel("${\lambda}/{H}$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k) / \left \langle \overline{u'^2} \right \rangle$"), 'Interpreter', 'latex');
% Single legend for all plots
legend([p1, p2, p3], {'$\omega_t = 1.5 \mathrm{rpm}$', '$\omega_t = 3.0 \mathrm{rpm}$', '$\omega_t = 4.5 \mathrm{rpm}$'}, ...
    'FontSize', 22, 'Interpreter', 'latex');
hold off;
saveas(pre_psd_01h, '01h-pre-PSD-concate-lamb.eps', 'epsc');
saveas(pre_psd_01h, '01h-pre-PSD-concate-lamb.svg', 'svg');
savefig(pre_psd_01h, '01h-pre-PSD-concate-lamb.fig');

% Create a figure for the PSD plot at 05h
pre_psd_05h = figure('Position', [10 10 1000 618]);
p1 = plot((data_05h(1).U ./ data_05h(1).f) / data_05h(1).H, data_05h(1).f .* data_05h(1).pxx ./ urms_05h(1)^2, 'Color', colors.green, 'LineWidth', 3);
hold on;
p2 = plot((data_05h(2).U ./ data_05h(2).f) / data_05h(2).H, data_05h(2).f .* data_05h(2).pxx ./ urms_05h(2)^2, 'Color', colors.blue, 'LineWidth', 3);
p3 = plot((data_05h(3).U ./ data_05h(3).f) / data_05h(3).H, data_05h(3).f .* data_05h(3).pxx ./ urms_05h(3)^2, 'Color', colors.red, 'LineWidth', 3);

xlim([1e-2 1e2]); % ylim([0 4e-5]);
grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex'); %set(gca, 'YScale', 'log'); 
set(xlabel("${\lambda}/{H}$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k) / \left \langle \overline{u'^2} \right \rangle$"), 'Interpreter', 'latex');% Single legend for all plots
legend([p1, p2, p3], {'$\omega_t = 1.5 \mathrm{rpm}$', '$\omega_t = 3.0 \mathrm{rpm}$', '$\omega_t = 4.5 \mathrm{rpm}$'}, ...
    'FontSize', 22, 'Interpreter', 'latex');
hold off;
saveas(pre_psd_05h, '05h-pre-PSD-concate-lamb.eps', 'epsc');
saveas(pre_psd_05h, '05h-pre-PSD-concate-lamb.svg', 'svg');
savefig(pre_psd_05h, '05h-pre-PSD-concate-lamb.fig');

% Create a figure for the PSD plot at 09h
pre_psd_09h = figure('Position', [10 10 1000 618]);
p1 = plot((data_09h(1).U ./ data_09h(1).f) / data_09h(1).H, data_09h(1).f .* data_09h(1).pxx ./ urms_09h(1)^2, 'Color', colors.green, 'LineWidth', 3);
hold on;
p2 = plot((data_09h(2).U ./ data_09h(2).f) / data_09h(2).H, data_09h(2).f .* data_09h(2).pxx ./ urms_09h(2)^2, 'Color', colors.blue, 'LineWidth', 3);
p3 = plot((data_09h(3).U ./ data_09h(3).f) / data_09h(3).H, data_09h(3).f .* data_09h(3).pxx ./ urms_09h(3)^2, 'Color', colors.red, 'LineWidth', 3);

xlim([1e-2 1e2]); % ylim([0 4e-5]);
grid on; set(gca, 'XScale', 'log'); set(gca, 'FontSize', 22, 'TickLabelInterpreter', 'latex'); %set(gca, 'YScale', 'log'); 
set(xlabel("${\lambda}/{H}$"), 'Interpreter', 'latex'); 
set(ylabel("$kS_{uu}(k) / \left \langle \overline{u'^2} \right \rangle$"), 'Interpreter', 'latex');
% Single legend for all plots
legend([p1, p2, p3], {'$\omega_t = 1.5 \mathrm{rpm}$', '$\omega_t = 3.0 \mathrm{rpm}$', '$\omega_t = 4.5 \mathrm{rpm}$'}, ...
    'FontSize', 22, 'Interpreter', 'latex');
hold off;
saveas(pre_psd_09h, '09h-pre-PSD-concate-lamb.eps', 'epsc');
saveas(pre_psd_09h, '09h-pre-PSD-concate-lamb.svg', 'svg');
savefig(pre_psd_09h, '09h-pre-PSD-concate-lamb.fig');