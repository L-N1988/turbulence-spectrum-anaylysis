% Example synthetic PIV data (replace with real data)
clc; clear; close all;
[X, Y] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
% Synthetic turbulence example (add random fluctuations to sinusoidal base)
u = sin(2*pi*X) .* cos(2*pi*Y) + 0.5 * randn(size(X));  % Add Gaussian noise
v = -cos(2*pi*X) .* sin(2*pi*Y) + 0.5 * randn(size(X)); % Add Gaussian noise
dx = X(1,2) - X(1,1);            % Grid spacing in x
dy = Y(2,1) - Y(1,1);            % Grid spacing in y
filter_size = 2.0;                % Gaussian filter width (sigma)

Pi = compute_energy_flux(u, v, dx, dy, filter_size);

% --- Step 1: Compute Energy Flux at Different filter_size ---
filter_size = 1:0.5:100;
Pi_center = zeros(1,length(filter_size));
Pi_quarter = zeros(1,length(filter_size));
eps_center = zeros(1,length(filter_size));
eps_quarter = zeros(1,length(filter_size));
for ii = 1:length(filter_size)
    [tmp_Pi, tmp_eps] = compute_energy_flux(u, v, dx, dy, filter_size(ii));
    Pi_center(ii) = tmp_Pi(50,50);
    Pi_quarter(ii) = tmp_Pi(25,25);
    eps_center(ii) = tmp_eps(50,1);
    eps_quarter(ii) = tmp_eps(25,1);
end

% --- Step 2: Plot Results ---
figure;
plot(filter_size, Pi_center, 'Linestyle', '--', 'Linewidth', 2);
hold on;
plot(filter_size, Pi_quarter, 'Linestyle', '--', 'LineWidth', 2);
plot(filter_size, eps_center, 'Linestyle', '-.', 'LineWidth', 2);
plot(filter_size, eps_quarter, 'Linestyle', '-.', 'LineWidth', 2);
legend('Pi center', 'Pi quarter', 'eps center', 'eps quarter');
xlabel('Filter Size');
ylabel('Energy Flux at (50,50) and (25,25)');
title('Energy Flux vs. Filter Size');

% --- Step 2: Plot Results ---
figure;

% Subplot 1: U-velocity
subplot(2, 2, 1);
imagesc(u);
colorbar;
title('U-Velocity');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis equal tight;

% Subplot 2: V-velocity
subplot(2, 2, 2);
imagesc(v);
colorbar;
title('V-Velocity');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis equal tight;

% Subplot 3: SGS Energy Flux (Pi)
subplot(2, 2, 3);
imagesc(Pi);
colorbar;
title('SGS Energy Flux (\Pi)');
xlabel('x (pixels)');
ylabel('y (pixels)');
axis equal tight;

% Subplot 4: Histogram of Pi (optional)
subplot(2, 2, 4);
histogram(Pi(:), 50, 'Normalization', 'pdf');
title('PDF of \Pi');
xlabel('\Pi');
ylabel('Probability');
grid on;

