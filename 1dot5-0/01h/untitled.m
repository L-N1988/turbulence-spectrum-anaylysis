% Sample data
x = linspace(0, 10, 100);
y = sin(x);

% Plot the data
figure;
plot(x, y);
xlabel('X-axis');
ylabel('Y-axis');
title('Equal Spacing of Axis Tick Labels');

% Ensure the tick label spacing is the same on both axes
axis equal;

% Optionally, adjust the aspect ratio of the plot if needed
pbaspect([1 1 1]);

% Adjust the axes properties to ensure equal tick spacing
ax = gca;
ax.DataAspectRatio = [1 1 1];

% Display the figure
grid on;
