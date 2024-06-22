function PSD(u, Fs, tile)
[pxx, f] = pwelch(u, [], [], [], Fs);
figure('Position', [10 10 1000 618]);
plot(f, pxx);

grid on; 
set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log'); 
set(gca, 'FontSize', 16);
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex'); 
set(ylabel("$S_{uu}(f) (\rm m^2/s)$"), 'Interpreter', 'latex');
set(title(tile), 'Interpreter', 'latex');