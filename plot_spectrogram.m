function customColormap = plot_spectrogram(s, f, t, L, freq_range)
%==========================================================================
% plot_spectrogram
%
% INPUTS:
%   t           - Time vector (s)
%   f           - Frequency vector (Hz)
%   s           - Spectrogram matrix
%   L           - Window length
%   freq_range  - Frequency range [f_min, f_max] for display
%
% OUTPUT:
%   customColormap - The colormap matrix used for plotting
%==========================================================================

% Limit frequency range
idx = f >= freq_range(1) & f <= freq_range(2);
f = f(idx);
s = s(idx, :);

% Figure setup
set(gcf, 'Color', 'w');
imagesc(t, f, abs(s/L) * 2);
axis xy;
hold on;

% Custom colormap
n = 256;
red1   = linspace(1, 0, n/3)';
green1 = linspace(1, 0, n/3)';
blue1  = linspace(1, 1, n/3)';

red2   = linspace(0, 1, n/1.5)';
green2 = linspace(0, 0, n/1.5)';
blue2  = linspace(1, 0, n/1.5)';

red   = [red1; red2];
green = [green1; green2];
blue  = [blue1; blue2];
customColormap = [red, green, blue];
colormap(customColormap);

% Axes formatting
xlabel('Time(s)', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel('Frequency(Hz)', 'FontSize', 12, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 12, 'LineWidth', 1,'TickDir', 'in', 'Box', 'on');
colorbar;
colorbar;
hold off;
end