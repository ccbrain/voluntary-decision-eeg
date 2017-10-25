function fig = plot_spectrogram( T, F, S, plottitle)
%PLOT_SPECTROGRAM Plotting spectrogram S
%   T -time vector
%   F - freq vector
%   S - spectrogram matrix
%   tit - plot title
if nargin < 4
    plottitle = 0;
end
fig = figure;
hold on
h = surf(T, F, S, 'EdgeColor','none');
view([0 90])
axis tight
xlabel('Time')
ylabel('Frequency')
ylim([min(F), max(F)])

z_max = max(max(get(h,'Zdata')));
line([0 0], [min(F)*0.1 max(F)*1.1], z_max*ones(1,2), 'color', 'red',...
    'LineStyle', '--','LineWidth', 1)
colorbar;
if plottitle
    title(plottitle)
end
hold off

end

