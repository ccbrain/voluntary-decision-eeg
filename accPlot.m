function accPlot(timex, datamat)
shadedErrorBar(timex ,mean(datamat,1),std(datamat,1)./sqrt(size(datamat,1)))
line([min(timex) max(timex)], [0.5 0.5], 'Color', 'r', 'LineStyle', '--');
line([0 0], [0 1], 'Color', 'r', 'LineStyle', '--');
ylim([0.9*min(mean(datamat,1)), 1.1*max(mean(datamat,1))])
xlabel('time');
ylabel('ACC');
end