function plot_erp(M,  freqsampl, timeshift, leglabels , ylimits, rev)
%PLOT_ERP plot event related potentials
% M - matrix with data (k,n), where k is number of conditions and n number
% of samples
% freqsampl - sampling frequency [in data points]
% timeshift - [in seconds] how much shift x scale from the event
% leglabels - cell with legend labels - must agree with k.
% ylimits - ylmits
if nargin < 6
   rev = 1; 
end
plot((1:size(M,2))/freqsampl - timeshift, M')
if rev
    set(gca,'Ydir','reverse')
end
ylim(ylimits)
legend(leglabels{:}, 'Location', 'sw')

end

