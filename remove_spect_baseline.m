function [ newersp ] = remove_spect_baseline(tm, ersp, base_lims, type)
%REMOVE_SPECT_BASELINE remove baseline spectrum
%   tm -time vector
%   ersp - spectrogram matrix
%   base_lims - baseline limits (in the same time scale as tm vector)
%   type - type of baseline removal: minus for substraction, div for
%   division
    if nargin < 4
        type = 'minus';
    end
    greater = find(tm >= base_lims(1));
    lower = find(tm <= base_lims(2));
    base_indices = intersect(greater, lower);
    if length(base_indices) <= 1
        error(['0 or one indices found for baseline: ' num2str(base_lims) ...
            ' Check scale or increase base_lims range.'])
    end
    baseline = mean(ersp(:,base_indices),2);
    if strcmp(type, 'minus')
        newersp = ersp - repmat(baseline, [1, size(ersp,2)]);
    elseif strcmp(type, 'div')
       newersp = ersp ./ repmat(baseline, [1, size(ersp,2)]);
    else
        error('Not recognized type')
    end
end

