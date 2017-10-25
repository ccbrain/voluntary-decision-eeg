% Spectrograms analysis
% (look at scripts p4 and p5 before using this one)
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
how_locked = 'StimLocked';
load(['/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrum' how_locked])
DirOut = ['/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/figs/spectrograms/' how_locked '/Averaged/'];
%brain_region = 'central';
%sub_idx = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subjects = [1000,1001,1003:1014,1016:1022];

% Frontal (Fz, Fp1, Fp2, F1 –F8)           - 1, 5:14
% Central (Cz, C1-C6, T7, T8)              - 2, 15:22
% Posterior (Pz, P1-P6, Oz, O1-O2, T5, T6) - 3, 4, 23:32
electrodes.frontal = [1, 5:14];
electrodes.central = [2, 15:22];
electrodes.posterior = [3, 4, 23:32];
electrodes.all = 1:32;
regions = {'frontal', 'central', 'posterior', 'all'};
cond_fields = fields(spectrumFull);

for sub_idx = 1:length(subjects)
    for rg = 1:length(regions)
        for cond = 1: length(cond_fields)
            condition = cond_fields{cond};
            brain_region = regions{rg};
            selected_electrodes = electrodes.(brain_region);
            disp(sprintf('%s_%i_%s', condition, subjects(sub_idx), brain_region))
            spctr = squeeze(mean(spectrumFull.(condition)(sub_idx, selected_electrodes, : ,:), 2));

            fig = plot_spectrogram( tm, fq, spctr,...
                     sprintf('%s -- Subject: %i, region: %s', condition, subjects(sub_idx), brain_region));

            print(fig, [DirOut sprintf('%s_%i_%s', condition, subjects(sub_idx), brain_region)], '-dpng');
            close(fig);
        end
    end
end


%%
for cond = 1: length(cond_fields)
    condition = cond_fields{cond};
    brain_region = 'all';
    selected_electrodes = electrodes.(brain_region);

    spctr = squeeze(mean(spectrumFull.(condition)(:, selected_electrodes, : ,:), 2));
    spctr = squeeze(mean(spctr,1));
    fig = plot_spectrogram( tm, fq, spctr,...
             sprintf('GrandAvg %s region: %s', condition, brain_region));

    print(fig, ['/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/figs/spectrograms/' ...
        how_locked '/' sprintf('grandAvg_%s_%s', condition, brain_region)], '-dpng');
    close(fig);
end
