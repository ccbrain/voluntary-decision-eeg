% Frequency band analysis - in particular BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
how_locked = 'StimLocked';
load(['/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrum' how_locked])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nrSubjects = size(spectrumOverElectrodes.Equal100,1);
conditions = fieldnames(spectrumOverElectrodes);

srate = 250;
tm = tm/1000;
time_window = 50/1000; % seconds

beta_range  = [35 55];

electrodes.frontal = [1, 5:14];
electrodes.central = [2, 15:22];
electrodes.posterior = [3, 4, 23:32];
electrodes.all = 1:32;
regions = {'frontal', 'central', 'posterior', 'all'};

time_steps = tm(1):time_window:tm(end);

brain_region = 'frontal';

clear betaband
for cc = 1:length(conditions)
    cond=conditions{cc};%betaband.(cond) = zeros(nrSubjects, length(time_steps));
    for us = 1:nrSubjects
        for tt = 1:length(time_steps)-1
            betaband.(cond)(us,tt) = power_in_frequency_and_time(tm, fq,...
                squeeze(mean(spectrumFull.(cond)(us,electrodes.(brain_region),:,:),2)), ...
                [time_steps(tt) time_steps(tt+1)], beta_range);
        end
    end
end

time_steps_ax = time_steps;
time_steps_ax = time_steps_ax(2:end);

%%plots
colors = ['m'; 'b'; 'r'; 'g'; 'y'; 'k'];

cond_names = {'Equal100'; 'Equal80'; 'Equal20'};
figure;
hold on
for c = 1:length(cond_names)
    errorbar(time_steps_ax, mean(betaband.(cond_names{c}),1), ...
        std(betaband.(cond_names{c}),1)./sqrt(nrSubjects), ['--' colors(c)]) % avg over users
end
title('Equals')
legend(cond_names)
hold off

cond_names = {'Control100'; 'Control80'; 'Control20'};
figure;
hold on
for c = 1:length(cond_names)
    errorbar(time_steps_ax, mean(betaband.(cond_names{c}),1), ...
        std(betaband.(cond_names{c}),1)./sqrt(nrSubjects), ['--' colors(c)]) % avg over users
end
title('Controls')
legend(cond_names)
hold off
