load /cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedDataStim.mat ;

viz = 1; subplrows = 3; subplcols = 7; % vizualization settings

subjects = [1000,1001,1003:1014,1016:1022]; 

% [2 3 4 15 16 23 24 29 30] Cz Pz Oz C1 C2 P1 P2 O1 O2
% [2 3 4] Cz Pz Oz
sel_electrodes = 1:32;%[2 3 4 15 16 23 24 29 30];

load('Data/chanlocs');

conditions = fields(FullData);

selected_conditions = {'Equal100', 'Equal80', 'Equal20'};%, 'Control100', 'Control80', 'Control20'};

subjGrandAvg = zeros(length(subjects), 32, 350);
%subjGrandAvg = zeros(length(subjects), 32, 875);
ss = 1;
for subj = subjects
    grandAvg = zeros(32,350);
    for ci = 1:length(conditions)
        cond = conditions{ci};
        if ~any(strcmp(selected_conditions, cond))
           continue 
        end
        trialAvg = squeeze(mean(FullData.(cond).(['u' num2str(subj)])(sel_electrodes,:,:),3));
        grandAvg = grandAvg + trialAvg;
    end
    subjGrandAvg(ss,:,:) = grandAvg ./ length(selected_conditions);
    ss = ss + 1;
end

for ci = 1:length(selected_conditions)
    cond = selected_conditions{ci};
    ii = 1;
    hold on
    if viz; hFig = figure; end;
    ss = 0;
    Vall = [];
    for subj = subjects
        ss = ss + 1;
        cond = conditions{ci};
        % average over conditions
        trialAvg = squeeze(mean(FullData.(cond).(['u' num2str(subj)])(sel_electrodes,:,:),3));
        % average whole data
        %trialAvg = squeeze(subjGrandAvg(ss,:,:));
        % LRP
        [U,S,V] = svd(trialAvg');
        
        for tr = 1:size(FullData.(cond).(['u' num2str(subj)]),3)
            trialData = squeeze(FullData.(cond).(['u' num2str(subj)])(sel_electrodes,:,tr));
            firstComp = trialData' * V(:,1);
            SingleERPs.(cond).(['u' num2str(subj)])(tr,:) = firstComp;
        end
        Vall = [Vall V(:,1)];
        if viz
            set(0,'CurrentFigure',hFig)
            subplot(subplrows, subplcols, ii)
            topoplot(V(:,1), chanlocs32)
            title(subj)
            ii = ii + 1;
        end
    end
    hold off
    if viz; saveas(hFig,['svd_comps_' cond '.png']); end;
    if viz
        hFig = figure;
        topoplot(abs(mean(Vall,2)), chanlocs32)
        colorbar
        saveas(hFig,['svd_comps_avg_' cond '.png']) 
    end;
end

save('Data/singleerp_stim_a','SingleERPs')

for ch=1:32
    subplot(32,1,ch)
    plot(FullData.Equal100.u1005(ch,:,10))
    ylim([-60 60])
end

tt = (1:350)/250 + epoch_length(1);
hold on
plot(tt, SingleERPs.Equal100.u1020(50,:))
xlim([-0.1,0.8])
set(gca, 'YDir','reverse')
hold off
