% Script for epoch cutting and Response Locked Analysis
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleaning_flag = 0;
% Frontal (Fz, Fp1, Fp2, F1 –F8)           - 1, 5:14
% Central (Cz, C1-C6, T7, T8)              - 2, 15:22
% Posterior (Pz, P1-P6, Oz, O1-O2, T5, T6) - 3, 4, 23:32
selected_electrodes = [2, 15:22];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirIn = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/AfterICA';
BehDirIn = '/cubric/collab/ccbrain/data/Raw_Data and subjects/Behavioural_data/';
DirOutEvents = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EventsRespLocked';
DirOutEpochs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EpochsRespLocked';

subjects = [1000,1001,1003:1014,1016:1022]; 

% 1:Equal100, 2:Equal80, 3:Equal20, 4:Equal80_20, 5:NotEqual100vs80, 
% 6:NotEqual100vs20, 7:NotEqual100vs80&20, 8:NotEqual80vs20, 9:Control
ConditionName = {'Equal100', 'Equal80', 'Equal20', 'Equal80_20', ...
    'NotEqual100vs80', 'NotEqual100vs20', 'NotEqual100vs80_20', ...
    'NotEqual80vs20', 'Control100', 'Control80', 'Control20'};

% To extract responses
Trig.E1 = {'212' '221'};
Trig.E2 = {'234' '243'};
Trig.E3 = {'256' '265'};
Trig.E4 = {'234' '243' '256' '265' };
Trig.E5 = {'213'  '214'  '223'  '224'  '231'  '232'  '241'  '242'};
Trig.E6 = {'215'  '216'  '225'  '226'  '251'  '252'  '261'  '262'};
Trig.E7 = {'213'  '214'  '215'  '216'  '223'  '224'  '225'  '226'  '231'  '232'  '241'  '242'  '251'  '252'  '261'  '262'};
Trig.E8 = {'235'  '246'  '253'  '264' };
Trig.E9 = {'217' '227'  '271' '272' };
Trig.E10 = {'237' '247' '273' '274' };
Trig.E11 = {'257' '267' '275' '276' };


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for sub_idx = 1:length(subjects)
    disp(sprintf('Analysis of subject: %i -------------------', subjects(sub_idx)))
    EEG = pop_loadset('filename',[num2str(subjects(sub_idx)) '.set'],'filepath', DirIn);
    EEG = eeg_checkset( EEG );

    % channels VEOG HEOG no needed anymore
    EEG = pop_select( EEG,'nochannel',{'VEOG' 'HEOG'});
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG,1,'overwrite', 'on','gui','off');
    
    EEG = pop_eegfiltnew(EEG, 35,1,826,0,[],0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG,1,'overwrite','on','gui','off'); 
    EEG = eeg_checkset( EEG );
    eeglab redraw
    % (?) cleaning artifacts here 
    if cleaning_flag
        channels_names = {};
        for ch = 1:length(EEG.chanlocs)
            channels_names = [channels_names EEG.chanlocs(ch).labels];
        end

        [ ALLEEG EEG CURRENTSET ] = artif_remove_interpolate(ALLEEG, EEG, channels_names);
    end
    % Load behavioural data to prepare a matrix with the event info
    load([BehDirIn, num2str(subjects(sub_idx)), '/result/', num2str(subjects(sub_idx)), '_EEG_Reward_1.mat'])
    A = ResponseArray;
    load([BehDirIn, num2str(subjects(sub_idx)), '/result/', num2str(subjects(sub_idx)), '_EEG_Reward_2.mat'])
    B = ResponseArray;
    load([BehDirIn, num2str(subjects(sub_idx)), '/result/', num2str(subjects(sub_idx)),  '_EEG_Reward_3.mat'])
    C = ResponseArray;
    load([BehDirIn, num2str(subjects(sub_idx)), '/result/', num2str(subjects(sub_idx)), '_EEG_Reward_4.mat'])
    D = ResponseArray;
    BehavData = [A;B;C;D];


    % Type info of each stimulus

    index_boundary = strncmp('boundary', { EEG.event.type }, 10);
    index_boundary = find(index_boundary==1); 
    for b = 1:length(index_boundary)
        EEG.event(1,index_boundary(b)).type = 0;
    end
    eeglab redraw

    z = 1;
    for iz = 1:length(EEG.event)
        if EEG.event(1,iz).type == '2'
            EEG.event(1,iz).type =  str2double([num2str(2), num2str(BehavData(z,16)), ...
                                                num2str(BehavData(z,17))]);
            EEG.event(1,iz).reactiontime = BehavData(z,9);
            z = z + 1;
        elseif EEG.event(1,iz).type == '1'
            EEG.event(1,iz).type = 1;
            EEG.event(1,iz).reactiontime = 0;
        end
    end
    clear iz z index_boundary
    
    % Clean events with empty reaction time
    vec2rm = [];
    for iz = 1:length(EEG.event)
        if isnan(EEG.event(1,iz).reactiontime)
            vec2rm = [vec2rm, iz];
        end
    end
    if length(vec2rm) > 0
       EEG.event(vec2rm) = []; 
    end
    clear iz vec2rm 
    % Dataset for each condition
    EEG = pop_saveset( EEG, 'filename', ...
        [num2str(subjects(sub_idx)), '.set'],'filepath', [DirOutEvents, '/']);
    for j= 1:length(ConditionName)
        field = ConditionName{j};
        EEG = pop_loadset('filename',[num2str(subjects(sub_idx)),'.set'],'filepath',[DirOutEvents, '/']);
        EEG = eeg_checkset( EEG );
        epoch_length = [-2.5 1];
        EEG = pop_epoch( EEG,   Trig.(['E' num2str(j)]) , epoch_length, 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );

        % Baseline calculation and removal
        Baseline = [];
        for k= 1:length(EEG.epoch)
            % Calculating baseline with eventlatency of the stimulus
            % (converting it first from ms in timeframes and then
            % substracting the time period of the fixation to get the
            % values for a baseline correction of -600 to -500)
            % 250Hz -> 1 Timeframes + 4 ms
            % -------------------------------------------------------------
            rt = EEG.epoch(1,k).eventreactiontime{length(EEG.epoch(1,k).eventreactiontime)};
            %Baseline = [Baseline; ((3000+q)/4) - 150, ((3000+q)/4) - 125];
            Baseline = [Baseline; round((-epoch_length(1)-rt-0.1)*EEG.srate), ...
                                  round((-epoch_length(1)-rt)*EEG.srate)];
        end

        X=[];

        % substraction of the baseline mean (X) from all of the channels
        % (DataEEG)
        for k= 1:length(EEG.epoch)
            x = Baseline(k,1);
            y = Baseline(k,2);
            X = [X,mean(EEG.data(:, x:y, k), 2)];
            for ii= 1:length(EEG.data)
                DataEEG = EEG.data(:,ii,k);
                DataEEG = DataEEG - X(:,k);
                EEG.data(:, ii, k) = DataEEG;
            end 
        end

        EEG = pop_saveset( EEG, 'filename',[num2str(subjects(sub_idx)), '_', field, '.set'],'filepath', [DirOutEpochs,'/']);
        EEG = eeg_checkset( EEG );
        FullData.(field).(['u' num2str(subjects(sub_idx))])(:,:,:) = EEG.data;
        Mean.(field) = mean(EEG.data,3);
        GlobalMean.(field)(sub_idx,:) = mean(Mean.(field),1);
        SelectedMean.(field)(sub_idx,:) = mean(Mean.(field)(selected_electrodes, :),1);
        WithElectrodesMean.(field)(sub_idx,:,:) = Mean.(field);
        GlobalVar.(field)(sub_idx,:) = var(Mean.(field), 1);
        GlobalStd.(field)(sub_idx,:) = std(Mean.(field), 1);
    end

end

clear ii iz j k ch

save('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedDataResp.mat', ...
    'subjects', 'epoch_length','GlobalMean', 'FullData', 'WithElectrodesMean');
%% Figures

% use that only when you load data from file:
% for sub_idx = 1:length(subjects)
%     for j= 1:length(ConditionName)
%         field = ConditionName{j};
%         Mean.(field) = mean(FullData.(field).(['u' num2str(subjects(sub_idx))])(:,:,:),3);
%         GlobalMean.(field)(sub_idx,:) = mean(Mean.(field),1);
%         SelectedMean.(field)(sub_idx,:) = mean(Mean.(field)(selected_electrodes, :),1);
%         GlobalStd.(field)(sub_idx,:) = std(Mean.(field), 1);
%     end
% end

WhatToPlot = GlobalMean;
tttt = -epoch_length(1) ;
% Plot mean of the three equal conditions in one graph 
tmpA = [mean(WhatToPlot.Equal100(:,:),1);
    mean(WhatToPlot.Equal80(:,:),1);
    mean(WhatToPlot.Equal20(:,:),1);
    mean(WhatToPlot.Control100(:,:),1);
    mean(WhatToPlot.Control80(:,:),1);
    mean(WhatToPlot.Control20(:,:),1)];

limits = 1.1*[min(min(tmpA)) max(max(tmpA))];
clear tmpA

figure;
suptitle('equal conditions (resp-locked)');
subplot(3,1,1);
A = [mean(WhatToPlot.Equal100(:,:),1);
    mean(WhatToPlot.Control100(:,:),1)];
plot_erp(A, EEG.srate, tttt, {'Equal100', 'Control100'}, limits);

subplot(3,1,2);
A = [mean(WhatToPlot.Equal80(:,:),1);
    mean(WhatToPlot.Control80(:,:),1)];
plot_erp(A, EEG.srate, tttt, {'Equal80', 'Control80'}, limits);

subplot(3,1,3);
A = [mean(WhatToPlot.Equal20(:,:),1);
    mean(WhatToPlot.Control20(:,:),1)];
plot_erp(A, EEG.srate, tttt, {'Equal20', 'Control20'}, limits);

% Plot mean of the three not-equal conditions in one graph 
tmpB = [mean(WhatToPlot.NotEqual100vs80(:,:),1);
    mean(WhatToPlot.NotEqual100vs20(:,:),1);
    mean(WhatToPlot.NotEqual80vs20(:,:),1);
    mean(WhatToPlot.Control100(:,:),1);
    mean(WhatToPlot.Control80(:,:),1);
    mean(WhatToPlot.Control20(:,:),1)];
limits = 1.1*[min(min(tmpB)) max(max(tmpB))];
clear tmpB

figure;
suptitle('not-equal conditions (resp-locked)');
subplot(3,1,1);
A = mean(WhatToPlot.NotEqual100vs80(:,:),1);
plot_erp(A, EEG.srate, tttt, {'NotEqual100vs80'}, limits);

subplot(3,1,2);
A = mean(WhatToPlot.NotEqual100vs20(:,:),1);
plot_erp(A, EEG.srate, tttt, {'NotEqual100vs20'}, limits);

subplot(3,1,3);
A = mean(WhatToPlot.NotEqual80vs20(:,:),1);
plot_erp(A, EEG.srate, tttt, {'NotEqual80vs20'}, limits);


% Plot global field power (GFP) of the three equal conditions in one graph 
tmpC= [mean(sqrt(GlobalVar.Equal100(:,:)),1);
    mean(sqrt(GlobalVar.Equal80(:,:)),1);
    mean(sqrt(GlobalVar.Equal20(:,:)),1);
    mean(sqrt(GlobalVar.Control100(:,:)),1);
    mean(sqrt(GlobalVar.Control80(:,:)),1);
    mean(sqrt(GlobalVar.Control20(:,:)),1)];
limits = 1.1*[min(min(tmpC)) max(max(tmpC))];
clear tmpC

figure;
suptitle('GFP of the three equal conditions (resp-locked)')
subplot(3,1,1);
A = [mean(sqrt(GlobalVar.Equal100(:,:)),1);
    mean(sqrt(GlobalVar.Control100(:,:)),1)];
plot_erp(A, EEG.srate, tttt, {'Equal100', 'Control100'}, limits, 0);

subplot(3,1,2);
A = [mean(sqrt(GlobalVar.Equal80(:,:)),1);
    mean(sqrt(GlobalVar.Control80(:,:)),1)];
plot_erp(A, EEG.srate, tttt, {'Equal80', 'Control80'}, limits, 0);

subplot(3,1,3);
A = [mean(sqrt(GlobalVar.Equal20(:,:)),1);
    mean(sqrt(GlobalVar.Control20(:,:)),1)];
plot_erp(A, EEG.srate, tttt, {'Equal20', 'Control20'}, limits , 0);

%% Comparison Triplot
A = [mean(WhatToPlot.Equal100(:,:),1);
    mean(WhatToPlot.Equal80(:,:),1);
    mean(WhatToPlot.Equal20(:,:),1)];
C= [mean(GlobalStd.Equal100(:,:),1);
    mean(GlobalStd.Equal80(:,:),1);
    mean(GlobalStd.Equal20(:,:),1)];
C = C./sqrt(size(GlobalStd.Equal20,1));
limits = 1.8*[min(min(A)) max(max(A))];


figure
hold on
h1 = plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(1,:)), 'Color', 'b', 'LineWidth', 1.2);
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(1,:))+squeeze(C(1,:)), '--', 'Color', 'b');
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(1,:))-squeeze(C(1,:)), '--', 'Color', 'b');
h2 = plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(2,:)), 'Color', 'r', 'LineWidth', 1.2);
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(2,:))+squeeze(C(2,:)), '--', 'Color', 'r');
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(2,:))-squeeze(C(2,:)), '--', 'Color', 'r');
h3 = plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(3,:)), 'Color', 'g', 'LineWidth', 1.2);
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(3,:))+squeeze(C(3,:)), '--', 'Color', 'g');
plot((1:size(A,2))/EEG.srate - tttt, squeeze(A(3,:))-squeeze(C(3,:)), '--', 'Color', 'g');
set(gca,'Ydir','reverse')
xlabel('time')
ylim(limits)
legend([h1 h2 h3], {'Eq100', 'Eq80', 'Eq20'}, 'Location', 'sw')
hold off
xlabel('time')
