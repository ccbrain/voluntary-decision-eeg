% Script for epoch cutting and Stimuli Locked Analysis
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cleaning_flag = 0;
% Frontal (Fz, Fp1, Fp2, F1 –F8)           - 1, 5:14
% Central (Cz, C1-C6, T7, T8)              - 2, 15:22
% Posterior (Pz, P1-P6, Oz, O1-O2, T5, T6) - 3, 4, 23:32
selected_electrodes = [1 5:14];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirIn = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/AfterICA';
BehDirIn = '/cubric/collab/ccbrain/data/Raw_Data and subjects/Behavioural_data/';
DirOutEvents = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EventsStimLocked';
DirOutEpochs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EpochsStimLocked';

subjects = [1000,1001,1003:1014,1016:1022]; 

% 1:Equal100, 2:Equal80, 3:Equal20, 4:Equal80_20, 5:NotEqual100vs80, 
% 6:NotEqual100vs20, 7:NotEqual100vs80&20, 8:NotEqual80vs20, 9:Control
ConditionName = {'Equal100', 'Equal80', 'Equal20', 'Equal80_20', ...
    'NotEqual100vs80', 'NotEqual100vs20', 'NotEqual100vs80_20', ...
    'NotEqual80vs20', 'Control100', 'Control80', 'Control20'};

% For extracting Stimuli
Trig.E1 = {'12' '21'};
Trig.E2 = {'34' '43'};
Trig.E3 = {'56' '65'};
Trig.E4 = {'34' '43' '56' '65' };
Trig.E5 = {'13'  '14'  '23'  '24'  '31'  '32'  '41'  '42'};
Trig.E6 = {'15'  '16'  '25'  '26'  '51'  '52'  '61'  '62'};
Trig.E7 = {'13'  '14'  '15'  '16'  '23'  '24'  '25'  '26'  '31'  '32'  '41'  '42'  '51'  '52'  '61'  '62'};
Trig.E8 = {'35'  '46'  '53'  '64' };
Trig.E9 = {'17' '27'  '71' '72' };
Trig.E10 = {'37' '47' '73' '74' };
Trig.E11 = {'57' '67' '75' '76' };


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
        if  EEG.event(1,iz).type == '1'
            EEG.event(1,iz).type =  str2double([num2str(BehavData(z,16)), ...
                                                num2str(BehavData(z,17))]); 
            z = z + 1;
        elseif EEG.event(1,iz).type == '2'
            EEG.event(1,iz).type = 2;
        end
    end
    clear z index_boundary

    % Dataset for each condition
    EEG = pop_saveset( EEG, 'filename', ...
        [num2str(subjects(sub_idx)), '.set'],'filepath', [DirOutEvents, '/']);
    for j= 1:length(ConditionName)
        field = ConditionName{j};
        EEG = pop_loadset('filename',[num2str(subjects(sub_idx)),'.set'],'filepath',[DirOutEvents, '/']);
        EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG,   Trig.(['E' num2str(j)]) , [-0.6 2], 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        EEG = pop_rmbase( EEG, [-600    -500]);
        EEG = eeg_checkset( EEG );
        EEG = pop_saveset( EEG, 'filename',[num2str(subjects(sub_idx)), '_', field, '.set'],'filepath', [DirOutEpochs,'/']);
        Mean.(field) = mean(EEG.data,3);
        GlobalMean.(field)(sub_idx,:) = mean(Mean.(field),1);
        SelectedMean.(field)(sub_idx,:) = mean(Mean.(field)(selected_electrodes, :),1);
        GlobalVar.(field)(sub_idx,:) = var(Mean.(field));
    end

end

clear ii iz j k ch

save('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedData.mat', ...
    'subjects', 'GlobalMean','GlobalVar'); 

%% Figures

WhatToPlot = SelectedMean;

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
suptitle('equal conditions (stim-locked)');
subplot(3,1,1);
A = [smooth(mean(WhatToPlot.Equal100(:,:),1))';
    smooth(mean(WhatToPlot.Control100(:,:),1))'];
plot_erp(A, EEG.srate, 0.6, {'Equal100', 'Control100'}, limits);

subplot(3,1,2);
A = [smooth(mean(WhatToPlot.Equal80(:,:),1))';
    smooth(mean(WhatToPlot.Control80(:,:),1))'];
plot_erp(A, EEG.srate, 0.6, {'Equal80', 'Control80'}, limits);

subplot(3,1,3);
A = [smooth(mean(WhatToPlot.Equal20(:,:),1))';
    smooth(mean(WhatToPlot.Control20(:,:),1))'];
plot_erp(A, EEG.srate, 0.6, {'Equal20', 'Control20'}, limits);

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
suptitle('not-equal conditions (stim-locked)');
subplot(3,1,1);
A = smooth(mean(WhatToPlot.NotEqual100vs80(:,:),1))';
plot_erp(A, EEG.srate, 0.6, {'NotEqual100vs80'}, limits);

subplot(3,1,2);
A = smooth(mean(WhatToPlot.NotEqual100vs20(:,:),1))';
plot_erp(A, EEG.srate, 0.6, {'NotEqual100vs20'}, limits);

subplot(3,1,3);
A = smooth(mean(WhatToPlot.NotEqual80vs20(:,:),1))';
plot_erp(A, EEG.srate, 0.6, {'NotEqual80vs20'}, limits);


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
suptitle('GFP of the three equal conditions (stim-locked)')
subplot(3,1,1);
A = [mean(sqrt(GlobalVar.Equal100(:,:)),1);
    mean(sqrt(GlobalVar.Control100(:,:)),1)];
plot_erp(A, EEG.srate, 0.6, {'Equal100', 'Control100'}, limits, 0);

subplot(3,1,2);
A = [mean(sqrt(GlobalVar.Equal80(:,:)),1);
    mean(sqrt(GlobalVar.Control80(:,:)),1)];
plot_erp(A, EEG.srate, 0.6, {'Equal80', 'Control80'}, limits, 0);

subplot(3,1,3);
A = [mean(sqrt(GlobalVar.Equal20(:,:)),1);
    mean(sqrt(GlobalVar.Control20(:,:)),1)];
plot_erp(A, EEG.srate, 0.6, {'Equal20', 'Control20'}, limits , 0);
