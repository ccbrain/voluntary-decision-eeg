% Script for epoch cutting and Stimuli Locked Analysis
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirIn = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/AfterICA';
BehDirIn = '/cubric/collab/ccbrain/data/Raw_Data and subjects/Behavioural_data/';
DirOutEvents = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EventsStimLocked';
DirOutEpochs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EpochsStimLocked';
DirSpectFigs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/figs/spectrograms/StimLocked';

subjects = [1000,1001,1003:1014,1016:1022]; 

% 1:Equal100, 2:Equal80, 3:Equal20, 4:Equal80_20, 5:NotEqual100vs80, 
% 6:NotEqual100vs20, 7:NotEqual100vs80&20, 8:NotEqual80vs20, 9:Control
ConditionName = {'Equal100', 'Equal80', 'Equal20', 'Equal80_20', ...
    'NotEqual100vs80', 'NotEqual100vs20', 'NotEqual100vs80_20', ...
    'NotEqual80vs20', 'Control100', 'Control80', 'Control20'};

% For extracting Stimuli
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
    
    eeglab redraw

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
            z = z + 1;
        elseif EEG.event(1,iz).type == '1'
            EEG.event(1,iz).type = 1;
        end
    end
    clear z index_boundary

    % Dataset for each condition
    EEG = pop_saveset( EEG, 'filename', ...
        [num2str(subjects(sub_idx)), '.set'],'filepath', [DirOutEvents, '/']);
    for j= 1:length(ConditionName)
        field = ConditionName{j};
        if strmatch('Not',field)
           continue; 
        end
                field = ConditionName{j};
        EEG = pop_loadset('filename',[num2str(subjects(sub_idx)),'.set'],'filepath',[DirOutEvents, '/']);
        EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG,   Trig.(['E' num2str(j)]) , [-3 1], 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        
        for nr_electrode = 1:length(EEG.chanlocs)
            fig = figure; 
            [ersp itc powbase tm fq] = pop_newtimef( EEG, 1, nr_electrode, [-600  1996], [0] , ...
                'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                'caption', [field ' electr:' EEG.chanlocs(nr_electrode).labels], ...
                'baseline', [0], 'freqs', [0 50], 'plotitc' , 'off', 'plotphase', 'off', 'padratio', 1);
            print(fig, [DirSpectFigs '/' field '/' [field '_' num2str(subjects(sub_idx)) '_' EEG.chanlocs(nr_electrode).labels]], '-dpng');
            close(fig);
            if nr_electrode == 1
                erspGlobal = ersp;
            else
                erspGlobal = erspGlobal + ersp;
            end
        end
        erspGlobal = erspGlobal / length(EEG.chanlocs);
        spectrumOverElectrodes.(field)(sub_idx, :, :) = erspGlobal;
    end
    save('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrumRespLocked', spectrumOverElectrodes, tm, fq)
end

clear ii iz j k ch
