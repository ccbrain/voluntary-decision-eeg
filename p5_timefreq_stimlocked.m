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

    for j= 1:length(ConditionName)
        field = ConditionName{j};
        if strmatch('Not',field)
           continue; 
        end
        EEG = pop_loadset('filename',[num2str(subjects(sub_idx)),'.set'],'filepath',[DirOutEvents, '/']);
        EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG,   Trig.(['E' num2str(j)]) , [-0.6 2], 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        for nr_electrode = 1:length(EEG.chanlocs)
            fig = figure; 
            [ersp{nr_electrode} itc powbase tm fq] = pop_newtimef( EEG, 1, nr_electrode, [-500  800], [0] , ...
                'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                'caption', [field ' electr:' EEG.chanlocs(nr_electrode).labels], ...
                'baseline', [-500 -300], 'freqs', [0 75], ...
                'plotitc' , 'off', 'plotphase', 'off', 'padratio', 1); %'plotersp', 'off',
            print(fig, [DirSpectFigs '/' field '/' [field '_' num2str(subjects(sub_idx)) '_' EEG.chanlocs(nr_electrode).labels]], '-dpng');
            close(fig);
            if nr_electrode == 1
                erspGlobal = ersp{nr_electrode};
            else
                erspGlobal = erspGlobal + ersp{nr_electrode};
            end
            spectrumFull.(field)(sub_idx, nr_electrode, :, :) = ersp{nr_electrode};
        end
        erspGlobal = erspGlobal / length(EEG.chanlocs);
        spectrumOverElectrodes.(field)(sub_idx, :, :) = erspGlobal;
    end
    save('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrumStimLocked', ...
        'spectrumFull', 'spectrumOverElectrodes', 'tm', 'fq')
end

clear nr_electrode sub_idx j field ersp
