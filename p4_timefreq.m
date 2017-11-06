% Script for epoch cutting and Stimuli Locked Analysis
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirIn = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/AfterICA';
BehDirIn = '/cubric/collab/ccbrain/data/Raw_Data and subjects/Behavioural_data/';
DirOutEvents = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EventsRespLocked';
DirOutEpochs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/EpochsRespLocked';
DirSpectFigs = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/figs/spectrograms/RespLocked';

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
BASELINE_LENGTH = 0.3; %s

for sub_idx = 1:length(subjects)
    disp(sprintf('Analysis of subject: %i -------------------', subjects(sub_idx)))

    for j= 1:length(ConditionName)
        field = ConditionName{j};
        if strmatch('Not',field)
           continue; 
        end
                field = ConditionName{j};
        EEG = pop_loadset('filename',[num2str(subjects(sub_idx)),'.set'],'filepath',[DirOutEvents, '/']);
        EEG = eeg_checkset( EEG );
        epoch_length = [-1 1];
        EEG = pop_epoch( EEG,   Trig.(['E' num2str(j)]) , epoch_length, 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        
        for nr_electrode = 1:length(EEG.chanlocs)
%             [ersp itc powbase tm fq] = pop_newtimef( EEG, 1, nr_electrode, [-3000   996], [0] , ...
%                 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
%                 'caption', [field ' electr:' EEG.chanlocs(nr_electrode).labels], ...
%                 'baseline',[-600 -500], 'freqs', [0 75], 'plotitc' , 'off', ...
%                 'plotphase', 'off', 'padratio', 1); %'plotersp', 'off',
            for epoch_ix = 1:length(EEG.epoch)
                [erspepoch fq tm] = specgram(EEG.data(nr_electrode,:,epoch_ix),round(EEG.srate/8),EEG.srate,[],20);
                rt = EEG.epoch(1,epoch_ix).eventreactiontime{length(EEG.epoch(1,epoch_ix).eventreactiontime)};
                if -rt < -1
                    baseline_range = [-0.9 -0.6];
                elseif -rt-BASELINE_LENGTH < -1
                    baseline_range = [-0.9 -rt+0.2];
                else
                    baseline_range = [-rt-BASELINE_LENGTH -rt];
                end
                erspepoch = abs(erspepoch);
                tm = tm + epoch_length(1);
                erspepoch = remove_spect_baseline(tm, erspepoch, baseline_range, 'div');
                if epoch_ix == 1
                    ersp = erspepoch;
                else
                    ersp = ersp + erspepoch;
                end
            end
            ersp = ersp / length(EEG.epoch);
            ersp = 10*log10(ersp);
            fig = plot_spectrogram(tm, fq, ersp, [field ' electr:' EEG.chanlocs(nr_electrode).labels],1);
            print(fig, [DirSpectFigs '/' field '/' [field '_' num2str(subjects(sub_idx)) '_' EEG.chanlocs(nr_electrode).labels]], '-dpng');
            close(fig);
            if nr_electrode == 1
                erspGlobal = ersp;
            else
                erspGlobal = erspGlobal + ersp;
            end
            spectrumFull.(field)(sub_idx, nr_electrode, :, :) = ersp;
            end
        erspGlobal = erspGlobal / length(EEG.chanlocs);
        spectrumOverElectrodes.(field)(sub_idx, :, :) = erspGlobal;
    end

    save('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrumRespLocked', ...
        'spectrumFull', 'spectrumOverElectrodes', 'tm', 'fq')
end

clear nr_electrode sub_idx j field ersp
