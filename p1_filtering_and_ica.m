% Script for filtering the data and ICA channels removal done manually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all
DirIn = '/cubric/collab/ccbrain/data';                      % Path of where the data is stored (with raw data)
DirOut = '/cubric/collab/ccbrain/data/Prepro_250Hz';        % Path of where preprocessed data will be saved 
DirAfterIca = '/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/AfterICA025';
subject = [1000,1001,1003:1014,1016:1022]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick number of subject here:
%i = 15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

EEG = pop_loadset('filename',[num2str(subject(i)), '.set'],'filepath', ...
    [DirIn,'/Raw_Data and subjects/Trimmed EEG Data/']);
EEG = eeg_checkset( EEG );

% EOG channels manipulation
% As a result we get 2 horizontal and 2 vertical eyechannel to obtain 
% 2 eyechannels (VEOG:33, HEOG:34; Mastoids: 35,36)
tmp = EEG.data;
tmp(33,:) = tmp(33,:) - tmp(34,:);
tmp(34,:) = tmp(35,:) - tmp(36,:);
tmp(35,:) = tmp(37,:);
tmp(36,:) = tmp(38,:);
tmp = tmp(1:36,:);
EEG.data = tmp;
clear tmp

% Put template of 36 channel locations 
EEG.chanlocs= readlocs([DirIn,'/Chan_Template/chan36.elp']);
EEG = eeg_checkset( EEG );

%rereference to mastoids
EEG = pop_reref( EEG, [35 36] );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );

% filtering
% 0.5-100
EEG = pop_eegfiltnew(EEG, 0.25, 100, 38000, 0,[],0);
% notch filter 58-62
EEG = pop_eegfiltnew(EEG, 48, 52, 12000, 1,[],0);

EEG = pop_resample( EEG, 250);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 

% Run ICA on channels without taking into acount eyechannels

ica_channels = [1:32]; % Specifies the channel range for the ICA analysis
EEG = pop_runica(EEG, 'chanind', ica_channels);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
eeglab redraw


EOG = EEG.data(33:34,:);
% Compute temporal course of ICA components (without eyechannels: VEOG-33, HEOG-34)
DataICA = EEG.icaweights*EEG.data(1:32,:);

RHO_VEOG =[];
RHO_HEOG =[];
for k =1:32
    RHO_VEOG(k) = corr(DataICA(k,:)',EOG(1,:)');
end
for k =1:32
    RHO_HEOG(k) = corr(DataICA(k,:)',EOG(2,:)');
end

figure;
subplot(211)
bar(abs(RHO_VEOG))
ylim([0,1])
title('VEOG')
subplot(212)
bar(abs(RHO_HEOG))
ylim([0,1])
title('HEOG')

[value_HEOG,xHEOG] = sort(abs(RHO_HEOG),2, 'descend');
[value_VEOG,xVEOG] = sort(abs(RHO_VEOG),2, 'descend');

OEEG = pop_subcomp(EEG, unique([xHEOG(1:3) xVEOG(1:3)]));
OEEG = eeg_checkset( OEEG );
pop_saveset( OEEG, 'filename',[num2str(subject(i)), '.set'],...
    'filepath', [DirAfterIca, '/']);
close all
clear EEG OEEG xHEOG xVEOG RHO_HEOG RHO_VEOG