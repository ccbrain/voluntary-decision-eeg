function [ ALLEEG EEG CURRENTSET ] = artif_remove_interpolate(ALLEEG,  EEG , standardchannels)
%ARTIF_REMOVE_INTERPOLATE automatically remove artifacts using clean_artifacts
% function and interpolate removed channels.
% ALLEEG - ALLEEG from eeglab
% EEG - EEGlab structure
% standardchannels - cell with channel names  

EEG = clean_artifacts(EEG, 'WindowCriterion', 'off','highpass_band','off', 'LineNoiseCriterion', 'off', 'ChannelCriterionMaxBadTime', 0.6); 

% Spherical Interpolation of removed channels

canales='';
for j=1:EEG.nbchan
    canales=[canales,EEG.chanlocs(1,j).labels,','];
end
channels = regexp(canales, ',', 'split');
canales  = find(ismember(standardchannels,channels)==1);
chan     = find(ismember(standardchannels,setdiff(standardchannels, channels))==1);
InterpolChan = standardchannels(chan);
EEG = pop_interp(EEG, ALLEEG(1).chanlocs([chan]), 'spherical');
EEG = eeg_checkset( EEG );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'overwrite','on','gui','off'); 

end

