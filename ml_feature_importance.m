% LDA classification of ERPs in timewise manner among trials
clear all
%%%%%%%%%%%%%
% Stimuli Locked or Response Locked
StimLocked = 1;
if StimLocked
    load('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedDataStim.mat')
    timeshift = epoch_length(1);
else
    load('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedDataResp.mat')
    timeshift = epoch_length(1);
end
%%%%%%%%%%%%%
clear GlobalMean GlobalVar GlobalStd

srate = 250;
timex = (1:size(FullData.Control100.u1000,2))/srate + timeshift;

%% Setup configuration struct for LDA

cfg_clf = mv_get_classifier_param('lda');


%%%%%%%%%%%%%%%% Hyperparams
smoothing_window = 0.02; % sec
timestamps       = [0.35];%[0.1 0.15 0.2 0.25 0.3 0.35];
%%%%%%%%%%%%%%%

%% Loop with model training among the subjects
for sub_idx = 1 : length(subjects) %sub_idx = 1;
subject_code = ['u' num2str(subjects(sub_idx))];
disp(['----------------------- ' subject_code])

tmpA = permute(FullData.Control100.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare100 = double([tmpA; tmpB]);
edata.compare100 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare100, 3);
% hold on; plot(squeeze(edata.compare100(2,4,:))); plot(squeeze(xw(2,4,:))); hold off;
labels.compare100= ones(1, size(edata.compare100,1));
labels.compare100(size(tmpA,1):size(edata.compare100,1)) = 2;

tmpA = permute(FullData.Control80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare80 = double([tmpA; tmpB]);
edata.compare80 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare80, 3);
labels.compare80= ones(1, size(edata.compare80,1));
labels.compare80(size(tmpA,1):size(edata.compare80,1)) = 2;

tmpA = permute(FullData.Control20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal20.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare20 = double([tmpA; tmpB]);
edata.compare20 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare20, 3);
labels.compare20= ones(1, size(edata.compare20,1));
labels.compare20(size(tmpA,1):size(edata.compare20,1)) = 2;

tmpA = permute(FullData.Equal80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare80100 = double([tmpA; tmpB]);
edata.compare80100 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare80100, 3);
labels.compare80100 = ones(1, size(edata.compare80100,1));
labels.compare80100(size(tmpA,1):size(edata.compare80100,1)) = 2;

tmpA = permute(FullData.Equal20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare20100 = double([tmpA; tmpB]);
edata.compare20100 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare20100, 3);
labels.compare20100 = ones(1, size(edata.compare20100,1));
labels.compare20100(size(tmpA,1):size(edata.compare20100,1)) = 2;

tmpA = permute(FullData.Equal20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.compare2080 = double([tmpA; tmpB]);
edata.compare2080 = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.compare2080, 3);
labels.compare2080 = ones(1, size(edata.compare2080,1));
labels.compare2080(size(tmpA,1):size(edata.compare2080,1)) = 2;

avg_w = zeros(length(timestamps), size(edata.compare100,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare100(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare100);
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare100(sub_idx, :) = avg_w;
%compare80
avg_w = zeros(length(timestamps), size(edata.compare80,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare80(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare80);
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare80(sub_idx, :) = avg_w;
%compare20
avg_w = zeros(length(timestamps), size(edata.compare20,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare20(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare20);
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare20(sub_idx, :) = avg_w;
% compare80100
avg_w = zeros(length(timestamps), size(edata.compare80100,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare80100(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare80100);
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);

all_weights.compare80100(sub_idx, :) = avg_w;

% compare20100
avg_w = zeros(length(timestamps), size(edata.compare20100,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare20100(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare20100);
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare20100(sub_idx, :) = avg_w;

% compare2080
avg_w = zeros(length(timestamps), size(edata.compare2080,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare2080(:,:,tt_idx));
    cc = train_lda(cfg_clf, dd,labels.compare2080);
    % From S. Haufe 2014:
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare2080(sub_idx, :) = avg_w;
end
save(['Data/ml_feature_importance_lda'], 'all_weights', 'timex', 'cfg_clf')

%%
% Plotting the resuts
load('Data/chanlocs');
figure;
topoplot(mean(all_weights.compare100,1), chanlocs32)
title('Equal 100 vs Control 100')
figure;
topoplot(mean(all_weights.compare80,1), chanlocs32)
title('Equal 80 vs Control 80')
figure;
topoplot(mean(all_weights.compare80,1), chanlocs32)
title('Equal 20 vs Control 20')
figure;
topoplot(mean(all_weights.compare80100,1), chanlocs32)
title('Equal 80 vs Equal 100')
figure;
topoplot(mean(all_weights.compare20100,1), chanlocs32)
title('Equal 20 vs Equal 100')
figure;
topoplot(mean(all_weights.compare2080,1), chanlocs32)
title('Equal 20 vs Equal 80')