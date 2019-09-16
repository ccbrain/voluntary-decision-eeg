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

cfg_LDA            =  [];
cfg_LDA.classifier = 'lda';
cfg_LDA.param      = struct('lambda','auto');

cfg_SVM = mv_get_classifier_param('svm');
cfg_SVM.classifier = 'svm';
cfg_SVM.param.C          = 0.1;%[0.01 0.1 0.5];

cfg_LR            = [];
cfg_LR.classifier = 'logreg';
cfg_LR.param      = struct('lambda','auto' );

cfg_libSVM = mv_get_classifier_param('libsvm');
cfg_libSVM.classifier = 'libsvm';
cfg_libSVM.C          = 0.1;%[0.01 0.1 0.5];
cfg_libSVM.kernel_type = 0; % 0 linear, 2 rbf
%%%%%%%%%%%%%%%%%%%%%%%%%%Pick classifier
cfg_clf = cfg_SVM;

%%%%%%%%%%%%%%%% Hyperparams
smoothing_window = 0.02; % sec
timestamps       = [0.2, 0.35, 0.5];
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

tmpA = permute(FullData.pref.(subject_code),[3 1 2]);
tmpB = permute(FullData.notpref.(subject_code),[3 1 2]);
tmpA = addEvery2ndRows(tmpA);
tmpB = addEvery2ndRows(tmpB);
edata.comparePref = double([tmpA; tmpB]);
edata.comparePref = applyToDimension(@(xx) gaussfilt(timex, xx, smoothing_window), edata.comparePref, 3);
labels.comparePref = ones(1, size(edata.comparePref,1));
labels.comparePref(size(tmpA,1):size(edata.comparePref,1)) = 2;

avg_w = zeros(length(timestamps), size(edata.compare100,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.compare100(:,:,tt_idx));
    cc = train_svm(cfg_clf, dd,labels.compare100);
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
    cc = train_svm(cfg_clf, dd,labels.compare80);
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
    cc = train_svm(cfg_clf, dd,labels.compare20);
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
    cc = train_svm(cfg_clf, dd,labels.compare80100);
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
    cc = train_svm(cfg_clf, dd,labels.compare20100);
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
    cc = train_svm(cfg_clf, dd,labels.compare2080);
    % From S. Haufe 2014:
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.compare2080(sub_idx, :) = avg_w;

% comparePREF
avg_w = zeros(length(timestamps), size(edata.comparePref,2));
for tt = 1:length(timestamps)
    tx = find(timex >= timestamps(tt));
    tt_idx = tx(1);
    dd = squeeze(edata.comparePref(:,:,tt_idx));
    cc = train_svm(cfg_clf, dd,labels.comparePref);
    % From S. Haufe 2014:
    activation = transformToInterpretable(dd,cc);
    avg_w(tt,:) = activation;%cc.w;
end
avg_w = mean(avg_w,1);
all_weights.comparePref(sub_idx, :) = avg_w;

end
%save(['Data/ml_feature_importance_lsvm'], 'all_weights', 'timex', 'cfg_clf')

%%
% Plotting the resuts
load('Data/chanlocs');
fields = fieldnames(all_weights);
val.max = [];
val.min = [];
for i = 1:length(fields)
    val.max = [val.max max(abs(mean(all_weights.(fields{i}))))];
    val.min = [val.min min(abs(mean(all_weights.(fields{i}))))];
end

topoplot(abs(mean(all_weights.comparePref,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])

figure;
subplot(131)
topoplot(abs(mean(all_weights.compare100,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 100 vs Control 100')
subplot(132)
topoplot(abs(mean(all_weights.compare80,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 80 vs Control 80')
subplot(133)
topoplot(abs(mean(all_weights.compare80,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Control 20')
figure;
subplot(131)
topoplot(abs(mean(all_weights.compare80100,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 80 vs Equal 100')
subplot(132)
topoplot(abs(mean(all_weights.compare20100,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Equal 100')
subplot(133)
topoplot(abs(mean(all_weights.compare2080,1)), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Equal 80')

%%
fields = fieldnames(all_weights);
val.max = [];
val.min = [];
for i = 1:length(fields)
    val.max = [val.max max(mean(all_weights.(fields{i})))];
    val.min = [val.min min(mean(all_weights.(fields{i})))];
end

figure;
subplot(131)
topoplot(mean(all_weights.compare100,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 100 vs Control 100')
subplot(132)
topoplot(mean(all_weights.compare80,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 80 vs Control 80')
subplot(133)
topoplot(mean(all_weights.compare80,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Control 20')
figure;
subplot(131)
topoplot(mean(all_weights.compare80100,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 80 vs Equal 100')
subplot(132)
topoplot(mean(all_weights.compare20100,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Equal 100')
subplot(133)
topoplot(mean(all_weights.compare2080,1), chanlocs32, 'maplimits', [min(val.min), max(val.max)])
title('Equal 20 vs Equal 80')
