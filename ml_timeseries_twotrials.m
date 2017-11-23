% LDA classification of ERPs in timewise manner among trials
clear all
%%%%%%%%%%%%%
% Stimuli Locked (1) or Response Locked (0)
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

%% Setup configuration struct for LDA and Logistic Regression

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
cfg_clf = cfg_libSVM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg_clf.CV         = 'kfold';
cfg_clf.K          = 5;
cfg_clf.repeat     = 8;
cfg_clf.metric     = 'acc';
cfg_clf.balance    = 'undersample';

%%%%%%%%%%%%%%%% Hyperparams
smoothing_window = 0.05; % sec
n_comp           = 6;
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

% Xdat = reduceInAllTimepoints(edata.compare100, n_comp);
% [auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare100);
% all_auc.compare100(sub_idx, :) = auc_clf;
% Xdat = reduceInAllTimepoints(edata.compare80, n_comp);
% [auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare80);
% all_auc.compare80(sub_idx, :) = auc_clf;
% Xdat = reduceInAllTimepoints(edata.compare20, n_comp);
% [auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare20);
% all_auc.compare20(sub_idx, :) = auc_clf;
Xdat = reduceInAllTimepoints(edata.compare80100, n_comp);
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare80100);
all_auc.compare80100(sub_idx, :) = auc_clf;
Xdat = reduceInAllTimepoints(edata.compare20100, n_comp);
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare20100);
all_auc.compare20100(sub_idx, :) = auc_clf;
Xdat = reduceInAllTimepoints(edata.compare2080, n_comp);
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, Xdat, labels.compare2080);
all_auc.compare2080(sub_idx, :) = auc_clf;
end
%save(['Data/ml_timeseries_2nd_gauss_results_' cfg_clf.classifier], 'all_auc', 'timex', 'cfg_clf')

%%
% Plotting the resuts
% resp lcoked: ['Data/ml_timeseries_resp_2nd_gauss_results_' cfg_clf.classifier]
%load(['Data/ml_timeseries_2nd_gauss_results_' cfg_clf.classifier])
y_limits = [0.4 0.7];
figure;
subplot(131);
accPlot(timex, all_auc.compare100)
title('Equal 100 vs Control 100')
ylim(y_limits);
subplot(132);
accPlot(timex, all_auc.compare80)
title('Equal 80 vs Control 80')
ylim(y_limits);
subplot(133);
accPlot(timex, all_auc.compare20)
title('Equal 20 vs Control 20')
ylim(y_limits);
figure;
subplot(131);
accPlot(timex, all_auc.compare80100)
title('Equal 80 vs Equal 100')
ylim(y_limits);
subplot(132);
accPlot(timex, all_auc.compare20100)
title('Equal 20 vs Equal 100')
ylim(y_limits);
subplot(133);
accPlot(timex, all_auc.compare2080)
title('Equal 20 vs Equal 80')
ylim(y_limits);
