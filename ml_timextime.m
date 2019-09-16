% Time generalisation: classification of ERPs in timewise manner
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
cfg_clf = cfg_SVM;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg_clf.CV         = 'kfold';
cfg_clf.K          = 5;
cfg_clf.repeat     = 3;
cfg_clf.metric     = 'acc';
cfg_clf.balance    = 'undersample';
cfg_clf.normalise  = 'demean';

%% Loop with model training among the subjects
for sub_idx = 1 : length(subjects) %sub_idx = 1;
subject_code = ['u' num2str(subjects(sub_idx))];
disp(['----------------------- ' subject_code])

tmpA = permute(FullData.Control100.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare100 = double([tmpA; tmpB]);
labels.compare100= ones(1, size(edata.compare100,1));
labels.compare100(size(tmpA,1):size(edata.compare100,1)) = 2;

tmpA = permute(FullData.Control80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
edata.compare80 = double([tmpA; tmpB]);
labels.compare80= ones(1, size(edata.compare80,1));
labels.compare80(size(tmpA,1):size(edata.compare80,1)) = 2;

tmpA = permute(FullData.Control20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal20.(subject_code),[3 1 2]);
edata.compare20 = double([tmpA; tmpB]);
labels.compare20= ones(1, size(edata.compare20,1));
labels.compare20(size(tmpA,1):size(edata.compare20,1)) = 2;

tmpA = permute(FullData.Equal80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare80100 = double([tmpA; tmpB]);
labels.compare80100 = ones(1, size(edata.compare80100,1));
labels.compare80100(size(tmpA,1):size(edata.compare80100,1)) = 2;

tmpA = permute(FullData.Equal20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare20100 = double([tmpA; tmpB]);
labels.compare20100 = ones(1, size(edata.compare20100,1));
labels.compare20100(size(tmpA,1):size(edata.compare20100,1)) = 2;

tmpA = permute(FullData.Equal20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
edata.compare2080 = double([tmpA; tmpB]);
labels.compare2080 = ones(1, size(edata.compare2080,1));
labels.compare2080(size(tmpA,1):size(edata.compare2080,1)) = 2;

[acc_clf, result_acc] = mv_classify_timextime(cfg_clf, edata.compare100, labels.compare100);
all_auc.compare100(sub_idx, :, :) = acc_clf;
[auc_clf, result_clf] = mv_classify_timextime(cfg_clf, edata.compare80, labels.compare80);
all_auc.compare80(sub_idx, :, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_timextime(cfg_clf, edata.compare20, labels.compare20);
all_auc.compare20(sub_idx, :, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_timextime(cfg_clf, edata.compare80100, labels.compare80100);
all_auc.compare80100(sub_idx, :, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_timextime(cfg_clf, edata.compare20100, labels.compare20100);
all_auc.compare20100(sub_idx, :, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_timextime(cfg_clf, edata.compare2080, labels.compare2080);
all_auc.compare2080(sub_idx, :, :) = auc_clf;
end
save(['Data/ml_timextime_results_' cfg_clf.classifier], 'all_auc', 'timex', 'cfg_clf')

%%
load(['Data/ml_timextime_results_' cfg_clf.classifier], 'all_auc', 'timex', 'cfg_clf')
figure;
subplot(131)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare100(:, :, :),1)));
colormap jet
title('Accuracy: Equal 100 vs Control 100')
subplot(132)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare80(:, :, :),1)));
colormap jet
title('Accuracy: Equal 80 vs Control 80')
subplot(133)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare20(:, :, :),1)));
colormap jet
title('Accuracy: Equal 20 vs Control 20')
figure
subplot(131)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare80100(:, :, :),1)));
colormap jet
title('Accuracy: Equal 80 vs Equal 100')
subplot(132)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare20100(:, :, :),1)));
colormap jet
title('Accuracy: Equal 20 vs Equal 100')
subplot(133)
cfg_plot= [];
cfg_plot.x   = timex;
cfg_plot.y   = cfg_plot.x;
mv_plot_2D(cfg_plot, squeeze(mean(all_auc.compare2080(:, :, :),1)));
colormap jet
title('Accuracy: Equal 20 vs Equal 80')
