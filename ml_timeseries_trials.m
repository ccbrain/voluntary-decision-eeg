% LDA classification of ERPs in timewise manner
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

%cfg_clf =  [];
%cfg_clf.classifier = 'lda';
%cfg_clf.param      = struct('lambda','auto');

cfg_clf = mv_get_classifier_param('svm');

cfg_clf.CV         = 'kfold';
cfg_clf.K          = 5;
cfg_clf.repeat     = 8;
cfg_clf.metric     = 'acc';
cfg_clf.balance    = 'undersample';

cfg_LR =  cfg_clf;
cfg_LR.classifier = 'logreg';
cfg_LR.param      = struct('lambda','auto' );

%% Loop with model training among the subjects
for sub_idx = 1 : length(subjects) %sub_idx = 1;
subject_code = ['u' num2str(subjects(sub_idx))];
disp(['----------------------- ' subject_code])

tmpA = permute(FullData.Control100.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare100 = [tmpA; tmpB];
labels.compare100= ones(1, size(edata.compare100,1));
labels.compare100(size(tmpA,1):size(edata.compare100,1)) = 2;

tmpA = permute(FullData.Control80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
edata.compare80 = [tmpA; tmpB];
labels.compare80= ones(1, size(edata.compare80,1));
labels.compare80(size(tmpA,1):size(edata.compare80,1)) = 2;

tmpA = permute(FullData.Control20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal20.(subject_code),[3 1 2]);
edata.compare20 = [tmpA; tmpB];
labels.compare20= ones(1, size(edata.compare20,1));
labels.compare20(size(tmpA,1):size(edata.compare20,1)) = 2;

tmpA = permute(FullData.Control80.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare80100 = [tmpA; tmpB];
labels.compare80100 = ones(1, size(edata.compare80100,1));
labels.compare80100(size(tmpA,1):size(edata.compare80100,1)) = 2;

tmpA = permute(FullData.Control20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal100.(subject_code),[3 1 2]);
edata.compare20100 = [tmpA; tmpB];
labels.compare20100 = ones(1, size(edata.compare20100,1));
labels.compare20100(size(tmpA,1):size(edata.compare20100,1)) = 2;

tmpA = permute(FullData.Control20.(subject_code),[3 1 2]);
tmpB = permute(FullData.Equal80.(subject_code),[3 1 2]);
edata.compare2080 = [tmpA; tmpB];
labels.compare2080 = ones(1, size(edata.compare2080,1));
labels.compare2080(size(tmpA,1):size(edata.compare2080,1)) = 2;

[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare100, labels.compare100);
all_auc.compare100(sub_idx, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare80, labels.compare80);
all_auc.compare80(sub_idx, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare20, labels.compare20);
all_auc.compare20(sub_idx, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare80100, labels.compare80100);
all_auc.compare80100(sub_idx, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare20100, labels.compare20100);
all_auc.compare20100(sub_idx, :) = auc_clf;
[auc_clf, result_clf] = mv_classify_across_time(cfg_clf, edata.compare2080, labels.compare2080);
all_auc.compare2080(sub_idx, :) = auc_clf;
end
save('Data/ml_timeseries_results_svm_linear', 'all_auc', 'timex', 'cfg_clf')

%%
% Plotting the resuts
% load('Data/ml_timeseries_results')
accPlot(timex, all_auc.compare100)
title('Equal 100 vs Control 100')
accPlot(timex, all_auc.compare80)
title('Equal 80 vs Control 80')
accPlot(timex, all_auc.compare20)
title('Equal 20 vs Control 20')
accPlot(timex, all_auc.compare80100)
title('Equal 80 vs Equal 100')
accPlot(timex, all_auc.compare20100)
title('Equal 20 vs Equal 100')
accPlot(timex, all_auc.compare80100)
title('Equal 20 vs Equal 80')
