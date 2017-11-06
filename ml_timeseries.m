% LDA classification of ERPs in timewise manner
clear all
%%%%%%%%%%%%%
% Stimuli Locked or Response Locked
StimLocked = 1;
if StimLocked
    load('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedData.mat')
    timeshift = -0.4;
else
    load('/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/GlobalAveragedDataRespLocked.mat')
    timeshift = -2.5;
end
%%%%%%%%%%%%%
clear GlobalMean GlobalVar GlobalStd

edata.compare100 = [WithElectrodesMean.Equal100; WithElectrodesMean.Control100];
edata.compare80 = [WithElectrodesMean.Equal80; WithElectrodesMean.Control80];
edata.compare20 = [WithElectrodesMean.Equal20; WithElectrodesMean.Control20];
edata.compare80100 = [WithElectrodesMean.Equal80; WithElectrodesMean.Equal100];
edata.compare20100 = [WithElectrodesMean.Equal100; WithElectrodesMean.Equal20];
edata.compare2080 = [WithElectrodesMean.Equal100; WithElectrodesMean.Equal80];

labels= ones(1, size(edata.compare100,1));
labels(size(edata.compare100,1)/2:size(edata.compare100,1)) = 2;
srate = 250;
timex = (1:size(edata.compare100,3))/srate + timeshift;

%% Setup configuration struct for LDA and Logistic Regression

cfg_LDA =  [];
cfg_LDA.CV         = 'kfold';
cfg_LDA.K          = 5;
cfg_LDA.repeat     = 10;
cfg_LDA.classifier = 'lda';
cfg_LDA.param      = struct('lambda','auto');
cfg_LDA.metric     = 'acc';
cfg_LDA.balance = 'undersample';

cfg_LR =  cfg_LDA;
cfg_LR.classifier = 'logreg';
cfg_LR.param      = struct('lambda','auto' );


param = mv_get_classifier_param('svm');
param.C = 1;

cfg_SVM = cfg_LDA;
cfg_SVM.classifier = 'svm';
cfg_SVM.param = param;

%% Run classification Eq100 vs Ctr100
[auc_LR, result_LR] = mv_classify_across_time(cfg_SVM, edata.compare100, labels);

[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare100, labels);

mv_plot_result(result_LDA, timex)
title('Eq100 vs Contr100')

%% Run classification Eq80 vs Ctr80

[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare80, labels);

mv_plot_result(result_LDA, timex)
title('Eq80 vs Ctr80')

%% Run classification Eq20 vs Contr20

[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare20, labels);

mv_plot_result(result_LDA, timex)
title('Eq20 vs Contr20')

%% Run classification Eq100 vs Eq80
[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare80100, labels);

mv_plot_result(result_LDA, timex)

title('Eq100 vs Eq80')

%% Run classification Eq100 vs Eq20
cfg_LDA.CV = 'none'
[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare20100, labels);

mv_plot_result(result_LDA, timex)
title('Eq100 vs Eq20')

%% Run classification Eq20 vs Eq80
[auc_LDA, result_LDA] = mv_classify_across_time(cfg_LDA, edata.compare2080, labels);

mv_plot_result(result_LDA, timex)
title('Eq20 vs Eq80')
