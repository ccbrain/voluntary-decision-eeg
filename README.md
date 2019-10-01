# Voluntary Decision project

This is code with the analysis for the project: *Breaking deadlocks: reward probability and spontaneous preference affected voluntary decision processes*.

## EEG analysis

You should read this repository in the following way:

- `p1_filtering_and_ica.m` - is the code used for cleaning and filtering the raw EEG data.
- `p2_events.m` - is the code for cutting the data in Response-Locked condition and making ERP plots.
- `p3_stimlocked.m` - is the code for cutting the data in Stimuli-Locked condition and making ERP plots.
- `p4_timefreq.m` and `p5_timefreq_stimlocked.m` contain our time-frequency analysis (not included in the final study).
- `ml_timeseries.m` and the rest of files with `ml_` prefix relate to the classification analysis between experimental condition in this study.
- `ml_feature_importance.m` script was used to determine feature importance of L-SVM classifiers.
- `single_trial_erp.m` - script for single trial SVD decomposition;

The raw data for this study you minght find online: TODO.

## Behavioural analysis

We used the above results to compute ERP correlated with bahavioural model of decision making. The R/STAN script for this you might find here:
```
Behavavioral_Modelling/
    Behavioral_modelling.R
    summarySE2.R
```

The folder above contains also data for the behavioural model fitting step.

## Data
Data for this project we uploaded here: .
