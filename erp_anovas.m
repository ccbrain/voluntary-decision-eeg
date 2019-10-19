%Change for loop if you change dataset
load('GlobalAveragedDataStim.mat')

levels = {20, 80, 100};
conds  = {'Equal', 'Control'};
nr_subjects = 21;

p_conds_time = [];
p_levels_time = [];
F_conds_time = [];
F_levels_time = [];

for time_point = 1:350
    cond_col = {};
    lev_col = [];
    eeg_col = [];
    subj_col = [];
    ii=0;
    datamat = zeros(nr_subjects, length(conds), length(levels));

    e_cc = 0;
    for cc = conds
        e_cc = e_cc + 1;
        e_ll = 0;
        for ll = levels
            e_ll = e_ll + 1;
            cond = cc{1};
            lev  = ll{1};
            lev_col = [lev_col; repmat(lev, [nr_subjects 1])];
            subj_col = [subj_col 1:nr_subjects'];
            eeg_col = [eeg_col; GlobalMean.([cond num2str(lev)])(:,time_point)];
            for i=1:nr_subjects
               ii=ii+1;
               cond_col{ii} = cond;
               datamat(i,e_cc, e_ll) =  GlobalMean.([cond num2str(lev)])(i,time_point);
            end
        end
    end

    antable = table(subj_col', eeg_col, cond_col', lev_col);

    tbl = simple_mixed_anova(datamat);

    p_conds = tbl{3,5}
    p_levels = tbl{5,5}
    p_conds_time = [p_conds_time p_conds];
    p_levels_time = [p_levels_time p_levels];
    F_conds_time = [F_conds_time tbl{3,4}];
    F_levels_time = [F_levels_time tbl{5,4}];
end

%vox_threshold = 0.01;
%cpconds = clusterPermute2(p_conds_time, F_conds_time, vox_threshold);
%cplevels = clusterPermute2(p_levels_time, F_levels_time, vox_threshold);
cpconds = fdr_bh(p_conds_time);
cplevels = fdr_bh(p_levels_time);

xtime = (1:350)/250-0.4;

xtime(logical(cpconds))

