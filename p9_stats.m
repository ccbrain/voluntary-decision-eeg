% Frequency band analysis - in particular BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
how_locked = 'StimLocked';
brain_region = 'central';
load(['/cubric/collab/ccbrain/data/Scripts/eeg_analysis2/Data/spectrum' how_locked])
permuteP_alpha  = 0.05;
vox_threshold   = 0.05;
srate           = 250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
electrodes.frontal = [1, 5:14];
electrodes.central = [2, 15:22];
electrodes.posterior = [3, 4, 23:32];
electrodes.all = 1:32;
regions = {'frontal', 'central', 'posterior', 'all'};

cond1 = 'Equal80';
cond2 = 'Equal100';

% Average over electrodes from particular brain region

dat1 = squeeze(mean(spectrumFull.(cond1)(:,electrodes.(brain_region),:,:),2));
dat2 = squeeze(mean(spectrumFull.(cond2)(:,electrodes.(brain_region),:,:),2));

% Ttest
[~, p, ~, stats] = ttest(dat1, dat2);
p = squeeze(p);
tstat = squeeze(stats.tstat);
%plot_spectrogram(tm,fq, p);

%cc = clusterPermute2D(p, tstat, permuteP_alpha, vox_threshold);

%plot_spectrogram(tm,fq, cc.Mask);

%%
% Permutation test

% dat_tot = [dat1;dat2];
% permuteP_alpha=0.05;
% Nrep = 10000;
% Nrsubj = size(dat_tot,1);
% dat_resampled = zeros(Nrep, size(dat_tot,2),size(dat_tot,3));
% for rr = 1: Nrep
%     id1=randi(size(dat_tot,1),Nrsubj,1); id2=randi(size(dat_tot,1),Nrsubj,1);
%     dat_resampled(rr,:,:) = abs(mean(dat_tot(id1,:,:)) - mean(dat_tot(id2,:,:)));
% end
% 
% sidak_correction = 1 - (1-permuteP_alpha)^(1./(size(dat_tot,2)*size(dat_tot,3)));
% 
% dat_abs_diff = squeeze(abs(mean(dat1,1)-mean(dat2,1)));
% sigmat = zeros(size(dat_tot,2),size(dat_tot,3));
% for xf=1:size(dat_tot,2)
%     for xt=1:size(dat_tot,3)
%         prc = prctile(dat_resampled(:,xf,xt),100-round(sidak_correction*100));
%         if dat_abs_diff(xf,xt) > prc
%             sigmat(xf,xt) = 1;
%         end
%         
%         sigmatP(xf,xt) = (sum(dat_abs_diff(xf,xt)<squeeze(dat_resampled(:,xf,xt)))+1)/5001;
%     end
% end

figure;
subplot(221)
plot_spectrogram(tm,fq,squeeze(mean(dat1,1)), cond1, 2, [-1.5 1.5]);
subplot(222)
plot_spectrogram(tm,fq,squeeze(mean(dat2,1)), cond2, 2, [-1.5 1.5]);
subplot(223)
%plot_spectrogram(tm,fq,dat_abs_diff,'Abs difference',2);
p(find(p>=0.01)) = 1;
plot_spectrogram(tm,fq, p, 'Ttest p-val',2);
subplot(224)
%plot_spectrogram(tm,fq,sigmat,'Perm test stats',2);
plot_spectrogram(tm,fq, cc.Mask, 'Cluster correction',2);
