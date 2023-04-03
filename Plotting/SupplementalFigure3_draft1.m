%Supplemental Figure 3
%first load data to get the correct indexing of motifs; 
%EntropyFigure
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
fp = fig_params_vpa; 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\SupplementalFig3';
if ~exist(savedir)
    mkdir(savedir);
end

% Load refitting data
file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
% file_list = GrabFiles(['\w*','dynamics3','\w*'],0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
grp = isVPA(unique(mouse));
data_all = cellfun(@(x) load(x,'stats_refit'),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.stats_refit,data_all,'UniformOutput',0);
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});

%zscore
data_mouse = zscore(data_mouse,[],1);

%Build neurotype axis
[nt, Y_null, betas, betas_null, allstats,bias] = DistanceToHyperplane(data_mouse,grp'+1,5,0,0.4);
[betas, idx_of_final_motifs] = sort(betas,'ascend');


%% now load all the individual motifs 
motif_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
gp = general_params_vpa;
file_list = GrabFiles('\w*chunk\w*.mat', 0, {motif_dir});
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);

%load all the motifs
temp = cellfun(@(x) load(x,'W','nanpxs','stats_train','stats_test'),file_list,'UniformOutput',0);
trainpev = cellfun(@(x) x.stats_train.pev,temp,'UniformOutput',1);
trainpev = arrayfun(@(x) nanmean(trainpev(mouse==x)), unique(mouse),'UniformOutput',1);
testpev = cellfun(@(x) x.stats_test.pev,temp,'UniformOutput',1);
testpev = arrayfun(@(x) nanmean(testpev(mouse==x)), unique(mouse),'UniformOutput',1);
lambdaval = cellfun(@(x) x.stats_train.lambda, temp,'UniformOutput',1);
numbermotifs = cellfun(@(x) x.stats_train.n_motifs, temp,'UniformOutput',1);
W = cellfun(@(x) x.W, temp,'UniformOutput',0);
motif_id = arrayfun(@(n) repmat(mouse(n),size(W{n},2),1), 1:numel(W),'UniformOutput',0);%get the mouse id for each motif
motif_id = cat(1,motif_id{:});
nanpxs = cellfun(@(x) x.nanpxs, temp,'UniformOutput',0);

%Recondition W
for i = 1:numel(W)    
   temp = zeros(gp.pixel_dim(1)*gp.pixel_dim(1),size(W{i},2),size(W{i},3));
   idx = ones(size(temp,1),1);
   idx(nanpxs{i})=0; %NOTE, NAN Will only use the same across all pixels, zeros will smooth essentailly for clustering 
   temp(idx==1,:,:) = W{i};
   W{i} = temp;
end

W = cat(2,W{:});
nanpxs = find(nanvar(reshape(W,[size(W,1),size(W,2)*size(W,3)]),[],2)<=eps);
W(nanpxs,:,:) = [];

%load the basis motifs
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL\basismotifs.mat');

W_alligned = AllignW(W,core_comm_idx,lags,cluster_idx,lag_mat);
W_alligned = W_alligned(:,:,nanvar(squeeze(sum(W_alligned,1)),[],1)>eps);

%% plot the maximum xcorr within and across groups
N = size(W_basis,2);

%ignoring the noiseclusters
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];
    
for i = 1:N
   if sum(ismember(noise_clusters,i))==0
       %idx of current cluster   
       temp_mat = tcorr_mat(cluster_idx==i,cluster_idx==i);
       grp_idx = isVPA(motif_id(cluster_idx==i));
       %remove diagnol
       temp_mat(1:1+size(temp_mat,1):end) = NaN;   
       animal_idx = motif_id(cluster_idx==i);
       
       %between groups
       sal_vpa = temp_mat(grp_idx==0,grp_idx==1);
       sal_vpa_avg = fisherInverse(nanmean(fisherZ(sal_vpa(:))));
       sal_vpa_sem = fisherInverse(sem(fisherZ(sal_vpa(:))));

       %between all animals (but across animals)
       all_within = temp_mat;   
       %get a mask of the same animal connections
       temp_animal_idx = animal_idx;   
       unique_temp = unique(temp_animal_idx); 
       between_mat = cell(1,numel(unique_temp));
       for j = 1:numel(unique_temp)
           temp_idx = temp_animal_idx==unique_temp(j);
           temp = all_within(temp_idx,temp_idx==0);
           between_mat{j} = temp(:);
       end
       all_within = cat(1,between_mat{:});   
       all_within_avg = fisherInverse(nanmean(fisherZ(all_within(:))));
       all_within_sem = fisherInverse(sem(fisherZ(sal_vpa(:))));
       
%        perm_stat = NaN(1,1000);
%        temp_grp = cat(1,ones(numel(sal_vpa(:)),1),2*ones(numel(all_within(:)),1));
%        for j = 1:1000
%           temp = cat(1,sal_vpa(:),all_within(:));        
%           temp = temp(randperm(numel(temp),numel(temp)));
%           perm_stat(j) = nanmean(temp(temp_grp==1))-nanmean(temp(temp_grp==2)); 
%        end
%        diff_stat = nanmean(sal_vpa(:))-nanmean(all_within(:));
%        if diff_stat>0
%            pval = sum([perm_stat,diff_stat]>=diff_stat)/(1000+1);
%        else
%            pval = sum([perm_stat,diff_stat]<=diff_stat)/(1000+1);
%        end
             
       pval = ranksum(sal_vpa(:),all_within(:));

       %plot
       figure; hold on;
       %sal       
       b = bar(1,all_within_avg,'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
       b.CData(1,:) = [0.3 0.3 0.3];
       errorbar(1,all_within_avg,all_within_sem,'LineWidth',1.25,'Color','k');
       %both
       b = bar(2,sal_vpa_avg,'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
       b.CData(1,:) = [0.6 0.6 0.6];
       errorbar(2,sal_vpa_avg,sal_vpa_sem,'LineWidth',1.25,'Color','k');
       
       title(sprintf('Motif %d pval=%0.2g',idx_of_final_motifs(motif_list(i)),pval)); 
       %format
       set(gca,'xlim',[0.25, 2.75],'units','centimeters','position',[4,3,1,2])
       set(gca,'ylim',[0.2 0.8]);
       fp.FormatAxes(gca)   
       
   end
end

%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','Individualmotifdiff',savedir,1); close all

%% Plot the average number of discovred motifs 
figure; hold on; 
pos = rand(20,1)/2-0.25;
b = bar(1,nanmean(trainpev),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.5 0.5 0.5];
% errorbar(1,nanmean(trainpev),sem(trainpev(:)),'LineWidth',1.25,'Color','k');
plot(pos+1,trainpev,'.','markersize',fp.m_markersize,'color',[0.5 0.5 0.5])
b = bar(2,nanmean(testpev),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.5 0.5 0.5];
% errorbar(2,nanmean(testpev),sem(testpev(:)),'LineWidth',1.25,'Color','k');
plot(pos+2,testpev,'.','markersize',fp.m_markersize,'color',[0.5 0.5 0.5])    

% %permuation test comparing them 
% perm_stat = NaN(1,1000);
% temp_grp = cat(1,ones(numel(trainpev(:)),1),2*ones(numel(testpev(:)),1));
% for j = 1:1000
%   temp = cat(1,trainpev(:),testpev(:));        
%   temp = temp(randperm(numel(temp),numel(temp)));
%   perm_stat(j) = nanmean(temp(temp_grp==1))-nanmean(temp(temp_grp==2)); 
% end
% diff_stat = nanmean(trainpev(:))-nanmean(testpev(:));
% if diff_stat>0
%    pval = sum([perm_stat,diff_stat]>=diff_stat)/(1000+1);
% else
%    pval = sum([perm_stat,diff_stat]<=diff_stat)/(1000+1);
% end

avg_train = nanmean(trainpev);
ci_train = bootci(1000,@nanmean,trainpev);
avg_test = nanmean(testpev);
ci_test = bootci(1000,@nanmean,testpev);


pval = ranksum(trainpev(:),testpev(:));

%format
title(sprintf('pval %0.2g tailed',pval));
set(gca,'xlim',[0.5, 2.5],'ylim',[0 1],'units','centimeters','position',[3.25 2 1.5 2])
fp.FormatAxes(gca)

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','traintestpev',savedir,1); close all

%% plot distribution of lambda values chosen
figure; hold on; 
histogram(lambdaval,'BinWidth',0.001,'EdgeAlpha',0,'FaceColor',[0.15 0.15 0.15],'FaceAlpha',0.3)
yval = get(gca,'ylim');
line([nanmean(lambdaval),nanmean(lambdaval)],yval,'color','r','linestyle','-','linewidth',2);
set(gca,'ylim',yval);
ylabel('# Shuffles')
xlabel('\lambda','Interpreter','tex')
set(gca,'units','centimeters','position',[3.25 2 2.25 2.25]);
fp.FormatAxes(gca)
set(gca,'xscale','linear')

%% Plot the number of motifs discovered per epoch
figure; hold on; 
histogram(numbermotifs,'BinWidth',1,'EdgeAlpha',0,'FaceColor',[0.15 0.15 0.15],'FaceAlpha',0.3)
yval = get(gca,'ylim');
line([nanmean(numbermotifs),nanmean(numbermotifs)],yval,'color','r','linestyle','-','linewidth',2);
set(gca,'ylim',yval);
ylabel('# Shuffles')
xlabel('# Motifs')
set(gca,'units','centimeters','position',[3.25 2 1.5 2]);
fp.FormatAxes(gca)
set(gca,'xscale','linear')

%% Plot the duration of each motif
temp = squeeze(nanmean(W_alligned,1));    
figure; hold on; 
shadedErrorBar(1:size(temp,2),nanmean(temp,1),sem(temp,1),'lineprops',{'color',[0.25 0.25 0.25]},'transparent',1,'patchSaturation',0.075);
plot(nanmean(temp),'linewidth',2,'color','k')
yval = get(gca,'ylim');
plot([8,8+13],[yval(2)-0.005,yval(2)-0.005],'linewidth',2,'color','r')
set(gca,'units','centimeters','position',[3.25 2 1.5 2],'xtick',[1,26]);
fp.FormatAxes(gca)


temp = squeeze(nanmean(W,1));    
figure; hold on; 
shadedErrorBar(1:size(temp,2),nanmean(temp,1),sem(temp,1),'lineprops',{'color',[0.25 0.25 0.25]},'transparent',1,'patchSaturation',0.075);
plot(nanmean(temp),'linewidth',2,'color','k')
yval = get(gca,'ylim');
plot([8,8+13],[yval(2)-0.005,yval(2)-0.005],'linewidth',2,'color','r')
set(gca,'units','centimeters','position',[3.25 2 1.5 2],'xtick',[1,13]);
fp.FormatAxes(gca)
%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','Validations',savedir,1); close all

%% Plot the Reorderd Correlation Matrix
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];

order = NaN(1,numel(idx_of_final_motifs));
for i = 1:numel(idx_of_final_motifs)
   order(i) = find(motif_list==idx_of_final_motifs(i));    
end
%add the noise clusters onto the end
order = [order, find(motif_list==0)];

Plot_OrderedSimilarityMatrix(tcorr_mat,cluster_idx,'order',order,'ytextoffset',-500,'xtextoffset',-500,'fig_position',[680   558   560   420],'range',[0.2 0.8]);
title(sprintf('%d,',idx_of_final_motifs))
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','CorrMatrix_Reordered',savedir,1); close all





