%Supplemental Figure 3
%first load data to get the correct indexing of motifs; 
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
data_mouse_nozscore = data_mouse;

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
brainsize = cellfun(@(x) 4624-numel(x.nanpxs),temp,'UniformOutput',1);
brainsize = arrayfun(@(x) nanmean(brainsize(mouse==x)), unique(mouse),'UniformOutput',1);
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
       
       title(sprintf('Motif %d pval=%0.2g',find(idx_of_final_motifs==motif_list(i)),pval)); 
       %format
       set(gca,'xlim',[0.25, 2.75],'units','centimeters','position',[4,3,1,2])
       set(gca,'ylim',[0.2 0.8]);
       fp.FormatAxes(gca)   
       
   end
end

%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','Individualmotifdiff',savedir,1); close all

%% Compare the average motif per animal to the basis motif across animals
N = size(W_basis,2);

%ignoring the noiseclusters
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];
    
%loop through motif clusters
rho = NaN(max(motif_list),numel(unique(motif_id)));
n_contrib = NaN(max(motif_list),numel(unique(motif_id)));
mouse_idx = unique(motif_id);
mouse_avg_motif={};
rng('default')
for i = 1:N
   if sum(ismember(noise_clusters,i))==0
       temp_motifs = W_alligned(:,cluster_idx==i,:); 
       basis_flat = squeeze(W_basis(:,motif_list(i),:));%use W_basis because it's made using the corr community
       basis_flat(nanpxs,:) = [];       
       animal_idx = motif_id(cluster_idx==i);
       
       mouse_list = unique(animal_idx);
%        mouse_avg_motif = NaN(size(W_alligned,1),numel(mouse_list),size(W_alligned,3));
       for cur_mouse = 1:numel(mouse_list)           
          temp_idx = find(mouse_idx == mouse_list(cur_mouse));
          %get that mouse's alligned motifs
          temp = temp_motifs(:,animal_idx==mouse_list(cur_mouse),:);
          temp_flat = zeros(size(temp,1)*size(temp,3),size(temp,2));
          for j = 1:size(temp,2)
             temp_flat(:,j) = reshape(squeeze(temp(:,j,:)),size(temp,1)*size(temp,3),1);
          end
          %get the average correlation to the basis (do it per motif then average, because averaging first may lead to more smoothing for motifs that use a motif more          
          rho(motif_list(i),temp_idx) = fisherInverse(nanmean(fisherZ(corr(temp_flat,basis_flat(:)))));
          [~,max_idx] = min(corr(temp_flat,basis_flat(:)));
          %get the number of motifs contributed by an animal
          n_contrib(motif_list(i),temp_idx) = size(temp_flat,2);
          %get a random representative for each animal. Could also use 'representative' with the mean rho val.
          rand_idx = find(animal_idx==mouse_list(cur_mouse));
          
%           mouse_avg_motif{motif_list(i)}(:,temp_idx,:) = temp_motifs(:,rand_idx(randperm(numel(rand_idx),1)),:);              %random
          mouse_avg_motif{motif_list(i)}(:,temp_idx,:) = temp_motifs(:,rand_idx(max_idx),:);              %max or min
%           temp_flat = mouse_avg_motif(:,temp_idx,:);          
       end                     
   end
end

% convert
rho_z = fisherZ(rho);
%% plotting
group_id = isVPA(mouse_idx);
%loop through motifs
for i = 1:size(rho_z,1)    
    figure('position',[681   216   950   763]); hold on; 
    
    subplot(2,2,1); hold on;
    %plot bar and jitter comparing between groups
    temp_vpa = rho_z(i,group_id==1);
    temp_sal = rho_z(i,group_id==0);    
    pos = rand(20,1)/2-0.25;
    b = bar(1,nanmean(temp_vpa),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
    b.CData(1,:) = fp.c_vpa;
    plot(pos(group_id==1)+1,temp_vpa,'.','markersize',fp.m_markersize,'color',fp.c_vpa)
    b = bar(2,nanmean(temp_sal),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
    b.CData(1,:) = fp.c_sal;
    plot(pos(group_id==0)+2,temp_sal,'.','markersize',fp.m_markersize,'color',fp.c_sal)    
    pval = ranksum(temp_vpa(:),temp_sal(:));

    %format
    fp.FormatAxes(gca)   
    title({'',sprintf('pval %0.2g',pval)},'fontsize',10);
    axis square
%     set(gca,'xlim',[0.5, 2.5],'ylim',[0 2.5],'units','centimeters'); axis square
     

    %compare the fit to the average motif with the PEV
    subplot(2,2,2); hold on; 
    plot(data_mouse_nozscore(:,i),rho_z(i,:)','.','markersize',fp.m_markersize,'color',[0.5 0.5 0.5]);
    [rho,pval] = corr(data_mouse_nozscore(:,i),rho_z(i,:)','rows','complete'); 
    mdl = fitlm(data_mouse_nozscore(:,i),rho_z(i,:)','RobustOpts','on');
    rhorob = mdl.Rsquared.Ordinary; 
    pvalrob = table2array(mdl.Coefficients(2,4));            
    fp.FormatAxes(gca) 
    ylabel('RhoZ','fontsize',10);
    xlabel('PEV','fontsize',10);
    title({'PEV vs RhoZ',sprintf('r=%0.2g p=%0.2g, r-rob=%0.2g p-rob=%0.2g',rho,pval,rhorob,pvalrob)},'fontsize',10);
    axis square
    
    
    %compare the number of motifs in each cluster with the PEV
    subplot(2,2,3); hold on; 
    plot(data_mouse_nozscore(:,i),n_contrib(i,:)','.','markersize',fp.m_markersize,'color',[0.5 0.5 0.5]);
    [rho,pval] = corr(data_mouse_nozscore(:,i),n_contrib(i,:)','rows','complete'); 
    fp.FormatAxes(gca)   
    mdl = fitlm(data_mouse_nozscore(:,i),n_contrib(i,:)','RobustOpts','on');
    rhorob = mdl.Rsquared.Ordinary; 
    pvalrob = table2array(mdl.Coefficients(2,4));    
    title({'PEV vs n_contrib',sprintf('r=%0.2g p=%0.2g, r-rob=%0.2g p-rob=%0.2g',rho,pval,rhorob,pvalrob)},'fontsize',10);
    ylabel('n_contrib','fontsize',10);
    xlabel('PEV','fontsize',10);    
    axis square
      
    
    %get the partial correlation number of contributers
    subplot(2,2,4); hold on; 
    %predict rho_z from n_contribution
    lm = fitlm(n_contrib(i,:)',rho_z(i,:)');
    rY = lm.Residuals.Raw;
    
    %predict data_mouse from n_contrib
    lm = fitlm(n_contrib(i,:)',data_mouse(:,i));
    rX = lm.Residuals.Raw;
    
    %plot the residuals
    plot(rX,rY,'.','markersize',fp.m_markersize,'color',[0.5 0.5 0.5]);
    [rho,pval] = corr(rX,rY,'rows','complete');    
    mdl = fitlm(rX,rY,'RobustOpts','on');
    rhorob = mdl.Rsquared.Ordinary; 
    pvalrob = table2array(mdl.Coefficients(2,4));    
    axis square
    fp.FormatAxes(gca) 
    ylabel('resid RhoZ','fontsize',10);
    xlabel('resid PEV','fontsize',10);       
    title({'pCorr RhoZ vs PEV',sprintf('r=%0.2g p=%0.2g, r-rob=%0.2g p-rob=%0.2g',rho,pval,rhorob,pvalrob)},'fontsize',10);
    
    
    sgtitle(sprintf('Motif %d',find(idx_of_final_motifs==i)));    
end
%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','indi_cluster_corr_vs_basis_PEV',savedir,1); 


%% plot the average or min or max or median motif for each mouse for each motif
 
% %reconstruct
% for cur_motif = 1:numel(mouse_avg_motif)
%     for i = 1:size(mouse_avg_motif{cur_motif},2)
%        temp = zeros(68*68,size(mouse_avg_motif{cur_motif},3));    
%        temp(~ismember(1:size(temp,1),nanpxs),:) = squeeze(mouse_avg_motif{cur_motif}(:,i,:));
%        temp = reshape(temp,[68,68,27]);
%        for j =1:27
%           imagesc(temp(:,:,j),[0 0.05]); title(sprintf('m%d,f%d',i,j)); pause();
%        end
%     end    
% end

%% Force clustering within each motif into two to see how it splits between groups

N = size(W_basis,2);

%ignoring the noiseclusters
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];
    
%loop through motif clusters
mouse_idx = unique(motif_id);
clustID = zeros(16,4);
rng('default')
for i = 1:N
   if sum(ismember(noise_clusters,i))==0
       temp = W_alligned(:,cluster_idx==i,:); 
       temp_flat = zeros(size(temp,1)*size(temp,3),size(temp,2));
       for j = 1:size(temp,2)
          temp_flat(:,j) = reshape(squeeze(temp(:,j,:)),size(temp,1)*size(temp,3),1);
       end
       idx = kmeans(temp_flat',2,'Replicates',10);
       %split animal ids by group   
       animal_idx = motif_id(cluster_idx==i);
       clustID(motif_list(i),:) = [sum(isVPA(unique(animal_idx(idx==1)))==0),sum(isVPA(unique(animal_idx(idx==1)))==1),sum(isVPA(unique(animal_idx(idx==2)))==0),sum(isVPA(unique(animal_idx(idx==2)))==1)];
       
       figure; hold on; 
       b = bar([1,2,3.5,4.5],clustID(motif_list(i),:),'FaceColor','flat','FaceAlpha',0.5,'EdgeColor',[1 1 1]);
       b.CData = [fp.c_sal;fp.c_vpa;fp.c_sal;fp.c_vpa];
       plot([2.75,2.75],[0 11],'linestyle','--','color','k','linewidth',1)
       set(gca,'xtick',[1.5,4],'xticklabel',[1,2],'ytick',[0 11])
       title(sprintf('%d',motif_list(i)),'fontweight','normal');
       fp.FormatAxes(gca)
       set(gca,'xlim',[0.5, 5],'ylim',[0 11],'units','centimeters','position',[3.25 2 1.5 2])
   end
end

% get the ratio
figure; hold on; 
clustRatio = (clustID(:,1)./clustID(:,2)+clustID(:,3)./clustID(:,4))/2;
pos = rand(16,1)/2-0.25;
%plot bar and jitter comparing between groups    
b = bar(1,nanmean(clustRatio),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.25 0.25 0.25];
plot(pos+1,clustRatio,'.','markersize',fp.m_markersize,'color',[0.1 0.1 0.1])
line([0 2],[9/11,9/11],'linewidth',2,'color',[0.25 0.25 0.25],'linestyle',':')
set(gca,'xlim',[0.25 1.75],'xtick','','units','centimeters','position',[4 2 2.5 4]);
ylim([0.5 1.2])
ylabel({'SAL/VPA','Ratio'})
title('Kmeans Ratios'); 
fp.FormatAxes(gca)

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','ForcedClusterSplit',savedir,1); close all

%% Compare statistics about basis motifs relative to the neurotype
N = size(W_basis,2);

%ignoring the noiseclusters
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];
    
%loop through motif clusters
halfpeak_width = NaN(max(motif_list),1);
rho = NaN(max(motif_list),1);
rng('default')
for i = 1:N
   if sum(ismember(noise_clusters,i))==0
       temp_motifs = W_alligned(:,cluster_idx==i,:); 
       basis_flat = squeeze(W_basis(:,motif_list(i),:));%use W_basis because it's made using the corr community
       basis_flat(nanpxs,:) = [];       
       basis_flat_avg = nanmean(basis_flat);
       basis_flat_avg = basis_flat_avg-min(basis_flat_avg); %center at zero
       
       %get motif half peak duration      
       [pks,loc,halfpeak_width(motif_list(i)),~] = findpeaks(basis_flat_avg,'WidthReference','halfheight','Npeaks',1);
%        findpeaks(basis_flat_avg,'Annotate',"extents"); pause();
       halfpeak_width(motif_list(i)) = halfpeak_width(motif_list(i))*75/1000;

       %get spatial dissimilarity within a motif (at each peak)
       rho(motif_list(i)) = corr(basis_flat(:,find(basis_flat_avg>=(pks/2),1,'first')),basis_flat(:,find(basis_flat_avg>=(pks/2),1,'last')));       
   end
end
%reorder to match betas
rho = fisherZ(rho(idx_of_final_motifs));
halfpeak_width = halfpeak_width(idx_of_final_motifs);

%plot the comparison
figure; hold on;
lm = fitlm(betas,rho); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0.25 0.25 0.25]; %fit
p(3).Color = [0.25 0.25 0.25]; %bounds lower
p(4).Color = [0.25 0.25 0.25]; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

plot(betas,rho,'.','markersize',30,'color',[0.25 0.25 0.25])
legend off
set(gca,'xlim',[round(min(betas),1)-0.1,round(max(betas),1)+0.1])
set(gca,'ylim',[round(min(rho),1)-0.1,round(max(rho),1)+0.1]);
ylabel('Internal Spatial Rho')
xlabel('Neurotype Betas')
set(gca,'units','centimeters','position',[3.25 2 6 6]);
[rho_temp,pval_temp] = corr(betas,rho);
title(gca,{'Spatial Change vs Beta',sprintf('rho=%0.2g pval=%0.2g',rho_temp,pval_temp)});
fp.FormatAxes(gca)

%plot the comparison
figure; hold on;
lm = fitlm(betas,halfpeak_width); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0.25 0.25 0.25]; %fit
p(3).Color = [0.25 0.25 0.25]; %bounds lower
p(4).Color = [0.25 0.25 0.25]; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

plot(betas,halfpeak_width,'.','markersize',30,'color',[0.25 0.25 0.25])
legend off
set(gca,'xlim',[round(min(betas),1)-0.1,round(max(betas),1)+0.1])
set(gca,'ylim',[round(min(halfpeak_width),1)-0.1,round(max(halfpeak_width),1)+0.1])
ylabel('Half-Peak Width')
xlabel('Neurotype Betas')
set(gca,'units','centimeters','position',[3.25 2 6 6]);
[rho_temp,pval_temp] = corr(betas,halfpeak_width);
title(gca,{'Motif Width vs Beta',sprintf('rho=%0.2g pval=%0.2g',rho_temp,pval_temp)});
fp.FormatAxes(gca)
%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','BasicMotifStats',savedir,1); close all

%% Compare the number of nanpixels across mouse and groups with neurotype

figure; hold on;
lm = fitlm(nt',brainsize'); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0.25 0.25 0.25]; %fit
p(3).Color = [0.25 0.25 0.25]; %bounds lower
p(4).Color = [0.25 0.25 0.25]; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

plot(nt',brainsize','.','markersize',30,'color',[0.25 0.25 0.25])
legend off
set(gca,'xlim',[round(min(nt),1)-0.1,round(max(nt),1)+0.1])
set(gca,'ylim',[1700,1900])
ylabel('FOV (pixels)')
xlabel('Neurotype')
set(gca,'units','centimeters','position',[3.25 2 6 6]);
[rho_temp,pval_temp] = corr(nt',brainsize');
title(gca,{'FOV vs Neurotype',sprintf('rho=%0.2g pval=%0.2g',rho_temp,pval_temp)});
fp.FormatAxes(gca)

figure; hold on; 
pos = rand(20,1)/2-0.25;
b = bar(1,nanmean(brainsize(grp==0)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_sal;
% errorbar(1,nanmean(trainpev),sem(trainpev(:)),'LineWidth',1.25,'Color','k');
plot(pos(grp==0)+1,brainsize(grp==0),'.','markersize',fp.m_markersize,'color',fp.c_sal)
b = bar(2,nanmean(brainsize(grp==1)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_vpa;
% errorbar(2,nanmean(testpev),sem(testpev(:)),'LineWidth',1.25,'Color','k');
plot(pos(grp==1)+2,brainsize(grp==1),'.','markersize',fp.m_markersize,'color',fp.c_vpa)    

pval = ranksum(brainsize(grp==1),brainsize(grp==0));

%format
title(sprintf('p=%0.2g',pval),'fontweight','normal');
set(gca,'xlim',[0.5, 2.5],'ylim',[1700 1900],'units','centimeters','position',[3.25 2 1.5 2])
set(gca,'ytick',[1700,1900])
set(gca,'xticklabels',{'SAL','VPA'},'XTickLabelRotation',90)
fp.FormatAxes(gca)

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','BrainSize',savedir,1); close all

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

%% get the ratio of saline to vpa per core community 

N = size(W_basis,2);

%ignoring the noiseclusters
motif_list = [1,2,3,4,5,6,0,0,7,8,9,0,0,10,11,12,13,0,0,0,14,0,15,16];
    
%loop through motif clusters
core_comm_ratio = NaN(max(motif_list),1);
rng('default')
for i = 1:N
   if sum(ismember(noise_clusters,i))==0
       temp = isVPA(unique(motif_id(core_comm_idx{i}))); 
       core_comm_ratio(motif_list(i))=sum(temp)/numel(temp);       
   end
end

% get the ratio
figure; hold on; 
pos = rand(16,1)/2-0.25;
%plot bar and jitter comparing between groups    
b = bar(1,nanmean(core_comm_ratio),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = [0.25 0.25 0.25];
plot(pos+1,core_comm_ratio,'.','markersize',fp.m_markersize,'color',[0.1 0.1 0.1])
line([0 2],[11/20,11/20],'linewidth',2,'color',[0.25 0.25 0.25],'linestyle',':')
set(gca,'xlim',[0.25 1.75],'xtick','','units','centimeters','position',[4 2 2.5 4]);
ylim([0 1])
ylabel({'Fraction VPA'})
title('Core Community','fontweight','normal'); 
fp.FormatAxes(gca)


%plot relationship between accuracy and motif beta 
core_comm_ratio = core_comm_ratio(idx_of_final_motifs);
figure; hold on;
lm = fitlm(betas,core_comm_ratio); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0.25 0.25 0.25]; %fit
p(3).Color = [0.25 0.25 0.25]; %bounds lower
p(4).Color = [0.25 0.25 0.25]; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

plot(betas,core_comm_ratio,'.','markersize',30,'color',[0.25 0.25 0.25])
legend off
set(gca,'xlim',[round(min(betas),1)-0.1,round(max(betas),1)+0.1])
set(gca,'ylim',[round(min(core_comm_ratio),1)-0.1,round(max(core_comm_ratio),1)+0.1]);
ylabel('Fraction VPA')
xlabel('Neurotype Betas')
set(gca,'units','centimeters','position',[3.25 2 6 6]);
[rho_temp,pval_temp] = corr(betas,core_comm_ratio(:));
title(gca,{'CoreComm vs Beta',sprintf('rho=%0.2g pval=%0.2g',rho_temp,pval_temp)});
fp.FormatAxes(gca)
%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','CoreCommunityRatio',savedir,1); close all

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





