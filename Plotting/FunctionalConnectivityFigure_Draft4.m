%DynamicsFigure
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'))
addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'))
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis'))
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\'; 
save_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\FunctionalConnectivity';

if exist([save_dir filesep 'statsdiary.txt'],'file')
    delete([save_dir filesep 'statsdiary.txt']);
end
diary([save_dir filesep 'statsdiary.txt']);
diary on

%get the file names of everything (you are only loading the 'train data' so
%for to get all dataset you want both '_chunk' and '_reverse_chunk'.
fn = GrabFiles('chunk',0,{data_dir});
fp = fig_params_vpa;
close all;

% load phenotype and neurotype data
[pt,~,~,~,~,~,~,betas_pt,~,~,bias_pt] = Phenotype([data_dir,'\ALL'],'verbose',0,'num_classifiers',5,'holdout',0.5);
[nt,~, betas_nt, ~,allstats] = Neurotype([data_dir,'\ALL']);

%% Place ROIs in a grid across the cortex and make figure
temp = load(fn{101},'data_train','nanpxs');
temp = nanmax(conditionDffMat(temp.data_train',temp.nanpxs),[],3);
[x,y] = size(temp);
z=6;
roi = combvec(1+4:z:x,1:z:y)';
%remove nans
badidx = zeros(size(roi,1),1);
for i = 1:size(roi,1)
   if isnan(temp(roi(i,2),roi(i,1)))
       badidx(i)=1;
   end
end
roi(badidx==1,:)=[];
roi(15,:)=[]; %remove the extra point on right side
close; 
imagesc(temp)
hold on
plot(roi(:,1),roi(:,2),'.','color','r','MarkerSize',20)

%add lines showing the distance bins
% dist_bins = 1:10:51;
dist_bins = 0:11:44;
y = 1:1:5;
for i = 1:numel(dist_bins)-1
   plot([y(i),y(i)],[9,dist_bins(i+1)+9],'linewidth',2,'color','r')
end

saveCurFigs(gcf,'-dsvg','ROI',save_dir,0); close all
%% Grab the FC matrix for each epoch
data_fn = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data.mat';
if ~exist(data_fn,'file') %recompile
    fprintf('Recompiling functional connectivity data')
    fc_mat = zeros(size(roi,1),size(roi,1),numel(fn));
    for i = 1:numel(fn) %epoch loop
        temp = load(fn{i},'data_train','nanpxs');
        dff = conditionDffMat(temp.data_train',temp.nanpxs);
        [x,y,z] = size(dff);
        dff = reshape(dff,[x*y,z]);
        roi_trace = zeros(z,size(roi,1));    
        for j = 1:size(roi,1) %roi loop 
            mask = zeros(x,y);
    %         mask(roi{j}.position(1):roi{j}.position(1)+roi{j}.position(3),roi{j}.position(2):roi{j}.position(2)+roi{j}.position(4))=1;
            mask(roi(j,2)-1:roi(j,2)+1,roi(j,1)-1:roi(j,1)+1)=1;
    %         imagesc(mask); pause();
            mask = reshape(mask,[x*y,1]);
            roi_trace(:,j) = nanmean(dff(mask==1,:));        
        end
        fc_mat(:,:,i) = corr(roi_trace-nanmean(roi_trace));
    end
    save(data_fn,'fc_mat','roi','z');
else
    load(data_fn,'fc_mat','roi','z');
end


%% Average FC matrix across epochs per animal
rng('default')
mouse_num = MouseNumFromPath(fn,'Mouse-');
mice = unique(mouse_num);
grp = isVPA(mice);
avg_fc = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
for i = 1:numel(mice)
   idx = find(mouse_num == mice(i));   
   avg_fc(:,:,i) = nanmean(fisherZ(fc_mat(:,:,idx)),3);
end

%remove autocorr and cross hemi
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
avg_fc(repmat(mask,1,1,numel(mice))==1)=NaN;

%get the global difference between groups
all = squeeze(nanmean(avg_fc,[1,2]));
perm_stat = NaN(1,1000);
for j = 1:1000
    temp_group = grp(randperm(numel(grp)));
    perm_stat(j) = nanmean(all(temp_group==0))-nanmean(all(temp_group==1));
end
true_val = nanmean(all(grp==0))-nanmean(all(grp==1));
pval = sum([perm_stat,true_val]>=true_val)/numel([perm_stat,true_val]);

%get the distance between each roi
d_mat = pdist2(roi,roi);
d_mat(mask==1)=NaN;

%plot the correlation by distance for each group
% dist_bins = 1:10:51;
dist_bins = 0:11:44;
d_rho = cell(numel(mice),1);
for i = 1:numel(mice)
    temp = avg_fc(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    uniqd = unique(d);
    for j = 1:numel(uniqd)
        d_rho{i}(j) = nanmean(temp(d==uniqd(j)));
    end    
end
d_rho = cat(1,d_rho{:});
% bar plot version
vpa = d_rho(grp==1,:);
sal = d_rho(grp==0,:);

%statistics
fprintf('\nVPA functional connectivity by distance:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL functional connectivity by distance:')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n Distance between vpa and sal:')
nanmean(sal)-nanmean(vpa)

figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1:numel(dist_bins)-1,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5:numel(dist_bins)-0.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:numel(mice)
    if grp(i)==1
        x = (1.5:numel(dist_bins)-0.5)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_vpa,0.25])
    else
        x = (1:numel(dist_bins)-1)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_sal,0.25])
    end
end
ylim([0.5 1.75])
%add stats
for i = 1:numel(dist_bins)-1
    true_val = nanmean(sal(:,i))-nanmean(vpa(:,i)); 
    val_perm=[];
    for j = 1:1000
       grp_shuf = grp(randperm(20));       
       val_perm(j) = nanmean(d_rho(grp_shuf==0,i))-nanmean(d_rho(grp_shuf==1,i));
    end
    pval = sum(abs([val_perm,true_val])>=abs(true_val))/1001; 
    yval = nanmax(cat(1,vpa(:,i),sal(:,i)))+0.05;
    line([i,i+0.5],[yval,yval],'linewidth',1.5,'color',[0.25 0.25 0.25]);
%     text(i+0.25,yval,sprintf('p=%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom')
    text(i+0.25,yval,sprintf('%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',14)
end
set(gca,'xtick',1.25:4.25,'xticklabels',dist_bins(2:end),'ytick',[0.5:0.25:1.75])
xlim([0.70 4.8])
xlabel('Distance (pixels)')
ylabel('Connectivity (\rho_z)','Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.4 3.5 5],[])
fp.FormatAxes(gca);

% also plot a histogram of the distribution of distances
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
%get the distance between each roi, excluding cross hemispheres
d_mat = pdist2(roi,roi);
d_mat(mask==1)=[];
d = histcounts(d_mat,dist_bins); 
figure;  hold on;
bar([1,2],d(1:2),'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.5);
bar([3,4],d(3:4),'FaceColor',[0.05 0.05 0.05],'FaceAlpha',0.5);
set(gca,'xtick',1:5,'xticklabels',dist_bins(2:end))
xlim([0.25 4.75])
ylim([0 300])
xlabel('Distance (pixels)')
ylabel('# Connections','Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.5 3 1.5],[])
fp.FormatAxes(gca);

handles = get(groot, 'Children');
saveCurFigs(handles,'-dsvg','FunctionalConnectivity',save_dir,0); close all

%% Compare connectivity with neurotype, phenotype
% grp = isVPA(mice);
close all
L_type = {'long','short'};
for cur_type = 1:2
    rng('default')
    d_mat = pdist2(roi,roi);
    d_mat(mask==1)=NaN;
    fp = fig_params_vpa;
    d = discretize(d_mat,dist_bins);    
    if cur_type==1 %long
        d(d<3)=NaN;
        d(~isnan(d))=1;
    else
        d(d>2)=NaN;
        d(~isnan(d))=1;
    end
    %mask the fc
    temp = avg_fc;
    temp(repmat(isnan(d),1,1,numel(mice)))=NaN;
    connectivity = nanmean(reshape(temp,size(roi,1)^2,numel(mice)));        

    %plot the correlation between connectivity and Phenotype
    figure; hold on; 
    lm = fitlm(pt,connectivity); 
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit   
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper 
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(pt,connectivity,'.','markersize',30,'color',[0.5 0.5 0.5]) 
    legend off
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
    ylabel('Connectivity (\rho_z)','Interpreter','tex')
    xlabel('Phenotype')    
    set(gca,'units','centimeters','position',[3.25 2 4 5]);
    rho = corr(pt',connectivity');
    rho_perm=[];
    for j = 1:1000
       temp_conn = connectivity(randperm(20));
       rho_perm(j) = corr(pt',temp_conn');
    end
    %extend axes to limits
    if cur_type==1
        ylim([0.5 1.15])
        xlim([-3 3])
    end
        
%     ylim([floor(min(connectivity)*100)/100,ceil(max(connectivity)*100)/100])
    %given the previous plot, we would expect a negative correlation 
    pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);    
    fp.SetTitle(gca,[L_type{cur_type},sprintf('rho=%0.2g pval=%0.2g',rho,pval)]);
    fp.FormatAxes(gca)    


    figure; hold on; 
    lm = fitlm(nt,connectivity); 
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit   
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper 
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(nt,connectivity,'.','markersize',20,'color',[0.5 0.5 0.5]) 
    legend off
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
    ylabel('Connectivity')
    xlabel('Neurotype')    
    set(gca,'units','centimeters','position',[3.25 2 2.25 3.75]);
    rho = corr(nt',connectivity');
    rho_perm=[];
    for j = 1:1000
       temp_conn = connectivity(randperm(20));
       rho_perm(j) = corr(nt',temp_conn');
    end
    pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);
    if cur_type==1
        ylim([0.5 1.15])
        xlim([-2.25 2.25])
    end
%     ylim([floor(min(connectivity)*100)/100,ceil(max(connectivity)*100)/100])
    fp.SetTitle(gca,[L_type{cur_type},sprintf('rho=%0.2g pval=%0.2g',rho,pval)]);
    fp.FormatAxes(gca)        

    % plot the partial correlation between connectivity and Phenotype
    %plot the partial correlation
    lm = fitlm(nt,connectivity);
    pConnectivity = lm.Residuals.Raw;
%     lm = fitlm(nt,pt);
%     pPT = lm.Residuals.Raw;
    pPT = pt'; 
    figure; hold on; 
    lm = fitlm(pPT,pConnectivity); 
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit   
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper 
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(pPT,pConnectivity,'.','markersize',30,'color',[0.5 0.5 0.5]) 
    legend off
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
    ylabel('Residual Connectivity')
    xlabel('Phenotype')    
    set(gca,'units','centimeters','position',[3.25 2 4 5]);
    rho = corr(pPT,pConnectivity);
    rho_perm=[];
    for j = 1:1000
       temp_conn = pConnectivity(randperm(20));
       rho_perm(j) = corr(pPT,temp_conn);
    end
    pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);
    if cur_type==1
        ylim([-0.16 0.16])
        xlim([-2.75 2.75])
    end  
%     ylim([-0.2 .2])
%     ylim([floor(min(pConnectivity)*100)/100,ceil(max(pConnectivity)*100)/100])
%     xlim([floor(min(pPT)*100)/100,ceil(max(pPT)*100)/100])
    fp.SetTitle(gca,[L_type{cur_type},sprintf('rho=%0.2g pval=%0.2g',rho,pval)]);
    fp.FormatAxes(gca)           
    
    % plot the partial correlation between connectivity and neurotype
    %plot the partial correlation
    lm = fitlm(pt,connectivity);
    pConnectivity = lm.Residuals.Raw;
    lm = fitlm(pt,nt);
    pNT = lm.Residuals.Raw;
    figure; hold on; 
    lm = fitlm(pNT,pConnectivity); 
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit   
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper 
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(pNT,pConnectivity,'.','markersize',30,'color',[0.5 0.5 0.5]) 
    legend off
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
    ylabel('Residual Connectivity')
    xlabel('Residual Neurotype')    
    set(gca,'units','centimeters','position',[3.25 2 4 5]);
    rho = corr(pNT,pConnectivity);
    rho_perm=[];
    for j = 1:1000
       temp_conn = pConnectivity(randperm(20));
       rho_perm(j) = corr(pNT,temp_conn);
    end
    if cur_type==1
        ylim([-0.22 0.225])
        xlim([-1.5 2.])
    end    
    pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);
    fp.SetTitle(gca,[L_type{cur_type},sprintf('rho=%0.2g pval=%0.2g',rho,pval)]);
%     ylim([floor(min(pConnectivity)*100)/100,ceil(max(pConnectivity)*100)/100])
%     xlim([floor(min(pNT)*100)/100,ceil(max(pNT)*100)/100])
    fp.FormatAxes(gca)   

    handles = get(groot, 'Children');
    saveCurFigs(handles,'-dsvg',[L_type{cur_type},'FC_vs_Neuro_Pheno'],save_dir,0); close all
end %long/short loop

%% Compare original vs motif reconstructed connectivity
data_fn = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data_reconstructions.mat';
if ~exist(data_fn,'file') %recompile
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    fn = GrabFiles('chunk',0,{[data_dir,'\ALL']});
    fc_mat = zeros(size(roi,1),size(roi,1),numel(fn));
    for i = 1:numel(fn) %epoch loop
        %reconstruct data
        temp = load(fn{i},'w','H','bad_pxl');
        dff = NaN(numel(temp.bad_pxl),size(temp.H,2));    
        dff(temp.bad_pxl==0,:) = tensor_convolve(temp.w,temp.H);
        [x,z] = size(dff);
        x = sqrt(x);
        y = x; 
        roi_trace = zeros(z,size(roi,1));    
        for j = 1:size(roi,1) %roi loop 
            mask = zeros(x,y);
            mask(roi(j,2)-1:roi(j,2)+1,roi(j,1)-1:roi(j,1)+1)=1;
    %         imagesc(mask); pause();
            mask = reshape(mask,[x*y,1]);
            roi_trace(:,j) = nanmean(dff(mask==1,:));        
        end
        fc_mat(:,:,i) = corr(roi_trace-nanmean(roi_trace));
    end
    save(data_fn,'fc_mat','roi','z');
else
    load(data_fn,'fc_mat','roi','z');
end


%% Reconstruction
rng('default')
dist_bins = 0:11:44;
fn = GrabFiles('chunk',0,{[data_dir,'\ALL']});
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data_reconstructions.mat','fc_mat','roi','z');
%%REPEAT ON ORIG
mouse_num = MouseNumFromPath(fn,'Mouse-');
mice = unique(mouse_num);
avg_fc = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
for i = 1:numel(mice)
   idx = find(mouse_num == mice(i));   
   avg_fc(:,:,i) = nanmean(fisherZ(fc_mat(:,:,idx)),3);
end

%remove autocorr and cross hemi
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
avg_fc(repmat(mask,1,1,numel(mice))==1)=NaN;

%get the distance between each roi
d_mat = pdist2(roi,roi);
d_mat(mask==1)=NaN;

%plot the correlation by distance for each group
d_rho = cell(numel(mice),1);
for i = 1:numel(mice)
    temp = avg_fc(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    uniqd = unique(d);
    for j = 1:numel(uniqd)
        d_rho{i}(j) = nanmean(temp(d==uniqd(j)));
    end    
end
d_rho = cat(1,d_rho{:});
% bar plot version
vpa = d_rho(grp==1,:);
sal = d_rho(grp==0,:);

%statistics
fprintf('\nVPA Reconstruction functional connectivity by distance:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL Reconstructionfunctional connectivity by distance:')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n Distance Reconstruction between vpa and sal:')
nanmean(sal)-nanmean(vpa)

figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1:numel(dist_bins)-1,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5:numel(dist_bins)-0.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:numel(mice)
    if grp(i)==1
        x = (1.5:numel(dist_bins)-0.5)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_vpa,0.25])
    else
        x = (1:numel(dist_bins)-1)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_sal,0.25])
    end
end
ylim([0.5 2.5])
%add stats
for i = 1:numel(dist_bins)-1
    true_val = nanmean(sal(:,i))-nanmean(vpa(:,i)); 
    val_perm=[];
    for j = 1:1000
       grp_shuf = grp(randperm(20));       
       val_perm(j) = nanmean(d_rho(grp_shuf==0,i))-nanmean(d_rho(grp_shuf==1,i));
    end
    %based on the previous plot, we want to test if vpa<sal so tail
    pval = sum([val_perm,true_val]>=(true_val))/numel([val_perm,true_val]); 
    yval = nanmax(cat(1,vpa(:,i),sal(:,i)))+0.05;
    line([i,i+0.5],[yval,yval],'linewidth',1.5,'color',[0.25 0.25 0.25]);
%     text(i+0.25,yval,sprintf('p=%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom')
    text(i+0.25,yval,sprintf('%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',14)
end
set(gca,'xtick',1.25:numel(dist_bins)-1+.25,'xticklabels',dist_bins(2:end),'ytick',[0.5:0.5:2.5])
xlim([0.70 4.8])
xlabel('Distance (um)')
ylabel('Connectivity (\rho_z)','Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.4 5 3.75],[])
fp.FormatAxes(gca);
title('reconstruction','fontweight','normal');

handles = get(groot, 'Children');
saveCurFigs(handles,'-dsvg','FC_Reconstructions',save_dir,0); close all

%% look at dynamicness of motifs and long vs short functional connectivity
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL\basismotifs.mat','W_basis','noise_clusters');
W_basis(:,noise_clusters,:)=[];
[betas_sorted, idx_betas] = sort(betas_nt,'ascend');
W_basis=W_basis(:,idx_betas,:); 
n_mot = size(W_basis,2);
motif_dyn = NaN(1,n_mot);
avg_fc = zeros(size(roi,1),size(roi,1),n_mot);
for i =1:n_mot    
    w=squeeze(W_basis(:,i,:));
    %to avoid spurious inflation of distance in motifs with little activity and then a lot of activity 
    %we only want frames with activity above a threshold. Also, we don't
    %want to have this threshold to be depending on the size of the burst.
    %So take just the top 10 pixels and then keep the frames in which the
    %average activity is at least 10% of the max average activity.
    idx_top = nanmean(maxk(w,20)); 
    idx = [find(idx_top>(0.05*max(idx_top)),1,'first'),find(idx_top>(0.5*max(idx_top)),1,'last')];
    w = w(:,idx(1):idx(2));
    [x,z] = size(w); x = sqrt(x); y = x; 
    
    %FC 
    roi_trace = zeros(z,size(roi,1));    
    for j = 1:size(roi,1) %roi loop 
        mask = zeros(x,y);
        mask(roi(j,2)-1:roi(j,2)+1,roi(j,1)-1:roi(j,1)+1)=1;
        mask = reshape(mask,[x*y,1]);
        roi_trace(:,j) = nanmean(w(mask==1,:));        
    end    
    avg_fc(:,:,i) = fisherZ(corr(roi_trace-nanmean(roi_trace)));     

    %for dynamicness we care about the relative change within the motif... so normalize each frame to 0-->1    
    temp = normalize(w,1,'range');
    e = pdist2(temp',temp');
    mask = tril(ones(size(e)),-1);
    e(mask==0)=[];
    motif_dyn(i) = nanmean(e(:));     
end
%remove autocorr and cross hemi
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
avg_fc(repmat(mask,1,1,n_mot)==1)=NaN;

%get the distance between each roi
d_mat = pdist2(roi,roi);
d_mat(mask==1)=NaN;
fp = fig_params_vpa;

%long
d = discretize(d_mat,dist_bins);
d(d<3)=NaN;
d(~isnan(d))=1;
%mask the fc
long_fc = avg_fc;
long_fc(repmat(isnan(d),1,1,n_mot))=NaN;
long_conn= nanmean(reshape(long_fc,size(roi,1)^2,n_mot));

%total
d = discretize(d_mat,dist_bins);
d(d>2)=NaN;
d(~isnan(d))=1;
%mask the fc
short_fc = avg_fc;
short_fc(repmat(isnan(d),1,1,n_mot))=NaN;
short_conn= nanmean(reshape(short_fc,size(roi,1)^2,n_mot));
motif_connectivity = short_conn-long_conn;
   
%% Plot the motif dynamicsness for supplemental
% r = [12,8,16,9,11,15,4,6,1,14,13,3,2,5,7,10];    %(very) approximate manual ordering for original gut checking 
close all
rng('default')
figure; 
[~,a]=sort(motif_dyn,'ascend');
plot(motif_dyn(a),'linestyle',':','marker','.','markersize',20)
text([1:16]-0.25,motif_dyn(a)+0.75,string(a),'FontSize',14)
xlim([0.5 16.5]);
ylabel('Dynamicness')
set(gca,'xtick',[1,16],'XTickLabel',{'Less Dynamic','More Dynamic'});
fp.FigureSizing(gcf,[2.5 2.5 8 4],[])
fp.FormatAxes(gca);
title('Motif Dynamics','fontweight','normal');

figure; 
[~,a]=sort(motif_connectivity,'ascend');
plot(motif_connectivity(a),'linestyle',':','marker','.','markersize',20)
text([1:16]-0.25,motif_connectivity(a)+0.1,string(a),'FontSize',14)
xlim([0.5 16.5]);
ylabel('Long vs Short Connectivity')
set(gca,'xtick',[1,16],'XTickLabel',{'More Disconnected','Less Disconnected'});
fp.FigureSizing(gcf,[2.5 2.5 8 4],[])
fp.FormatAxes(gca);
title('Motif Dynamics','fontweight','normal');

% Plot the dynamicness as a function of betas
figure; hold on; 
lm = fitlm(betas_sorted',motif_dyn'); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(betas_sorted,motif_dyn,'.','markersize',30,'color',[0.5 0.5 0.5]) 
legend off
xlabel('Betas')
ylabel('Dynamicness')
set(gca,'units','centimeters','position',[3.25 2 4 5]);
rho = corr(betas_sorted,motif_dyn');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_dyn(randperm(16));
   rho_perm(j) = corr(betas_sorted,temp_conn');
end
ylim([0 12])
xlim([-0.6 .9])
pval = sum([rho_perm,rho]>=rho)/numel([rho_perm,rho]);
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca)  

% Plot the dynamicness by group
motif_group = betas_sorted>0; %positive betas = vpa
vpa = motif_dyn(motif_group==1);
sal = motif_dyn(motif_group==0);

%statistics
fprintf('\nVPA Dynamicness:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL Dynamicness :')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n Dynamicness between sal vpa:')
nanmean(sal)-nanmean(vpa)


figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:n_mot
    if motif_group(i)==1
        x = (1.5)+rand(1)/4-0.125;
        plot(x,motif_dyn(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
    else
        x = (1)+rand(1)/4-0.125;
        plot(x,motif_dyn(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
    end
end
perm_stat = NaN(1,1000);
for i = 1:1000
   temp_group = motif_group(randperm(numel(motif_group)));
   perm_stat(i) = nanmean(motif_dyn(temp_group==1))-nanmean(motif_dyn(temp_group==0));
end
%given the positive correlation, we would expect a tailed
true_val = nanmean(vpa)-nanmean(sal);
pval = sum([perm_stat,true_val]>=true_val)/numel([perm_stat,true_val]);
xlim([0.5 2]);
ylabel('Dynamicness')
fp.FigureSizing(gcf,[2.5 2.5 2 4],[])
fp.SetTitle(gca,sprintf('pval=%0.2g',pval));
fp.FormatAxes(gca);

% Plot the disconnectivity vs dynamicness
figure; hold on; 
lm = fitlm(motif_connectivity,motif_dyn); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(motif_connectivity,motif_dyn,'.','markersize',30,'color',[0.5 0.5 0.5]) 
text(motif_connectivity,motif_dyn,string(1:16),'FontSize',14)
legend off
ylabel('Dynamicness')
xlabel({'Long-Short','Disconnectivity'},'Interpreter','tex')
set(gca,'units','centimeters','position',[3.25 2 4 5]);
rho = corr(motif_connectivity',motif_dyn');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_connectivity(randperm(16));
   rho_perm(j) = corr(temp_conn',motif_dyn');
end
%testing if negatively correlated 
pval = sum([rho_perm,rho]>=rho)/numel([rho_perm,rho]);
ylim([0 12])
xlim([-0.1,1])
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca)  

% Plot the motif_connectivity as a function of betas
figure; hold on; 
lm = fitlm(betas_sorted,motif_connectivity); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(betas_sorted,motif_connectivity,'.','markersize',30,'color',[0.5 0.5 0.5]) 
legend off
xlabel('Betas (\rho_z)','Interpreter','tex')
ylabel({'Long-Short','Disconnectivity (\rho_z)'},'Interpreter','tex')
set(gca,'units','centimeters','position',[3.25 2 2.25 3.75]);
rho = corr(betas_sorted,motif_connectivity');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_connectivity(randperm(16));
   rho_perm(j) = corr(betas_sorted,temp_conn');
end
ylim([-1.25 0.25])
xlim([-0.6 .9])
% ylim([floor(min(motif_connectivity)*100)/100,ceil(max(motif_connectivity)*100)/100])
% xlim([floor(min(betas_sorted)*100)/100,ceil(max(betas_sorted)*100)/100])
pval = sum([rho_perm,rho]>=rho)/numel([rho_perm,rho]);
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca) 

%statistics
fprintf('\nVPA motif fc:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL motif fc :')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n motif fc between sal vpa:')
nanmean(sal)-nanmean(vpa)

% Plot the dynamicness by group
motif_group = betas_sorted>0; %positive betas = vpa
vpa = motif_connectivity(motif_group==1);
sal = motif_connectivity(motif_group==0);

figure; hold on; 
line([0.5 2],[0,0],'color','k','linewidth',1.5)
%do a bar plot for each distance with connected scatters
bar(1,nanmean(sal),'facecolor',[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5,nanmean(vpa),'facecolor',[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:n_mot
    if motif_group(i)==1
        x = (1.5)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color','k')
    else
        x = (1)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color','k')
    end
end
perm_stat = NaN(1,1000);
for i = 1:1000
   temp_group = motif_group(randperm(numel(motif_group)));
   perm_stat(i) = nanmean(motif_connectivity(temp_group==1))-nanmean(motif_connectivity(temp_group==0));
end
%given the positive correlation, we would expect a tailed
true_val = nanmean(vpa)-nanmean(sal);
pval = sum([perm_stat,true_val]>=true_val)/numel([perm_stat,true_val]);
set(gca,'xtick',[1 1.5])
xlim([0.5 2]);
ylim([-0.1,1])
ylabel({'Long-Range','Disconnectivity (\rho_z)'},'Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.5 2 5],[])
fp.SetTitle(gca,sprintf('pval=%0.2g',pval));
fp.FormatAxes(gca);

handles = get(groot, 'Children');
saveCurFigs(handles,'-dsvg','Motif_FCandDynamicness',save_dir,0); close all

%% opposite (long - short)
motif_connectivity = long_conn-short_conn;
close all
rng('default')
figure; 
[~,a]=sort(motif_dyn,'ascend');
plot(motif_dyn(a),'linestyle',':','marker','.','markersize',20)
text([1:16]-0.25,motif_dyn(a)+0.75,string(a),'FontSize',14)
xlim([0.5 16.5]);
ylabel('Dynamicness')
set(gca,'xtick',[1,16],'XTickLabel',{'Less Dynamic','More Dynamic'});
fp.FigureSizing(gcf,[2.5 2.5 8 4],[])
fp.FormatAxes(gca);
title('Motif Dynamics','fontweight','normal');

figure; 
[~,a]=sort(motif_connectivity,'ascend');
plot(motif_connectivity(a),'linestyle',':','marker','.','markersize',20)
text([1:16]-0.25,motif_connectivity(a)+0.1,string(a),'FontSize',14)
xlim([0.5 16.5]);
ylabel('Long vs Short Connectivity')
set(gca,'xtick',[1,16],'XTickLabel',{'More Disconnected','Less Disconnected'});
fp.FigureSizing(gcf,[2.5 2.5 8 4],[])
fp.FormatAxes(gca);
title('Motif Dynamics','fontweight','normal');

% Plot the dynamicness as a function of betas
figure; hold on; 
lm = fitlm(betas_sorted',motif_dyn'); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(betas_sorted,motif_dyn,'x','markersize',10,'color',[0.5 0.5 0.5],'linewidth',2) 
legend off
xlabel('Betas')
ylabel('Dynamicness')
set(gca,'units','centimeters','position',[3.25 2 4 5]);
rho = corr(betas_sorted,motif_dyn');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_dyn(randperm(16));
   rho_perm(j) = corr(betas_sorted,temp_conn');
end
ylim([0 12])
xlim([-0.6 .9])
pval = sum([rho_perm,rho]>=rho)/numel([rho_perm,rho]);
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca)  

% Plot the dynamicness by group
motif_group = betas_sorted>0; %positive betas = vpa
vpa = motif_dyn(motif_group==1);
sal = motif_dyn(motif_group==0);

%statistics
fprintf('\nVPA Dynamicness:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL Dynamicness :')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n Dynamicness between sal vpa:')
nanmean(sal)-nanmean(vpa)


figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:n_mot
    if motif_group(i)==1
        x = (1.5)+rand(1)/4-0.125;
        plot(x,motif_dyn(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
    else
        x = (1)+rand(1)/4-0.125;
        plot(x,motif_dyn(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
    end
end
perm_stat = NaN(1,1000);
for i = 1:1000
   temp_group = motif_group(randperm(numel(motif_group)));
   perm_stat(i) = nanmean(motif_dyn(temp_group==1))-nanmean(motif_dyn(temp_group==0));
end
%given the positive correlation, we would expect a tailed
true_val = nanmean(vpa)-nanmean(sal);
pval = sum([perm_stat,true_val]>=true_val)/numel([perm_stat,true_val]);
xlim([0.5 2]);
ylabel('Dynamicness')
fp.FigureSizing(gcf,[2.5 2.5 2 4],[])
fp.SetTitle(gca,sprintf('pval=%0.2g',pval));
fp.FormatAxes(gca);

% Plot the disconnectivity vs dynamicness
figure; hold on; 
lm = fitlm(motif_connectivity,motif_dyn); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(motif_connectivity,motif_dyn,'x','markersize',10,'color',[0.5 0.5 0.5],'linewidth',2) 
text(motif_connectivity,motif_dyn,string(1:16),'FontSize',14)
legend off
ylabel('Dynamicness')
set(gca,'units','centimeters','position',[3.25 3 4 5]);
rho = corr(motif_connectivity',motif_dyn');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_connectivity(randperm(16));
   rho_perm(j) = corr(temp_conn',motif_dyn');
end
%testing if negatively correlated 
pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);
ylim([0 12])
xlim([-1,.1])
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca)  
xlabel({'Long - Short Range','Connectivity (\rho_z)'},'Interpreter','tex','Fontsize',16)

% Plot the motif_connectivity as a function of betas
figure; hold on; 
line([0 0],[-1.25 0.25],'linestyle',':','color','k','linewidth',2)
lm = fitlm(betas_sorted,motif_connectivity); 
p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %points
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(betas_sorted,motif_connectivity,'x','markersize',10,'color',[0.5 0.5 0.5],'linewidth',2) 
legend off
xlabel('Betas (\rho_z)','Interpreter','tex')
set(gca,'units','centimeters','position',[3.25 2 2.25 3.75]);
rho = corr(betas_sorted,motif_connectivity');
rho_perm=[];
for j = 1:1000
   temp_conn = motif_connectivity(randperm(16));
   rho_perm(j) = corr(betas_sorted,temp_conn');
end
ylim([-1.25 0.25])
xlim([-0.6 .9])
% ylim([floor(min(motif_connectivity)*100)/100,ceil(max(motif_connectivity)*100)/100])
% xlim([floor(min(betas_sorted)*100)/100,ceil(max(betas_sorted)*100)/100])
pval = sum([rho_perm,rho]<=rho)/numel([rho_perm,rho]);
fp.SetTitle(gca,sprintf('rho=%0.2g pval=%0.2g',rho,pval));
fp.FormatAxes(gca) 
ylabel({'Long - Short Range','Connectivity (\rho_z)'},'Interpreter','tex','Fontsize',16)

%statistics
fprintf('\nVPA motif fc:')
nanmean(vpa)
ci = bootci(1000,@nanmean,vpa)

fprintf('\nSAL motif fc :')
nanmean(sal)
ci = bootci(1000,@nanmean,sal)

fprintf('\n motif fc between sal vpa:')
nanmean(sal)-nanmean(vpa)

% Plot the dynamicness by group
motif_group = betas_sorted>0; %positive betas = vpa
vpa = motif_connectivity(motif_group==1);
sal = motif_connectivity(motif_group==0);

figure; hold on; 
line([0.5 2],[0,0],'color','k','linewidth',1.5)
%do a bar plot for each distance with connected scatters
bar(1,nanmean(sal),'facecolor',[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5,nanmean(vpa),'facecolor',[0.5 0.5 0.5],'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:n_mot
    if motif_group(i)==1
        x = (1.5)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color','k')
    else
        x = (1)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color','k')
    end
end
perm_stat = NaN(1,1000);
for i = 1:1000
   temp_group = motif_group(randperm(numel(motif_group)));
   perm_stat(i) = nanmean(motif_connectivity(temp_group==1))-nanmean(motif_connectivity(temp_group==0));
end
%given the positive correlation, we would expect a tailed
true_val = nanmean(vpa)-nanmean(sal);
pval = sum([perm_stat,true_val]<=true_val)/numel([perm_stat,true_val]);
set(gca,'xtick',[1 1.5])
xlim([0.5 2]);
ylim([-1,.11])
ylabel({'Long-Range','Disconnectivity (\rho_z)'},'Interpreter','tex')
fp.FigureSizing(gcf,[4 2.5 2 5],[])
fp.SetTitle(gca,sprintf('pval=%0.2g',pval));
fp.FormatAxes(gca);
%%
handles = get(groot, 'Children');
saveCurFigs(handles,'-dsvg','Long_Minus_Short',save_dir,0); close all

diary off
%%




