%Supplemental Figure 2 % Load the original testing and training data (so just 'reverse'); 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\SupplementalFig2';
if ~exist(savedir)
    mkdir(savedir);
end

motif_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
gp = general_params_vpa;
file_list = GrabFiles('\w*reverse\w*.mat', 0, {motif_dir});
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);
grp = isVPA(unique(mouse)); 

%loop through each file and compute statistics
st_var = NaN(1,numel(file_list));
active_pxl = NaN(1,numel(file_list));
for i = 1:numel(file_list)
   if mod(i,10)==0; fprintf('\t\n Working on file %d of %d',i,numel(file_list)); end
   temp = load(file_list{i},'data_test','data_train');   
   %number of active pixels
   active_pxl(i) = size(temp.data_test,1);   
   %spatiotemporal variance across entire time
   st_var(i) = nanmean([nanvar(temp.data_test(:)),nanvar(temp.data_train(:))]);
end

tempact = active_pxl;
tempst = st_var; 

%average per mouse
st_var = arrayfun(@(x) nanmean(st_var(mouse==x)), unique(mouse),'UniformOutput',0);
st_var = cat(1,st_var{:});
active_pxl = arrayfun(@(x) nanmean(active_pxl(mouse==x)), unique(mouse),'UniformOutput',0);
active_pxl = cat(1,active_pxl{:});


%% Compare spatiotemporal variance between groups
fp = fig_params_vpa; 
figure; hold on; 
pos = rand(20,1)/2-0.25;
b = bar(1,nanmean(st_var(grp==0)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_sal;
% errorbar(1,nanmean(st_var(grp==1)),sem(st_var(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==0)+1,st_var(grp==0),'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
b = bar(2,nanmean(st_var(grp==1)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_vpa;
% errorbar(2,nanmean(st_var(grp==1)),sem(st_var(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==1)+2,st_var(grp==1),'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)    
pval = ranksum(st_var(grp==0),st_var(grp==1));

%format
title(sprintf('st var pval %0.2g',pval));
set(gca,'xlim',[0.5, 2.5],'units','centimeters','position',[3.25 2 2 4])
fp.FormatAxes(gca)



%% Compare number of active pixels between groups
fp = fig_params_vpa; 
figure; hold on; 
pos = rand(20,1)/2-0.25;
b = bar(1,nanmean(active_pxl(grp==0)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_sal;
% errorbar(1,nanmean(active_pxl(grp==1)),sem(active_pxl(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==0)+1,active_pxl(grp==0),'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
b = bar(2,nanmean(active_pxl(grp==1)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_vpa;
% errorbar(2,nanmean(active_pxl(grp==1)),sem(active_pxl(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==1)+2,active_pxl(grp==1),'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)    
pval = ranksum(active_pxl(grp==0),active_pxl(grp==1));

%format
title(sprintf('active pixels pval %0.2g',pval));
set(gca,'xlim',[0.5, 2.5],'units','centimeters','position',[3.25 2 2 4],'ylim',[1500 2000])
fp.FormatAxes(gca)


%% Load the histology data 
fp = fig_params_vpa; 
uiopen('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Density Bar Graph.fig',1)

set(gca,'units','centimeters','position',[3.25 2 4 4])
fp.FormatAxes(gca)


handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','OverallNeuroVPASAL',savedir,1); close all

%% Screen for epileptiform Events
%Load Preprocessed Data
fp = fig_params_vpa; 
dff_fn = load('AllOriginalDataFileList.mat','file_list');
dff_fn = dff_fn.file_list;
[~, fn] = cellfun(@(x) fileparts(x), dff_fn, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
mouse_uniq = unique(mouse); 
%just get the first recording from each mouse
for i = 1:numel(mouse_uniq)
   if mouse_uniq(i)==9027 || mouse_uniq(i)==9029 %avoid short recs not used in analysis
       idx = find(mouse == mouse_uniq(i), 2,'first');
       idx = idx(2);
   else
       idx = find(mouse == mouse_uniq(i), 1,'first');
   end
   temp = load(dff_fn{idx},'dff');
   %choose a new on if too short (<12min) recording
   if size(temp.dff,3) < gp.fps*11*60; fprintf('rec for %d too short',mouse_uniq(i)); 
   else
       stack = SpatialBin(temp.dff,2,[],1);
       stack(repmat(nanvar(stack,[],3)<=eps,[1,1,size(stack,3)])==1)=NaN;
       %mask with conservative mask
       temp = load('FigureMask.mat');
       mask = repmat(temp.mask,[1,1,size(stack,3)]);
       stack(mask==0)=NaN;   
       signal = nanmean(conditionDffMat(stack),2);
       [~,~,widths,prominences] = findpeaks(signal,13); 
       figure; plot(widths,prominences,'.','markersize',5,'color','k');
       xlim([0 0.75])
       ylim([0 10]);
       title(sprintf('Mouse %d', mouse_uniq(i)),'FontSize',16,'Fontweight','normal','FontName','Arial');
       xlabel('peak width (s)');
       ylabel('prominence (df/f)');
       fp.FormatAxes(gca)
   end
end

%%
handles = get(groot, 'Children');
for i = 1:numel(handles)
    set(0, 'CurrentFigure', handles(i))
    set(gca,'units','centimeters','position',[3.25 2 4 4],'xtick',[0, 0.25, 0.5, 0.75])
    rectangle(gca,'position',[0.05,5,0.15,4]);
end
fp.SaveFigs(handles,'-svg','DictalScreen',savedir,1); close all
%% compare variance in social recordings

motif_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social';
gp = general_params_vpa;
file_list = GrabFiles('\w*refit\w*.mat', 0, {motif_dir});
mouse = cellfun(@(x) load(x,'mouseID'), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) x.mouseID, mouse, 'UniformOutput',1); 
group = isVPA(mouse);
grp = isVPA(unique(mouse)); 

%loop through each file and compute statistics
st_var = NaN(1,numel(file_list));
active_pxl = NaN(1,numel(file_list));
for i = 1:numel(file_list)
   if mod(i,10)==0; fprintf('\t\n Working on file %d of %d',i,numel(file_list)); end
   temp = load(file_list{i},'data');   
   %number of active pixels
   active_pxl(i) = size(temp.data,1);   
   %spatiotemporal variance across entire time
   st_var_temp=[];
   for j = 1:size(temp.data,3)
       data = temp.data(:,:,j);
       st_var_temp(j) = nanvar(data(:));
   end
   st_var(i) = nanmean(st_var_temp);
end

tempact = active_pxl;
tempst = st_var; 

%average per mouse
st_var = arrayfun(@(x) nanmean(st_var(mouse==x)), unique(mouse),'UniformOutput',0);
st_var = cat(1,st_var{:});
active_pxl = arrayfun(@(x) nanmean(active_pxl(mouse==x)), unique(mouse),'UniformOutput',0);
active_pxl = cat(1,active_pxl{:});


%% Compare spatiotemporal variance between groups
fp = fig_params_vpa; 
figure; hold on; 
pos = rand(20,1)/2-0.25;
b = bar(1,nanmean(st_var(grp==0)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_sal;
% errorbar(1,nanmean(st_var(grp==1)),sem(st_var(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==0)+1,st_var(grp==0),'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
b = bar(2,nanmean(st_var(grp==1)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_vpa;
% errorbar(2,nanmean(st_var(grp==1)),sem(st_var(grp==1)),'LineWidth',1.25,'Color','k');
plot(pos(grp==1)+2,st_var(grp==1),'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)    
pval = ranksum(st_var(grp==0),st_var(grp==1));

%format
title(sprintf('SOCIAL st var pval %0.2g',pval));
set(gca,'xlim',[0.5, 2.5],'units','centimeters','position',[3.25 2 2 4])
fp.FormatAxes(gca)

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','SocialSTVariance',savedir,1); close all

















