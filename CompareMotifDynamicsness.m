function CompareMotifDynamicsness(data_dir)
%use varnames for any single unit statistics (e.g. pev)

file_list = GrabFiles('\w*chunk\w*',0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
mouse_grp = isVPA(unique(mouse));

fp = fig_params;
%load data
data_all = cellfun(@(x) load(x,'W'),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.W, data_all,'UniformOutput',0);

%% get the average spatiotemporal variance
data = ones(1,numel(data_all));
for i = 1:numel(data_all)
   temp = data_all{i};
   data(i) = nanmean(arrayfun(@(n) nanvar(temp(:,n,:),[],'all'),1:size(temp,2),'UniformOutput',1));   
end

figure; hold on;    
group_id = [1,0];
for j = 1:numel(group_id)        
    bar(j,nanmean(data(mouse_grp==group_id(j))),'facecolor',[0.75 0.75 0.75],'edgecolor','none')%bar plot    
    x = j-0.25+rand(sum(mouse_grp==group_id(j)),1)/2;
    plot(x,data(mouse_grp==group_id(j)),'.','markersize',20,'color',[0.1 0.1 0.1 0.5]);%jitter plot        
end

%statistics
p=ranksum(data(mouse_grp==1),data(mouse_grp==0));   
fp.FormatAxes(gca); 
title('avg variance')

%% getthe dissimilarity from the mean image

data = ones(1,numel(data_all));
for i = 1:numel(data_all)
   temp = data_all{i};       
   data(i) = fisherInverse(nanmean(arrayfun(@(n) nanmean(fisherZ(corr(squeeze(temp(:,n,:)),nanmean(squeeze(temp(:,n,:)),2)))),1:size(temp,2),'UniformOutput',1)));   
end

figure; hold on;    
group_id = [1,0];
for j = 1:numel(group_id)        
    bar(j,nanmean(data(mouse_grp==group_id(j))),'facecolor',[0.75 0.75 0.75],'edgecolor','none')%bar plot    
    x = j-0.25+rand(sum(mouse_grp==group_id(j)),1)/2;
    plot(x,data(mouse_grp==group_id(j)),'.','markersize',20,'color',[0.1 0.1 0.1 0.5]);%jitter plot        
end

%statistics
p=ranksum(data(mouse_grp==1),data(mouse_grp==0));   
fp.FormatAxes(gca); 
title('avg dissimilarity')





%comapre to the mean of that motif
