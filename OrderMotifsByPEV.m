function idx = OrderMotifsByPEV(data_dir,group)


%load mouse info
file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
data_all = cellfun(@(x) load(x,'stats_refit'),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.stats_refit,data_all,'UniformOutput',0);
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});

grp = isVPA(mouse);

%average by mouse       
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data = cat(1,data_mouse{:});

grp = isVPA(unique(mouse));

if group ==0 %order by saline
    [~,idx] = sort(nanmean(data(grp==0,:),1),'descend');
elseif group ==1 %order by VPA
    [~,idx] = sort(nanmean(data(grp==1,:),1),'descend');
elseif group ==2 %order by both
    [~,idx] = sort(nanmean(data,1),'descend');
else 
    error('unknown order type');
end

