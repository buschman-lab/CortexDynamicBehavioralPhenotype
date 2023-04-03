%Supplemental Figure 5; get the explained variance of the longer recordings

rec_names = {'431-10-17-2019','432-10-17-2019','432-10-18-2019'};
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\ExampleData';
fp = fig_params_vpa; 

% Load refitting data
pev = 1:numel(rec_names);
for i = 1:numel(rec_names)
    file_list = GrabFiles(['\w*',rec_names{i},'\w*'],0,{data_dir});
    data_all = cellfun(@(x) load(x,'stats_refit'),file_list,'UniformOutput',0);
    data_all = cellfun(@(x) x.stats_refit,data_all,'UniformOutput',0);
    data = cellfun(@(x) x.pev, data_all,'UniformOutput',0);
    pev(i) = nanmean(cat(1,data{:})); 
end

data_mouse = cat(1,data_mouse{:});
%zscore
data_mouse = zscore(data_mouse,[],1);

%Build neurotype axis
[nt, Y_null, betas, betas_null, allstats] = DistanceToHyperplane(data_mouse,grp'+1,5,0,0.4); 
[betas, idx] = sort(betas,'ascend');
AUC = cellfun(@(x) x.AUC, allstats(1:end-1),'UniformOutput',0);
AUC = cat(2,AUC{:});
% nanmean(AUC)





