function VPA_ClusterMotifs(motif_dir,save_dir,cur_group)

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis'));
end

fprintf('\nworking on group %d',cur_group)
if cur_group==2
    cur_group=[0,1];
end

gp = general_params_vpa;
file_list = GrabFiles('\w*chunk\w*.mat', 0, {motif_dir});
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);

%load all the motifs
temp = cellfun(@(x) load(x,'W','nanpxs'),file_list(ismember(group,cur_group)),'UniformOutput',0);
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

%initial clustering
fprintf('\n\t Generating Basis Motifs');
[W_basis, kval, ovr_q, cluster_idx, idx_knn, tcorr_mat, handles, lag_mat, lags, nanpxs] = ClusterW(W,gp,nanpxs);

%remove noise clusters and their representations in the original data
[noise_idx,noise_clusters,criteria,levels] = RemoveNoiseClusters(cluster_idx,motif_id,tcorr_mat);

save([save_dir filesep 'basismotifs.mat'],'W_basis', 'file_list','kval', 'ovr_q', 'cluster_idx', 'idx_knn', 'tcorr_mat','gp','lags','lag_mat','nanpxs',...
    'noise_idx','noise_clusters','criteria','levels','-v7.3')

saveCurFigs(handles,'-dpng','ClusteringMotifs',save_dir,0); close all;

fprintf('\n\t Done Saving Basis Motifs - Saved');

end













