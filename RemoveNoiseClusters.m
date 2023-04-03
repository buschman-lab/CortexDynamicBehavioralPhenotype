function [noise_idx,noise_clusters,criteria,levels] = RemoveNoiseClusters(cluster_idx,motif_id,tcorr_mat)

%must meet allof these criteria
levels = [10,0,0,0]; %number of mice, variance between clusters,average correlation, number of motifs

%Number of mice contributing to each motif
clust = unique(cluster_idx);
criteria = NaN(3,numel(clust));
for i = 1:numel(clust)
    criteria(1,i) = numel(unique(motif_id(cluster_idx==clust(i))));    
    temp = tcorr_mat(cluster_idx==clust(i),cluster_idx==clust(i));    
    temp(1:1+size(temp,1):end) = NaN;%remove the autocorrelation
    criteria(2,i) = nanvar(temp(:)); %variance in correlation
    criteria(3,i) = nanmean(temp(:)); %mean correlation
    criteria(4,i) = sum(cluster_idx==clust(i)); %number of contributing motifs
end

temp = arrayfun(@(n) criteria(n,:)>levels(n),1:numel(levels),'UniformOutput',0);
noise_clusters = find(any(cat(1,temp{:})==0,1)==1);

noise_idx = cell(1,numel(noise_clusters));
for i = 1:numel(noise_clusters)
   noise_idx{i} = cluster_idx==noise_clusters(i);
end

end



