%Supplemental Figure 3
motif_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
gp = general_params_vpa;
file_list = GrabFiles('\w*chunk\w*.mat', 0, {motif_dir});
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);

%load all the motifs
temp = cellfun(@(x) load(x,'W','nanpxs'),file_list,'UniformOutput',0);
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

%assign group to each motif
allgroup = isVPA(motif_id);
for i = 1:size(W_basis,2)
    if sum(ismember(noise_clusters,i)>0)
    else
    %load all motifs from the 
    temp = W_alligned(:,cluster_idx==i,:);
    tempid = allgroup(cluster_idx==i);
    minval = min(sum(tempid==0),sum(tempid==1));
    sal = temp(:,tempid==0,:);
    N = sum(tempid==0);
    sal = nanmean(sal(:,randperm(N,minval),:),3);
    vpa = temp(:,tempid==1,:);
    N = sum(tempid==1);
    vpa = nanmean(vpa(:,randperm(N,minval),:),3);
    vpa_flat = [];
    for j = 1:minval
        temp = squeeze(vpa(:,j,:));
        vpa_flat(j,:)=temp(:);
    end
    sal_flat = [];
    for j = 1:minval
        temp = squeeze(sal(:,j,:));
        sal_flat(j,:)=temp(:);
    end
    classid=cat(1,zeros(minval,1),ones(minval,1));
    rng('default');
    data_temp = cat(1,sal_flat,vpa_flat);
    data_temp = zscore(data_temp,0,2);
    [features, Observed, Shuffled, ~] = SVMClassifier_Binary([data_temp,classid+1],[],'holdout',0.4,'optimize',0,'kernel','linear','featureselect','anova','nshuf',1,'numfeatures',500);    

    figure; hold on
    temp = zeros(1,size(vpa_flat,2));
    temp(features)=1;
    temp = reshape(temp,size(vpa,1),1,size(vpa,3));
    temp = sum(temp,3);
    feature_map = zeros(numel(nanpxs)+numel(temp),1);
    feature_map(~ismember(1:numel(feature_map),nanpxs))=temp;
    feature_map = reshape(feature_map,68,68);
    imagesc(feature_map);
    title(sprintf('motif %d %0.2g',i,Observed.AUC));
    end
end

%% Let's just visually compare all of them
%assign group to each motif
allgroup = isVPA(motif_id);
for i = 1:size(W_basis,2)
    if sum(ismember(noise_clusters,i)>0)
    else
    %load all motifs from the 
    temp = W_alligned(:,cluster_idx==i,:);
    tempid = allgroup(cluster_idx==i);
    minval = min(sum(tempid==0),sum(tempid==1));
    sal = temp(:,tempid==0,:);
    N = sum(tempid==0);
    sal = (sal(:,randperm(N,minval),:));
    vpa = temp(:,tempid==1,:);
    N = sum(tempid==1);
    vpa = (vpa(:,randperm(N,minval),:));
    vpa_flat = [];
    for j = 1:minval
        temp = squeeze(vpa(:,j,:));
        vpa_flat(j,:)=temp(:);
    end
    vpa_flat = zscore(vpa_flat,0,2);
    sal_flat = [];
    for j = 1:minval
        temp = squeeze(sal(:,j,:));
        sal_flat(j,:)=temp(:);
    end
    figure; 
    sal_flat = zscore(vpa_flat,0,2);
    classid=cat(1,zeros(minval,1),ones(minval,1));
    temp = nanmean(vpa_flat,1)-nanmean(sal_flat,1);    
    temp = nanmean(sal_flat,1);
    temp = nanmean(squeeze(reshape(temp,size(vpa,1),1,size(vpa,3))),2);
    feature_map = zeros(numel(nanpxs)+numel(temp),1);
    feature_map(~ismember(1:numel(feature_map),nanpxs))=temp;
    feature_map = reshape(feature_map,68,68);
    imagesc(feature_map);
    colorbar
    title(sprintf('motif %d',i));
    end
end


%%
%Load 
rng('default');
[~, Observed, Shuffled, ~] = SVMClassifier_Binary([data,grp],[],'holdout',0.20,'optimize',0,'kernel','rbf','featureselect','none','nshuf',1000);

%Plot the ROC curve
fp = fig_params_vpa;
figure; hold on; 
plot(Observed(:).X(:,1),Observed(:).Y(:,1),'linewidth',2,'color','k')
arrayfun(@(n) plot(Shuffled(n).X(:,1),Shuffled(n).Y(:,1),'linewidth',1,'color',[0.75 0.75 0.75 0.5]), 1:numel(Shuffled),'UniformOutput',0);
ylabel('True Positive Rate');
xlabel('False Positive Rate');
fp.FormatAxes(gca);

%Plot the AUC comparisons
figure; hold on; 
histogram([Shuffled(:).AUC],'BinWidth',0.05,'FaceColor',[0.75 0.75 0.75],'EdgeColor','k')
yval = get(gca,'ylim');
line([Observed.AUC,Observed.AUC],yval,'linewidth',2,'color','r');
p = sum([Shuffled(:).AUC,Observed.AUC]>=Observed.AUC)/numel([Shuffled(:).AUC,Observed.AUC]);
fp.SetTitle(gca,sprintf('p = %0.2g',p));
xlabel('AUC');
ylabel('Shuffles');
fp.FormatAxes(gca);
