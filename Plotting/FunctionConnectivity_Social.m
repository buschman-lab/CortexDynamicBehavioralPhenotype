%Get the phenotype and neurotype axes
%Show that either nt_social is correatled with behavior or that the change
%in nt correlated with behavior, with those on the outskirts haven't the
%greatest change (normalize both to zero and 1 and subtract
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
[pt,~,~,~,~,~,~,betas_pt,~,~,bias_pt] = Phenotype([data_dir,'\ALL'],'verbose',0,'num_classifiers',5,'holdout',0.5);
[nt,~, betas_nt, ~,allstats] = Neurotype([data_dir,'\ALL']);
fn = GrabFiles('chunk',0,{ 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL'});
mice_solo = unique(MouseNumFromPath(fn,'Mouse-'));
%%
rng('default')
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social';
[nt_social,~, betas_nt_social, ~,allstats_social,mice_social] = Neurotype_Social([data_dir,'\ALL'],'num_classifiers',4,'holdout',0.2,'avg_method','median');
nt_solo = nt(ismember(mice_solo,mice_social));
pt_solo = pt(ismember(mice_solo,mice_social));
grp = isVPA(mice_social);
% delta_nt = (normalize(nt,'range',[0 1])-normalize(nt_social,'range',[0 1]));
close all ;plot(fitlm(nt_solo,nt_social));
figure; plot(fitlm(pt_solo,nt_social));
rho = corr(pt_solo',nt_social')
%shuffle pt within group
rho_perm = [];
for i = 1:1000
   idx = find(grp==1); 
   pt_temp = pt_solo; 
   pt_temp(idx) = pt_temp(idx(randperm(numel(idx))));
   idx = find(grp==0);
   pt_temp(idx) = pt_temp(idx(randperm(numel(idx))));
   rho_perm(i) = corr(pt_temp',nt_social');
end
pval = sum(abs([rho_perm,rho])>abs(rho))/1001
%%
%     0.1000    2.0000
%     0.1000    3.0000
%     0.1000    4.0000
%     0.2000    4.0000
%     0.3000    4.0000
%     0.1000    5.0000
%     0.2000    5.0000

%% there is a reversal in motifs
[betas_nt_sorted, idx] = sort(betas_nt,'ascend');
betas_social_sorted= betas_nt_social(idx);
figure('position',[680   480   560   498]); hold on; 
b = bar(betas_social_sorted,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
line([0,0],[-1,numel(betas_social_sorted)+1],'color','k','linewidth',1.5)
set(gca,'xtick',[1,4,8,12,16],'XTickLabelRotation',0)
set(gca,'units','centimeters','position',[3,2,5,2.5])
set(gca,'xlim',[0.5,numel(betas_social_sorted)+0.5]);
set(gca,'ylim',[round(min(betas_social_sorted),1)-0.1,round(max(betas_social_sorted),1)+0.1],'ytick',[round(min(betas_social_sorted),1)-0.1,round(max(betas_social_sorted),1)+0.1])
ylabel('\bf \beta')


%% get the change in betas and show the the largest change is in the ones with larger local-global
delta_beta = normalize(betas_social_sorted,'range',[-1 1])-normalize(betas_nt_sorted,'range',[-1 1]);
figure('position',[680   480   560   498]); hold on; 
b = bar(delta_beta,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
line([0,0],[-1,numel(delta_beta)+1],'color','k','linewidth',1.5)
set(gca,'xtick',[1,4,8,12,16],'XTickLabelRotation',0)
set(gca,'units','centimeters','position',[3,2,5,2.5])
set(gca,'xlim',[0.5,numel(delta_beta)+0.5]);
set(gca,'ylim',[round(min(delta_beta),1)-0.1,round(max(delta_beta),1)+0.1],'ytick',[round(min(delta_beta),1)-0.1,round(max(delta_beta),1)+0.1])
ylabel('\bf \beta')


% get the long vs short per motif
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL\basismotifs.mat','W_basis','noise_clusters');
W_basis(:,noise_clusters,:)=[];
[betas_sorted, idx] = sort(betas_nt,'ascend');
W_basis=W_basis(:,idx,:); 
avg_fc = [];
for i =1:n_mot    
    w=squeeze(W_basis(:,i,:));
    %to avoid spurious inflation of distance in motifs with little activity and then a lot of activity 
    %we only want frames with activity above a threshold. Also, we don't
    %want to have this threshold to be depending on the size of the burst.
    %So take just the top 10 pixels and then keep the 
    idx_top = nanmean(maxk(w,20));
    idx = [find(idx_top>(0.05*max(idx_top)),1,'first'),find(idx_top>(0.05*max(idx_top)),1,'last')];
    %we care about the relative change within the motif... so normalize
    %each frame to 0-->1    
    w = w(:,idx(1):idx(2));
    [x,z] = size(w);
    x = sqrt(x);
    y = x; 
    roi_trace = zeros(z,size(roi,1));    
    for j = 1:size(roi,1) %roi loop 
        mask = zeros(x,y);
        mask(roi(j,2)-1:roi(j,2)+1,roi(j,1)-1:roi(j,1)+1)=1;
        mask = reshape(mask,[x*y,1]);
        roi_trace(:,j) = nanmean(w(mask==1,:));        
    end    
    avg_fc(:,:,i) = fisherZ(corr(roi_trace-nanmean(roi_trace)));  
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
d(d<4)=NaN;
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

motif_connectivity = long_conn./short_conn;

motif_group = betas_social_sorted>0; %positive betas = vpa
motif_connectivity(motif_group==1)
vpa = motif_connectivity(motif_group==1);
sal = motif_connectivity(motif_group==0);

figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:n_mot
    if motif_group(i)==1
        x = (1.5)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
    else
        x = (1)+rand(1)/4-0.125;
        plot(x,motif_connectivity(i),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
    end
end

perm_stat = NaN(1,1000);
rng('default');
for i = 1:1000
   temp_group = motif_group(randperm(numel(motif_group)));
   perm_stat(i) = nanmean(motif_connectivity(temp_group==1))-nanmean(motif_connectivity(temp_group==0));
end
true_val = nanmean(vpa)-nanmean(sal);
pval_right = sum(perm_stat>true_val)/1000
pval_left = sum(perm_stat<true_val)/1000




%%

fn = GrabFiles('refit',0,{data_dir});
temp = load(fn{101},'data','nanpxs');
temp = nanmax(conditionDffMat(temp.data(:,:,1)',temp.nanpxs),[],3);
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
plot(roi(:,1),roi(:,2),'.','color','r')

%% Grab the FC matrix for each epoch for social
fc_mat = cell(1,numel(fn));
mouse = cell(1,numel(fn));
for i = 1:numel(fn) %epoch loop
    i
    temp = load(fn{i},'data','nanpxs','mouseID');
    for q = 1:size(temp.data,3) %loop through chunks
        dff = conditionDffMat(temp.data(:,:,q)',temp.nanpxs);
        [x,y,z] = size(dff);
        dff = reshape(dff,[x*y,z]);
        roi_trace = zeros(z,size(roi,1));    
        for j = 1:size(roi,1) %roi loop 
            mask = zeros(x,y);
            mask(roi(j,2)-1:roi(j,2)+1,roi(j,1)-1:roi(j,1)+1)=1;
            mask = reshape(mask,[x*y,1]);
            roi_trace(:,j) = nanmean(dff(mask==1,:));        
        end
        fc_mat{i}(:,:,q) = corr(roi_trace-nanmean(roi_trace));    
        mouse{i}(q) = temp.mouseID;
    end
end
fc_mat = cat(3,fc_mat{:});
mouse = cat(2,mouse{:});


save('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data_social.mat','fc_mat','roi','z','mouse');
load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data_social.mat','fc_mat','roi','z','mouse');
mice = unique(mouse);
grp = isVPA(mice);
fp = fig_params_vpa;

fn_solo = GrabFiles('chunk',0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\'});
solo_mouse = MouseNumFromPath(fn_solo,'Mouse-');

%% Average FC matrix across epochs per animal
solo = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data.mat','fc_mat','roi','z');
mice = unique(mouse);
avg_fc = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
avg_fc_solo = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
for i = 1:numel(mice)
   idx = find(mouse == mice(i));   
   avg_fc(:,:,i) = nanmedian(fisherZ(fc_mat(:,:,idx)),3);
   idx = find(solo_mouse==mice(i));
   avg_fc_solo(:,:,i) = nanmedian(fisherZ(solo.fc_mat(:,:,idx)),3);
end

%remove autocorr and cross hemi
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
avg_fc(repmat(mask,1,1,numel(mice))==1)=NaN;
avg_fc_solo(repmat(mask,1,1,numel(mice))==1)=NaN;

%get the distance between each roi
d_mat = pdist2(roi,roi);
d_mat(mask==1)=NaN;

%plot the correlation by distance for each group
dist_bins = 1:10:51;
d_rho_social = cell(numel(mice),1);
d_rho_solo = cell(numel(mice),1);
for i = 1:numel(mice)
    temp = avg_fc(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    uniqd = unique(d);
    for j = 1:numel(uniqd)
        d_rho_social{i}(j) = nanmean(temp(d==uniqd(j)));
    end  
    temp = avg_fc_solo(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    uniqd = unique(d);
    for j = 1:numel(uniqd)
        d_rho_solo{i}(j) = nanmean(temp(d==uniqd(j)));
    end       
end
d_rho_social = cat(1,d_rho_social{:});
d_rho_solo = cat(1,d_rho_solo{:});

%Since different spatiotemporal variances between environments, normalize
%to the average across distances per animal in both environments
d_rho_social = d_rho_social-nanmean(d_rho_social,2);
d_rho_solo = d_rho_solo-nanmean(d_rho_solo,2);

d_rho = abs(d_rho_social)-abs(d_rho_solo);
% bar plot version
vpa = d_rho(grp==1,:);
sal = d_rho(grp==0,:);
figure; hold on; 
%do a bar plot for each distance with connected scatters
bar(1:5,nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
bar(1.5:5.5,nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.4)
for i = 1:numel(mice)
    if grp(i)==1
        x = (1.5:5.5)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_vpa,0.25])
    else
        x = (1:5)+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
%         plot(x,d_rho(i,:),'linewidth',0.75,'marker','none','markersize',10,'color',[fp.c_sal,0.25])
    end
end
% ylim([0.5 1.75])
%add stats
for i = 1:5
    true_val = nanmean(sal(:,i))-nanmean(vpa(:,i)); 
    val_perm=[];
    for j = 1:1000
       grp_shuf = grp(randperm(20));       
       val_perm(j) = nanmean(d_rho(grp_shuf==0,i))-nanmean(d_rho(grp_shuf==1,i));
    end
    pval = sum(abs([val_perm,true_val])>abs(true_val))/1001;        
    yval = nanmax(cat(1,vpa(:,i),sal(:,i)))+0.05;
    line([i,i+0.5],[yval,yval],'linewidth',1.5,'color',[0.25 0.25 0.25]);
%     text(i+0.25,yval,sprintf('p=%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom')
    text(i+0.25,yval,sprintf('%0.2g',pval),'Rotation',45,'HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',12)
end
% set(gca,'xtick',1.25:5.25,'xticklabels',dist_bins(2:end),'ytick',[0.5:0.25:1.75])
xlim([0.70 5.8])
xlabel('Distance (um)')
ylabel('Connectivity (\rho_z)','Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.5 6 6],[])
fp.FormatAxes(gca);
%%
%conclusion - the long range disconnectivity stablizes in social environment,
%this is because either the VPA mice increase their long/short ratio or
%Saline mice decrease their long/short ratio (test this)
%this is paired with a change in the sampling of motifs
%we see that the change in long/short ratio across animals is correlated
%with the change in long/short ratio of motifs

%to do: 

%% figure showing long/short ratio of all mice 
solo = load('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FunctionalConnectivity\fc_data.mat','fc_mat','roi','z');
mice = unique(mouse);
avg_fc = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
avg_fc_solo = zeros(size(fc_mat,1),size(fc_mat,2),numel(mice));
for i = 1:numel(mice)
   idx = find(mouse == mice(i));   
   avg_fc(:,:,i) = nanmedian(fisherZ(fc_mat(:,:,idx)),3);
   idx = find(solo_mouse==mice(i));
   avg_fc_solo(:,:,i) = nanmedian(fisherZ(solo.fc_mat(:,:,idx)),3);
end

%remove autocorr and cross hemi
mask = triu(ones(size(roi,1),size(roi,1)),0);  
%remove cross - hemi corrs
hemiidx = roi(:,1)<35;
hemiidx=(hemiidx+1).*(1+hemiidx'); %2=cross hemi, 1 or 4 = within hemi
mask = mask+tril(hemiidx==2,-1);
avg_fc(repmat(mask,1,1,numel(mice))==1)=NaN;
avg_fc_solo(repmat(mask,1,1,numel(mice))==1)=NaN;

%get the distance between each roi
d_mat = pdist2(roi,roi);
d_mat(mask==1)=NaN;

%plot the correlation by distance for each group
dist_bins = 1:10:51;
d_rho_social = NaN(1,numel(mice));
d_rho_solo = NaN(1,numel(mice));
for i = 1:numel(mice)
    temp = avg_fc(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    d(d<4)=NaN;
    d(~isnan(d))=1;
    d_rho_long = nanmean(temp(d==1));
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    d(d>2)=NaN;
    d(~isnan(d))=1;
    d_rho_short = nanmean(temp(d==1));
    d_rho_social(i) = d_rho_long/d_rho_short;    
    
    temp = avg_fc_solo(:,:,i);
    good_idx = ~isnan(temp(:));
    temp = temp(good_idx);
    %discretize the distances and average within each
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    d(d<4)=NaN;
    d(~isnan(d))=1;
    d_rho_long = nanmean(temp(d==1));
    d = d_mat(good_idx);
    d = discretize(d,dist_bins);
    d(d>2)=NaN;
    d(~isnan(d))=1;
    d_rho_short = nanmean(temp(d==1));
    d_rho_solo(i) = d_rho_long/d_rho_short;    
end
%compare the short/long ratio across groups between environments
vpa = [d_rho_solo(grp==1)',d_rho_social(grp==1)'];
sal = [d_rho_solo(grp==0)',d_rho_social(grp==0)'];
d_rho = [d_rho_solo',d_rho_social'];
figure; hold on; 
%do a bar plot for each distance with connected scatters
bar([1,1.5],nanmean(sal),'facecolor',fp.c_sal,'facealpha',0.5,'EdgeColor','none','barwidth',0.5)
bar([2.5,3],nanmean(vpa),'facecolor',fp.c_vpa,'facealpha',0.5,'EdgeColor','none','barwidth',0.5)

for i = 1:numel(mice)
    if grp(i)==1
        x = [2.5,3]+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_vpa,0.5])
        plot(x,d_rho(i,:),'linewidth',1.5,'marker','none','markersize',10,'color',[fp.c_vpa,0.25])
    else
        x = [1,1.5]+rand(1)/4-0.125;
        plot(x,d_rho(i,:),'linestyle','none','marker','.','markersize',10,'color',[fp.c_sal,0.5])
        plot(x,d_rho(i,:),'linewidth',1.5,'marker','none','markersize',10,'color',[fp.c_sal,0.25])
    end
end
ylim([0.5 .75])
%add stats
pval_vpa = signrank(vpa(:,1)-vpa(:,2));
pval_sal = signrank(sal(:,1)-sal(:,2));
yval = nanmax(cat(1,vpa(:),sal(:)))+0.02;
line([1,1.5],[yval,yval],'linewidth',1.5,'color',[0.25 0.25 0.25]);
text(1.25,yval,sprintf('%0.2g',pval_sal),'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',12)
line([2.5,3],[yval,yval],'linewidth',1.5,'color',[0.25 0.25 0.25]);
text(2.75,yval,sprintf('%0.2g',pval_vpa),'Rotation',0,'HorizontalAlignment','center','VerticalAlignment','bottom','fontsize',12)
set(gca,'xtick',[1,1.5,2.5,3],'xticklabels',{'Solo','Paired','Solo','Paired'},'ytick',[0.5:0.125:1.75],'XTickLabelRotation',45)
xlim([0.70 3.3])
ylabel('Long vs Short Connectivity','Interpreter','tex')
fp.FigureSizing(gcf,[2.5 2.5 6 6],[])
fp.FormatAxes(gca);
%% Correlate Motifs with largest change in beta weights to local-global of that motif














