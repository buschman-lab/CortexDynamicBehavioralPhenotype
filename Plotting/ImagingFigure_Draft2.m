%ImagingFigure();
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
fp = fig_params_vpa; 

% Load refitting data
file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
grp = isVPA(unique(mouse));
data_all = cellfun(@(x) load(x,'stats_refit'),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.stats_refit,data_all,'UniformOutput',0);
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%zscore
data_mouse = zscore(data_mouse,[],1);

%Build neurotype axis
[nt, Y_null, betas, betas_null, allstats] = DistanceToHyperplane(data_mouse,grp'+1,5,0,0.4); 
[betas, idx] = sort(betas,'ascend');
AUC = cellfun(@(x) x.AUC, allstats(1:end-1),'UniformOutput',0);
AUC = cat(2,AUC{:});
% nanmean(AUC)

%% Order Motifs by their beta weights
%load the basis motifs
temp = load([data_dir filesep 'basismotifs.mat'],'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W = temp.W_basis;
W(:,temp.noise_clusters,:)=[];
W_noise = temp.W_basis(:,temp.noise_clusters,:);

[P,M,T] = size(W);
W = W(:,idx,:);

%% Plot each motif 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\Imaging99_lesssmooth_temp';
if ~exist(savedir)
    mkdir(savedir)
end
PlotMotifs(W(:,16,:),'SaveDir',savedir,'caxis',[0 98.5],'kernel',[1 1],'includeflow',0);   

%% Plot the average image of motifs positively/negatively 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\AverageBBA';
if ~exist(savedir)
    mkdir(savedir)
end
PlotMotifs(W(:,16,:),'SaveDir',savedir,'caxis',[0 98.5],'kernel',[1 1],'includeflow',0);   
%normalize between zero and 1 for each motif
W_norm = W;
for i = 1:size(W,2)
    temp = squeeze(W(:,i,:));
    temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
    W_norm(:,i,:) = temp; 
end

MP = nanmean(W_norm(:,[7,8,9,14,15,16],:),[2,3]);
PlotMotifs(MP,'SaveDir',savedir,'caxis',[0 97.5],'kernel',[1 1],'includeflow',0,'NameStr','MP');  
MN=nanmean(W_norm(:,[2,3,4,5,6],:),[2,3]);
PlotMotifs(MN,'SaveDir',savedir,'caxis',[0 97.5],'kernel',[1 1],'includeflow',0,'NameStr','MN');  
SP=nanmean(W_norm(:,[2,5,6,7,9,10],:),[2,3]);
PlotMotifs(SP,'SaveDir',savedir,'caxis',[0 97.5],'kernel',[1 1],'includeflow',0,'NameStr','SP');  
SN=nanmean(W_norm(:,[1,3,4,15,10,16],:),[2,3]);
PlotMotifs(SN,'SaveDir',savedir,'caxis',[0 97.5],'kernel',[1 1],'includeflow',0,'NameStr','SN');  



%% Plot the noise motifs

savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\Imaging99.5_smooth_noise';
if ~exist(savedir)
    mkdir(savedir)
end
PlotMotifs(W_noise,'SaveDir',savedir,'caxis',[0 99.5],'kernel',[1 1],'includeflow',0);  

%%
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\Neurotype';


%% Compare PEV Between Groups

data = cellfun(@(x) x.pev, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:})*100;

close all
rng('default');
sal_pos = 1+rand(sum(grp==0),1)/2-0.25;
vpa_pos = 2+rand(sum(grp==1),1)/2-0.25;
figure('position',[680   558   240   271]); hold on;      
vpa = data_mouse(grp==1);
sal = data_mouse(grp==0);
b = bar([nanmean(sal),nanmean(vpa)],'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);    
b.CData(1,:) = [fp.c_sal];
b.CData(2,:) = [fp.c_vpa];
plot(sal_pos,sal,'.','markersize',fp.m_markersize,'color',fp.c_sal)
plot(vpa_pos,vpa,'.','markersize',fp.m_markersize,'color',fp.c_vpa)      
errorbar(1,nanmean(sal),sem(sal),'LineWidth',2,'Color','k');
errorbar(2,nanmean(vpa),sem(vpa),'LineWidth',2,'Color','k');
p = ranksum(sal,vpa);
temp = max(cat(1,sal,vpa));
line([1,2],[temp+temp*0.025,temp+temp*0.025],'linewidth',1.5,'color','k');
text(1.5, temp+temp*0.05,sprintf('%0.2g',p),'Rotation',0,'FontSize',fp.sig_fontsize,'HorizontalAlignment','center');
set(gca,'Clipping','off','box','off');   
set(gca,'units','centimeters','position',[3 2 2.75 3.25])
set(gca,'xtick',[1,2],'xticklabel',{'SAL','VPA'},'xticklabelrotation',90,'xlim',[0.25 2.75])
ylabel({'Percent Explained','Variance'})
fp.FormatAxes(gca);
set(gca,'LineWidth',2)

%stats
nanmean(data_mouse(grp==1));
nanmean(data_mouse(grp==0));
ci = bootci(100,@nanmean,data_mouse(grp==1));
ci = bootci(100,@nanmean,data_mouse(grp==0));

%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','GroupedPEV',savedir,1); close all

%% plot the ROCs
% X = cellfun(@(x) x.X(:,1), allstats,'UniformOutput',0);
% X = cat(2,X{:});
% Y = cellfun(@(x) x.Y(:,1), allstats,'UniformOutput',0);
% Y = cat(2,Y{:});
AUC = cellfun(@(x) x.AUC, allstats,'UniformOutput',0);
AUC = cat(2,AUC{:});
% 
% figure; hold on; 
% plot(X,Y,'color',[0.5 0.5 0.5,0.3],'linewidth',1.5)
% plot(nanmean(X,2),nanmean(Y,2),'color','k','linewidth',2)
% xlabel('False Positive Rate')
% ylabel('True Positive Rate')
% title(sprintf('AUC = %0.2g',nanmean(AUC)),'FontWeight','normal');
% set(gca,'units','centimeters','position',[2 2 4 4]);
% fp.FormatAxes(gca)
% handles = get(groot, 'Children');
% fp.SaveFigs(handles,'-svg','neuralAUC',savedir,1); close all

%report AUC 
AUC = [];
for i = 1:numel(allstats)-1
   AUC(i) = allstats{i}.AUC
end

%% plot the betas
% figure('position',[680   480   560   498]); hold on; 
% b = barh(betas,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
% line([0,0],[-1,numel(betas)+1],'color','k','linewidth',1.5)
% set(gca,'ytick',[1,4,8,12,16],'YTickLabelRotation',0,'ydir','reverse','TickLength',[.05,.2])
% set(gca,'units','centimeters','position',[3,2,4,4.5],'YAxisLocation','left')
% set(gca,'ylim',[0.5,numel(betas)+0.5]);
% set(gca,'xlim',[round(min(betas),1)-0.1,round(max(betas),1)+0.1],'xtick',[round(min(betas),1)-0.1,round(max(betas),1)+0.1])
% xlabel('\bf \beta')

figure('position',[680   480   560   498]); hold on; 
b = bar(betas,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
line([0,0],[-1,numel(betas)+1],'color','k','linewidth',1.5)
set(gca,'xtick',[1,4,8,12,16],'XTickLabelRotation',0)
% set(gca,'units','centimeters','position',[3,2,2.5,2])
set(gca,'units','centimeters','position',[3,2,5,2.5]) %large
set(gca,'xlim',[0.5,numel(betas)+0.5]);
set(gca,'ylim',[round(min(betas),1)-0.1,round(max(betas),1)+0.1],'ytick',[round(min(betas),1)-0.1,round(max(betas),1)+0.1])
ylabel('\bf \beta')

fp.FormatAxes(gca);

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','neuralbetas_small',savedir,1); close all

%% visualize the spectrum of motifs
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%zscore
data_mouse = zscore(data_mouse,[],1);
%compute in PC space
[coef, score, latent, ~, explained, mu] = pca(data_mouse);

%Fit the hyperplane with classifier
partition = [];
partition.training = true(numel(grp),1); %train on entire dataset
partition.test = true(numel(grp),1); %doesn't matter as you aren't testing this classifier
rng('default'); %for reproducibility
[~, observed, ~, ~] = SVMClassifier_Binary([score(:,1:3),grp'+1],partition,'optimize',0,'kernel','linear','featureselect','none','nshuf',0,'verbose',0); 
mdl = observed.Classifier;
X = score(:,1:3);
%Gather support vectors from ClassificationSVM struct
sv =  mdl.SupportVectors;
%set step size for finer sampling
d =0.025;
%generate grid for predictions at finer sample rate
[x, y, z] = meshgrid(min(X(:,1)):d:max(X(:,1)),...
min(X(:,2)):d:max(X(:,2)), min(X(:,3)):d:max(X(:,3)));
xGrid = [x(:),y(:),z(:)];
%get scores, f
[ ~ , f] = predict(mdl,xGrid);
%reshape to same grid size as the input
f = reshape(f(:,2), size(x));
% Assume class labels are 1 and 0 and convert to logical
t = logical(grp');
%plot data points, color by class label
figure; hold on; 
plot3(score(grp==0,1),score(grp==0,2),score(grp==0,3),'.','markersize',30,'color',fp.c_sal)
plot3(score(grp==1,1),score(grp==1,2),score(grp==1,3),'.','markersize',30,'color',fp.c_vpa)
% load unscaled support vectors for plotting
% plot3(sv(:, 1), sv(:, 2), sv(:, 3), 'ro','markersize',9,'linewidth',1.5);
%plot decision surface
[faces,verts,~] = isosurface(x, y, z, f, 0, x);
patch('Vertices', verts, 'Faces', faces, 'FaceColor','r','edgecolor', 'none', 'FaceAlpha', 0.15);
grid on
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); 
set(gca,'units','centimeters','position',[3,2,6,6])
fp.FormatAxes(gca);
box on
campos([ 60.6485  -41.8919   16.3173])
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','pcaneurospectrum',savedir,1); close all

%% Compare the weightings of motifs across groups
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:})*100;
fp = fig_params_vpa; 
close all

rng('default');
sal_pos = rand(sum(grp==0),1)/4-0.125;
vpa_pos = rand(sum(grp==1),1)/4-0.125;
%plot them all on the same verticle axis
figure; hold on; 
for i = 1:size(data_mouse,2)%start at the second column  
    vpa = data_mouse(grp==1,i);
    sal = data_mouse(grp==0,i);
    b = bar([i,i+0.4],[nanmean(sal),nanmean(vpa)],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor',[1 1 1]);    
    b.CData(1,:) = [fp.c_sal_gray];
    b.CData(2,:) = [fp.c_vpa_gray];
    errorbar(i,nanmean(sal),sem(sal),'LineWidth',1.25,'Color',[0.5 0.5 0.5]);
    errorbar(i+0.4,nanmean(vpa),sem(vpa),'LineWidth',1.25,'Color',[0.5 0.5 0.5]);    
    plot(sal_pos+i,sal,'.','markersize',fp.m_markersize/1.75,'color',fp.c_sal_gray)
    plot(vpa_pos+i+0.4,vpa,'.','markersize',fp.m_markersize/1.75,'color',fp.c_sal_gray) 
end
set(gca,'units','centimeters','position',[2 2 6.75 2.75])
fp.FormatAxes(gca);
set(gca,'LineWidth',2)
handles = get(groot, 'Children');
% fp.SaveFigs(handles,'-svg','allmotifstests',savedir,1); close all
%% Version with no bars
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:})*100;
fp = fig_params_vpa; 
close all

rng('default');
sal_pos = rand(sum(grp==0),1)/4-0.125;
vpa_pos = rand(sum(grp==1),1)/4-0.125;
%plot them all on the same verticle axis
figure; hold on; 
for i = 1:size(data_mouse,2)%start at the second column  
    vpa = data_mouse(grp==1,i);
    sal = data_mouse(grp==0,i);
    plot(sal_pos+i,sal,'.','markersize',fp.m_markersize/1.75,'color',fp.c_sal)
    plot(vpa_pos+i+0.4,vpa,'.','markersize',fp.m_markersize/1.75,'color',fp.c_vpa)     
    errorbar(i,nanmean(sal),sem(sal),'LineWidth',1.25,'Color','k');
    errorbar(i+0.4,nanmean(vpa),sem(vpa),'LineWidth',1.25,'Color','k');    
end
set(gca,'units','centimeters','position',[2 2 6.75 2.75])
set(gca,'YAxisLocation','right')
xlim([0.25 17])
fp.FormatAxes(gca);
set(gca,'LineWidth',2)
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','allmotifstests_nobar',savedir,1); close all
%%
%Results of Individual Behavioral Tests
% rng('default');
% sal_pos = 1+rand(sum(grp==0),1)/2-0.25;
% vpa_pos = 2+rand(sum(grp==1),1)/2-0.25;
% for i = 1:size(data_mouse,2)%start at the second column       
%     figure('position',[680   558   240   271]); hold on;
%     vpa = data_mouse(grp==1,i);
%     sal = data_mouse(grp==0,i);
%     b = bar([nanmean(sal),nanmean(vpa)],'FaceColor','flat','FaceAlpha',0.5,'EdgeColor',[1 1 1]);    
%     b.CData(1,:) = [fp.c_sal_gray];
%     b.CData(2,:) = [fp.c_vpa_gray];
%     plot(sal_pos,sal,'.','markersize',fp.m_markersize/1.5,'color',fp.c_dot)
%     plot(vpa_pos,vpa,'.','markersize',fp.m_markersize/1.5,'color',fp.c_dot)      
%     errorbar(1,nanmean(sal),sem(sal),'LineWidth',1.5,'Color','k');
%     errorbar(2,nanmean(vpa),sem(vpa),'LineWidth',1.5,'Color','k');
%     [p,~] = ranksum(sal,vpa);
%     temp = max(cat(1,sal,vpa));
%     line([1,2],[temp+temp*0.025,temp+temp*0.025],'linewidth',1.5,'color','k');
%     text(1.5, temp+temp*0.05,sprintf('%0.2g',p),'Rotation',0,'FontSize',fp.sig_fontsize,'HorizontalAlignment','center');
%     set(gca,'Clipping','off','box','off');   
%     set(gca,'units','centimeters','position',[2 2 1.25 3])
%     set(gca,'xtick',[1,2],'xticklabel',{'SAL','VPA'},'xticklabelrotation',90,'xlim',[0.25 2.75])
%     ylabel(sprintf('Motif %d',i))
%     fp.FormatAxes(gca);
%     set(gca,'LineWidth',2)
% end    



%% Compare Phenotype and Neurotype
%Global Phenotype Axis
pt = Phenotype(data_dir,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Social Axis
pt_social = Phenotype(data_dir,'behavior_col',2:5,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Developmental Axis
pt_dev = Phenotype(data_dir,'behavior_col',6:8,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Motor Phenotype Axis
pt_motor = Phenotype(data_dir,'behavior_col',9:13,'verbose',0,'num_classifiers',5,'holdout',0.5);

%get the correlations
temp = cat(1,pt,pt_social,pt_motor,pt_dev);
rho = NaN(1,size(temp,2));
rho_sal = NaN(1,size(temp,2));
rho_vpa = NaN(1,size(temp,2));
for i = 1:size(temp,1)
    [rho(i),~] = corr(temp(i,:)',nt(:)); %whole correlation
    [rho_sal(i),~] = corr(temp(i,grp==0)',nt(grp==0)'); %sal correlation
    [rho_vpa(i),~] = corr(temp(i,grp==1)',nt(grp==1)'); %vpa correlation
end

%within group permutation of neurotype
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%zscore
data_mouse = zscore(data_mouse,[],1);

rng('default');
n_perm = 1000;
perm_data = NaN(size(data_mouse,1),size(data_mouse,2),n_perm);
for i = 1:n_perm
    temp_sal = data_mouse(grp==0,:);
    temp_sal = temp_sal(randperm(size(temp_sal,1),size(temp_sal,1)),:);
    temp_vpa = data_mouse(grp==1,:);    
    temp_vpa = temp_vpa(randperm(size(temp_vpa,1),size(temp_vpa,1)),:);
    perm_data(grp==0,:,i) = temp_sal;
    perm_data(grp==1,:,i) = temp_vpa;
end
rho_perm = NaN(4,n_perm);
rho_sal_perm = NaN(4,n_perm);
rho_vpa_perm = NaN(4,n_perm);
for i = 1:n_perm   
   if mod(i,10)==0; fprintf('\t\n Working on shuffle %d of %d',i,n_perm); end
   nt_perm = DistanceToHyperplane(perm_data(:,:,i),grp'+1,5,0,0.4);           
   rho_perm(1,i) = corr(pt(:),nt_perm'); %whole correlation
   rho_sal_perm(1,i) = corr(pt(grp==0)',nt_perm(grp==0)'); %sal correlation
   rho_vpa_perm(1,i)= corr(pt(grp==1)',nt_perm(grp==1)'); %vpa correlation

   rho_perm(2,i) = corr(pt_social(:),nt_perm'); %whole correlation
   rho_sal_perm(2,i) = corr(pt_social(grp==0)',nt_perm(grp==0)'); %sal correlation
   rho_vpa_perm(2,i)= corr(pt_social(grp==1)',nt_perm(grp==1)'); %vpa correlation
   
   rho_perm(3,i) = corr(pt_motor(:),nt_perm'); %whole correlation
   rho_sal_perm(3,i) = corr(pt_motor(grp==0)',nt_perm(grp==0)'); %sal correlation
   rho_vpa_perm(3,i)= corr(pt_motor(grp==1)',nt_perm(grp==1)'); %vpa correlation
   
   rho_perm(4,i) = corr(pt_dev(:),nt_perm'); %whole correlation
   rho_sal_perm(4,i) = corr(pt_dev(grp==0)',nt_perm(grp==0)'); %sal correlation
   rho_vpa_perm(4,i)= corr(pt_dev(grp==1)',nt_perm(grp==1)'); %vpa correlation   
end

%get the average values
rho_perm_avg = fisherInverse(nanmean(fisherZ(rho_perm(1,:))));
rho_perm_ci = fisherInverse(bootci(1000,@nanmean,fisherZ(rho_perm(1,:))));
rho_perm_diff = rho(1)-rho_perm_avg;

%p value of right tailed test
p_all = NaN(1,4);
p_sal = NaN(1,4);
p_vpa = NaN(1,4);
for i = 1:4
    p_all(i) = sum([rho(i),rho_perm(i,:)]>=rho(i))/(n_perm+1);
    p_sal(i) = sum([rho_sal(i),rho_sal_perm(i,:)]>=rho_sal(i))/(n_perm+1);
    p_vpa(i) = sum([rho_vpa(i),rho_vpa_perm(i,:)]>=rho_vpa(i))/(n_perm+1);
end

figure('position',[680   413   560   565]); hold on; 
lm = fitlm(pt,nt); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit   
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper 
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;

lm = fitlm(pt(grp==0),nt(grp==0)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = fp.c_sal; %fit
p(3).Color = fp.c_sal; %bounds lower
p(4).Color = fp.c_sal; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

lm = fitlm(pt(grp==1),nt(grp==1)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = fp.c_vpa; %fit
p(3).Color = fp.c_vpa; %bounds lower
p(4).Color = fp.c_vpa; %bound upper 
p(2).LineWidth = 1.5;
p(3).LineWidth = 1.;
p(4).LineWidth = 1.;

plot(pt(grp==0),nt(grp==0),'.','markersize',30,'color',fp.c_sal)
plot(pt(grp==1),nt(grp==1),'.','markersize',30,'color',fp.c_vpa)
legend off
% set(gca,'xlim',[round(min(pt),1)-0.1,round(max(pt),1)+0.1])
% set(gca,'ylim',[round(min(nt),1)-0.1,round(max(nt),1)+0.1])
set(gca,'xlim',[-3.4 3.4])
set(gca,'ylim',[-3.4 3.4])
set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
    'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
ylabel('Neurotype')
xlabel('Phenotype')
set(gca,'units','centimeters','position',[3.25 2 8.5 8.5]);
fp.SetTitle(gca,{sprintf('rho=%0.2g, rho_adj=%0.2g p=%0.2g rho=%0.2g \n p=%0.2g rho=%0.2g p=%0.2g',rho(1), rho(1)-nanmean(rho_perm(1,:)), p_all(1), rho_sal(1), p_sal(1), rho_vpa(1), p_vpa(1))});
fp.FormatAxes(gca)

% plot the shuffled distributions
label = {'all','social','motor','dev'};
for i = 1:numel(label)
   figure; hold on; 
   histogram(rho_perm(i,:),'BinWidth',0.01,'EdgeAlpha',0,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.25)
%    histogram(rho_perm(i,rho_perm(i,:)>=rho(i)),'BinWidth',0.01,'EdgeAlpha',0,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.5)
   yval = get(gca,'ylim');
   line([rho(i),rho(i)],yval,'color','r','linestyle','-','linewidth',2);
   line([nanmean(rho_perm(i,:)),nanmean(rho_perm(i,:))],yval,'color','b','linestyle',':','linewidth',1.5);
   set(gca,'ylim',yval);
   ylabel('Number Permutations')
   xlabel('Rho')
   set(gca,'units','centimeters','position',[3.25 2 4 2.25]);
   fp.FormatAxes(gca)
end

label = {'social','motor','dev'};
temp_all = cat(1,pt_social,pt_motor,pt_dev);
for i = 1:3
    temp = temp_all(i,:);
    figure; hold on; 
    lm = fitlm(temp,nt); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper
    
    
    lm = fitlm(temp(grp==0),nt(grp==0)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = fp.c_sal; %fit
    p(3).Color = fp.c_sal; %bounds lower
    p(4).Color = fp.c_sal; %bound upper 

    lm = fitlm(temp(grp==1),nt(grp==1)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = fp.c_vpa; %fit
    p(3).Color = fp.c_vpa; %bounds lower
    p(4).Color = fp.c_vpa; %bound upper 
    
    plot(temp(grp==0),nt(grp==0),'.','markersize',15,'color',fp.c_sal)
    plot(temp(grp==1),nt(grp==1),'.','markersize',15,'color',fp.c_vpa)
    legend off
    set(gca,'xlim',[round(min(temp),1)-0.1,round(max(temp),1)+0.1])
    set(gca,'ylim',[round(min(nt),1)-0.1,round(max(nt),1)+0.1])
    set(gca,'xtick',get(gca,'xlim'),'ytick',get(gca,'ylim'),...
        'xticklabel',{'SAL-like','VPA-like'},'yticklabel',{'SAL-like','VPA-like'});
    ylabel('Neurotype')
    xlabel(label{i})
%     [rho,~] = corr(temp(:),nt(:)); %correlation
    
    set(gca,'units','centimeters','position',[3 2 3.75 3.75]);
    fp.SetTitle(gca,{sprintf('rho=%0.2g, rho_adj=%0.2g, p=%0.2g',rho(i+1),rho(i+1)-nanmean(rho_perm(i+1,:)),p_all(i+1))});
    fp.FormatAxes(gca)
end


%% Randomly permute neurotypes within groups
%within group permutation of neurotype
nPerm =1000;
perm_stat = NaN(4,nPerm);
rng('default');
for i = 1:nPerm
   nt_perm = nt; 
   temp = nt(grp==0);
   nt_perm(grp==0)=temp(randperm(numel(temp),numel(temp)));
   temp = nt(grp==1);
   nt_perm(grp==1)=temp(randperm(numel(temp),numel(temp)));      
   perm_stat(1,i) = corr(pt(:),nt_perm(:));
   perm_stat(2,i) = corr(pt_social(:),nt_perm(:));
   perm_stat(3,i) = corr(pt_motor(:),nt_perm(:));
   perm_stat(4,i) = corr(pt_dev(:),nt_perm(:));
end

%get the correlations
temp = cat(1,pt,pt_social,pt_motor,pt_dev);
rho = NaN(1,size(temp,2));
rho_sal = NaN(1,size(temp,2));
rho_vpa = NaN(1,size(temp,2));
for i = 1:size(temp,1)
    [rho(i),~] = corr(temp(i,:)',nt(:)); %whole correlation
    [rho_sal(i),~] = corr(temp(i,grp==0)',nt(grp==0)'); %sal correlation
    [rho_vpa(i),~] = corr(temp(i,grp==1)',nt(grp==1)'); %vpa correlation
end

% plot the shuffled distributions
label = {'all','social','motor','dev'};
for i = 1:numel(label)
   figure; hold on; 
   histogram(perm_stat(i,:),'BinWidth',0.01,'EdgeAlpha',0,'FaceColor',[0.25 0.25 0.25],'FaceAlpha',0.25)
   yval = get(gca,'ylim');
   line([rho(i),rho(i)],yval,'color','r','linestyle','-','linewidth',2);
   line([nanmean(perm_stat(i,:)),nanmean(perm_stat(i,:))],yval,'color','b','linestyle',':','linewidth',1.5);
   set(gca,'ylim',yval);
   ylabel('Number Permutations')
   xlabel('Rho (permuted neuroaxis)')
   set(gca,'units','centimeters','position',[3.25 2 4 2.25]);
   fp.FormatAxes(gca)
end


% %get the average values
% rho_perm_within_avg = fisherInverse(nanmean(fisherZ(perm_stat(1,:))));
% rho_perm_within_ci = fisherInverse(bootci(1000,@nanmean,fisherZ(perm_stat(1,:))));
% rho_perm_within_diff = rho(1)-rho_perm_avg;

%p value of right tailed test
p_all_within = NaN(1,4);
for i = 1:4
    p_all_within(i) = sum(perm_stat(i,:)>=rho(i))/(n_perm+1);
end

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','neurotype_vs_phenotype_9_19_2021_Resizing',savedir,1); close all
%% Correlated each motif with the behavioral axis and plot
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%zscore
% data_mouse = zscore(data_mouse,[],1);
data_mouse_sorted = data_mouse(:,idx);

for i = 1:numel(idx)
    nt_temp = data_mouse_sorted(:,i)*100;
    figure; hold on; 
    lm = fitlm(pt,nt_temp); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper
    plot(pt(grp==0),nt_temp(grp==0),'.','markersize',30,'color',fp.c_sal)
    plot(pt(grp==1),nt_temp(grp==1),'.','markersize',30,'color',fp.c_vpa)
    legend off
    set(gca,'xlim',[round(min(pt),1)-0.1,round(max(pt),1)+0.1])
    set(gca,'ylim',[round(min(nt_temp),1)-0.1,round(max(nt_temp),1)+0.1])
    set(gca,'xtick',get(gca,'xlim'),...
        'xticklabel',{'SAL-like','VPA-like'});
    ylabel(sprintf('Motif %d',i))
    xlabel('Phenotype')
    [rho,~] = corr(pt(:),nt_temp(:)); %correlation    
    
    perm_stat = PermutationStatistic(cat(1,pt(:),nt_temp(:)),cat(1,zeros(size(pt(:))),ones(size(nt_temp(:)))), @(x,y) corr(x,y),1000);
    if rho<0
        p = sum(perm_stat<rho)/(1000+1);
    else
        p = 1- sum(perm_stat<rho)/(1000+1);
    end
    
    set(gca,'units','centimeters','position',[3 2 7 7]);
    fp.SetTitle(gca,{sprintf('and Neurotypes rho=%0.2g p=%0.2g',rho,p)});
    fp.FormatAxes(gca)
end

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','single_motifs_vspheno',savedir,1); close all


%% Correlation between each motif and each behavioral axis
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%zscore
data_mouse = zscore(data_mouse,[],1);
data_mouse_sorted = data_mouse(:,idx);

%Correlated each motif with each axis
rho_mat = [corr(data_mouse_sorted,pt'),corr(data_mouse_sorted,pt_social'),corr(data_mouse_sorted,pt_motor'),corr(data_mouse_sorted,pt_dev')]';

%make figure
figure('position',[680   739   560   239]); hold on; 
imagesc(rho_mat,[-0.3 0.3]); 
colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));

xlim([0.5,size(rho_mat,2)+0.5])
ylim([0.5,size(rho_mat,1)+0.5])
set(gca,'ytick',1:4,'yticklabel',{'\bfGlobal','Social','Motor','Dev'},'ydir','reverse')
set(gca,'xtick',[1,4,8,12,16])
ylabel('Phenotype');
xlabel('Motif');
set(gca,'units','centimeters','position',[3 2 8 4])
c = colorbar;
set(c,'ytick',[-0.3,0,0.3],'units','centimeters','position',[11.25    2    0.2    4])
fp.FormatAxes(gca)

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','motifsvspheno',savedir,1); close all
%% Plot the neurotype axis

sm = 20;
lg = 30;
rng('default')
figure('position',[680   558   669   420]); hold on;
x = rand(numel(nt),1)*0.5;
plot(nt(grp==0),x(grp==0)-0.5,'.','markersize',lg,'color',fp.c_sal);     
plot(nt(grp==1),x(grp==1)-0.5,'.','markersize',lg,'color',fp.c_vpa);       

line([0,0],[-1,4],'linewidth',2,'linestyle','--','color','k')
set(gca,'ylim',[-0.7,0.2],'ytick',[])
set(gca,'xtick',[min(get(gca,'xtick')),max(get(gca,'xtick'))],'xticklabel',{'SAL-like','VPA-like'})
set(gca,'units','centimeters','position',[3 2 6 4])
fp.FormatAxes(gca);

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','neuralaxis',savedir,1); close all























%% Compare the motif Loadigns Between Groups
% data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
% data = cat(1,data{:});    
% data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
% data_mouse = cat(1,data_mouse{:})*100;
% 
% close all
% figure('position',[ 290         558        1392         271]); hold on;     
% rng('default');
% sal_pos = rand(sum(grp==0),1)/2-0.25;
% vpa_pos = 1+rand(sum(grp==1),1)/2-0.25;
% pos = 1:3:size(data_mouse,2)*3;
% 
% for i = 1:size(data_mouse,2)     
%     vpa = data_mouse(grp==1,i);
%     sal = data_mouse(grp==0,i);
%     b = bar([pos(i),pos(i)+1],[nanmean(sal),nanmean(vpa)],'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);    
%     b.CData(1,:) = [fp.c_sal];
%     b.CData(2,:) = [fp.c_vpa];
%     plot(sal_pos+pos(i),sal,'.','markersize',fp.m_markersize,'color',fp.c_sal)
%     plot(vpa_pos+pos(i),vpa,'.','markersize',fp.m_markersize,'color',fp.c_vpa)      
%     errorbar(pos(i),nanmean(sal),sem(sal),'LineWidth',2,'Color','k');
%     errorbar(pos(i)+1,nanmean(vpa),sem(vpa),'LineWidth',2,'Color','k');
%     [p,~] = ranksum(sal,vpa);
%     temp = max(cat(1,sal,vpa));
%     line([pos(i),pos(i)+1],[temp+temp*0.025,temp+temp*0.025],'linewidth',1.5,'color','k');
%     text(pos(i)+0.5, temp+temp*0.125,sprintf('%0.2g',p),'Rotation',0,'FontSize',fp.sig_fontsize,'HorizontalAlignment','center');
% end
% 
% set(gca,'xtick',[],'xlim',[pos(1)-1, pos(end)+2],'ylim',[0 round(max(data_mouse(:)),-1)])
% ylabel({'Relative Percent','Explained Variance'})
% set(gca,'LineWidth',2)
% set(gca,'Clipping','off','box','off');   
% % set(gca,'units','centimeters','position',[3 2 3 4.5])
% fp.FormatAxes(gca);


% %% Compute motif complexity
% score = cell(1,M); 
% bad_pxl = cell(1,M);
% complexity = zeros(1,M);
% for i = 1:M   
%    warning('off'); temp = squeeze(W(:,i,:)); 
%    %remove non-active pixels (will inflate correlation)
%    bad_pxl{i} = nanvar(temp,[],2)<=eps;
%    temp(bad_pxl{i},:)=[];
%    [coef, score{i}, latent, ~, explained, mu] = pca(temp);
%    complexity(i) = find(cumsum(explained)>95,1,'first'); warning('on'); 
% end

%plot the scores for each; 
% for i = 1:M
%    for n = 1:complexity(i)
%        temp = NaN(numel(bad_pxl{i}),1);
%        temp(bad_pxl{i}==0) = score{i}(:,n);
%        imagesc(reshape(temp,[68 68]));
%        title(sprintf('pc %d',n)); axis off
%        colormap magma       
%        saveCurFigs(gcf,'-dpng',sprintf('Motif%dPC%d',i,n),pwd,0); close all;
%    end
% end
% data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
% file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
% % file_list = substring_filename(file_list,'nodynamics',opts.nodynamics);
% [imaged_mice, ~] = unique(MouseNumFromFileName(file_list));
% data_behavioral = GetBehavioralData();
% 
% % Load refitting data
% file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
% [mouse, ~] = MouseNumFromFileName(file_list);
% grp = isVPA(unique(mouse));
% data_all = cellfun(@(x) load(x,'stats_refit'),file_list,'UniformOutput',0);
% data_all = cellfun(@(x) x.stats_refit,data_all,'UniformOutput',0);
% data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
% data = cat(1,data{:});    
% data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
% data_mouse = cat(1,data_mouse{:});
% %zscore
% data_mouse = zscore(data_mouse,[],1);
% 
% nclass_pt = [1,5,10,15,20,25];
% nhold_pt = [0.1:0.05:0.5];
% nclass_nt = [1,5,10,15,20,25];
% nhold_nt = [0.1:0.05:0.5];
% paramsweep = combvec(nclass_pt,nhold_pt,nclass_nt,nhold_nt);
% %search parameter space
% rho = NaN(1,size(paramsweep,2));
% auc_pt = NaN(1,size(paramsweep,2));
% auc_nt = NaN(1,size(paramsweep,2));
% for j = 1:size(paramsweep,2)
%     if mod(j,10)==0
%        fprintf('\n\t working on parameter set %d',j);
%     end
%     [pt,~,~,~,~,~,allstats] = Phenotype_LOOCV(data_behavioral,imaged_mice,'verbose',0,'num_classifiers',paramsweep(1,j),'holdout',paramsweep(2,j));
%     allstats = [allstats{:}];
%     auc_pt(j) = nanmean([allstats.AUC]);
%     [nt,~,~,~,allstats] = DistanceToHyperplane(data_mouse,grp'+1,paramsweep(3,j),0,paramsweep(4,j)); 
%     allstats = [allstats{:}];
%     auc_nt(j) = nanmean([allstats.AUC]);
%     [rho(j),~] = corr(pt(:),nt(:)); %whole correlation    
% end
% 
% %%
% [~,b] = maxk(rho,10);
% paramsweep(:,b)
