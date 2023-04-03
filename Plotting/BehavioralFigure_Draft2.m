%Behavioral Figure

savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\BehavioralFigures';
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
warning off; [zscoreData, rawData, testName] = GetBehavioralData(); warning on;

fp = fig_params_vpa; 
%%
close all
%Results of Individual Behavioral Tests
group = rawData(:,end);
rng('default');
sal_pos = 1+rand(sum(group==0),1)/2-0.25;
vpa_pos = 2+rand(sum(group==1),1)/2-0.25;
for i = 2:numel(testName)+1 %start at the second column       
    figure('position',[680   558   240   271]); hold on;
    vpa = rawData(rawData(:,end)==1,i);
    sal = rawData(rawData(:,end)==0,i);
    plot(sal_pos,sal,'.','markersize',fp.m_markersize/1.5,'color',fp.c_sal)
    plot(vpa_pos,vpa,'.','markersize',fp.m_markersize/1.5,'color',fp.c_vpa)      
    errorbar(1,nanmean(sal),sem(sal),'LineWidth',1.5,'Color','k');
    errorbar(2,nanmean(vpa),sem(vpa),'LineWidth',1.5,'Color','k');
    [p,~] = ranksum(sal,vpa);
    temp = max(cat(1,sal,vpa));
    line([1,2],[temp+temp*0.025,temp+temp*0.025],'linewidth',1.5,'color','k');
    text(1.5, temp+temp*0.05,sprintf('%0.2g',p),'Rotation',0,'FontSize',fp.sig_fontsize,'HorizontalAlignment','center');
    set(gca,'Clipping','off','box','off');   
    set(gca,'units','centimeters','position',[2 2 1.25 3])
    set(gca,'xtick',[1,2],'xticklabel',{'SAL','VPA'},'xticklabelrotation',90,'xlim',[0.25 2.75])
    ylabel(testName{i-1})
    fp.FormatAxes(gca);
    set(gca,'LineWidth',2)
end    
%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','IndividualBehavioralTests',savedir,1); close all
%% Unsupervised projection 

rng('default'); %for reproducibility
%compute in PC space
[coef, score, latent, ~, explained, mu] = pca(zscoreData(:,1:end-1));

grp = zscoreData(:,end);
%Fit the hyperplane with classifier
partition = [];
partition.training = true(numel(grp),1); %train on entire dataset
partition.test = true(numel(grp),1); %doesn't matter as you aren't testing this classifier
[~, observed, ~, ~] = SVMClassifier_Binary([score(:,1:3),grp+1],partition,'holdout',0.5,'optimize',0,'kernel','linear','featureselect','none','nshuf',0,'verbose',0); 
mdl = observed.Classifier;
X = score(:,1:3);
group = zscoreData(:,end)+1;
%Gather support vectors from ClassificationSVM struct
sv =  mdl.SupportVectors;
%set step size for finer sampling
d =0.05;
%generate grid for predictions at finer sample rate
[x, y, z] = meshgrid(min(X(:,1)):d:max(X(:,1)),...
min(X(:,2)):d:max(X(:,2)), min(X(:,3)):d:max(X(:,3)));
xGrid = [x(:),y(:),z(:)];
%get scores, f
[ ~ , f] = predict(mdl,xGrid);
%reshape to same grid size as the input
f = reshape(f(:,2), size(x));
% Assume class labels are 1 and 0 and convert to logical
t = logical(group);
%plot data points, color by class label
figure; hold on; 
plot3(score(zscoreData(:,end)==0,1),score(zscoreData(:,end)==0,2),score(zscoreData(:,end)==0,3),'.','markersize',20,'color',fp.c_sal)
plot3(score(zscoreData(:,end)==1,1),score(zscoreData(:,end)==1,2),score(zscoreData(:,end)==1,3),'.','markersize',20,'color',fp.c_vpa)
% load unscaled support vectors for plotting
% plot3(sv(:, 1), sv(:, 2), sv(:, 3), 'ro','markersize',9,'linewidth',1.5);
%plot decision surface
[faces,verts,~] = isosurface(x, y, z, f, 0, x);
patch('Vertices', verts, 'Faces', faces, 'FaceColor','r','edgecolor', 'none', 'FaceAlpha', 0.15);
grid on
campos([-375.1951   29.9322   15.2842])
xlabel('PC1'); ylabel('PC2'); zlabel('PC3'); 
set(gca,'zlim',[-4 2],'units','centimeters','position',[3,2,6,6])
fp.FormatAxes(gca);
box on

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','unsupervised_behavioralaxis',savedir,1); close all

%%
%Global Phenotype Axis
[y_img,~,y,~,grp,~,~,betas,~,allstats,~,] = Phenotype(data_dir,'verbose',0,'num_classifiers',5,'holdout',0.5);
AUC_all_general  = allstats{end};
AUC = cellfun(@(x) x.AUC, allstats(1:end-1),'UniformOutput',0);
AUC_all = nanmean(cat(2,AUC{:}));
% y = (y-min(y))/(max(y)-min(y));

%Social Axis
[y_social_img,~,y_social,~,~,~,~,~,~,temp] = Phenotype(data_dir,'behavior_col',2:5,'verbose',0,'num_classifiers',5,'holdout',0.5);
AUC_social_general  = temp{end};
AUC = cellfun(@(x) x.AUC, temp(1:end-1),'UniformOutput',0);
AUC_social = nanmean(cat(2,AUC{:}));
% y_social = (y_social-min(y_social))/(max(y_social)-min(y_social));

%Developmental Axis
[y_dev_img,~,y_dev,~,~,~,~,~,~,temp] = Phenotype(data_dir,'behavior_col',6:8,'verbose',0,'num_classifiers',5,'holdout',0.5);
AUC_develop_general  = temp{end};
AUC = cellfun(@(x) x.AUC, temp(1:end-1),'UniformOutput',0);
AUC_develop = nanmean(cat(2,AUC{:}));
% y_dev = (y_dev-min(y_dev))/(max(y_dev)-min(y_dev));

%Motor Phenotype Axis
[y_motor_img,~,y_motor,~,~,~,~,~,~,temp] = Phenotype(data_dir,'behavior_col',9:13,'verbose',0,'num_classifiers',5,'holdout',0.5);
AUC_motor_general  = temp{end};
AUC = cellfun(@(x) x.AUC, temp(1:end-1),'UniformOutput',0);
AUC_motor = nanmean(cat(2,AUC{:}));

% y_motor = (y_motor-min(y_motor))/(max(y_motor)-min(y_motor));


%%
sm = 20;
lg = 30;
rng('default')
figure('position',[680   558   669   420]); hold on;
x = rand(numel(y),1)*0.5;
temp = cat(1,y,y_social,y_motor,y_dev);
plot(temp(:,grp==1),(x(grp==1)+[-0.5,1,2,3])','color',[fp.c_sal 0.2],'linewidth',1.5)
plot(temp(:,grp==2),(x(grp==2)+[-0.5,1,2,3])','color',[fp.c_vpa 0.2],'linewidth',1.5)
plot(y(grp==1),x(grp==1)-0.5,'.','markersize',lg,'color',fp.c_sal);    
plot(y_dev(grp==1),x(grp==1)+3,'.','markersize',sm,'color',fp.c_sal);   
plot(y_motor(grp==1),x(grp==1)+2,'.','markersize',sm,'color',fp.c_sal);   
plot(y_social(grp==1),x(grp==1)+1,'.','markersize',sm,'color',fp.c_sal);   
plot(y(grp==2),x(grp==2)-0.5,'.','markersize',lg,'color',fp.c_vpa);       
plot(y_dev(grp==2),x(grp==2)+3,'.','markersize',sm,'color',fp.c_vpa);   
plot(y_motor(grp==2),x(grp==2)+2,'.','markersize',sm,'color',fp.c_vpa);   
plot(y_social(grp==2),x(grp==2)+1,'.','markersize',sm,'color',fp.c_vpa);

%circle the imaged mice
x_temp = NaN(1,numel(y_img));
for i = 1:numel(y_img)
    x_temp(i) = x(find(y==y_img(i)));
end
plot(y_img,x_temp-0.5,'.','markersize',15,'color',[1 1 0]);  
plot(y_social_img,x_temp+1,'.','markersize',10,'color',[1 1 0]);  
plot(y_dev_img,x_temp+3,'.','markersize',10,'color',[1 1 0]);  
plot(y_motor_img,x_temp+2,'.','markersize',10,'color',[1 1 0]);  


line([0,0],[-1,4],'linewidth',2,'linestyle','--','color','k')
set(gca,'ylim',[-0.6,3.6],'ytick',[-0.5 1 2 3],'YTickLabel',{'Global','Social','Motor','Develop.'},'ydir','reverse','YTickLabelRotation',50)
set(gca,'xtick',[min(get(gca,'xtick')),max(get(gca,'xtick'))],'xticklabel',{'SAL-like','VPA-like'})
set(gca,'units','centimeters','position',[4,1,12,10])
fp.FormatAxes(gca);

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','behavioralaxis',savedir,1); close all


%% plot the betas
figure('position',[680   480   560   498]); hold on; 
barh(betas,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6])
line([0,0],[-1,numel(betas)+1],'color','k','linewidth',1.5)
set(gca,'ytick',[1:numel(betas)],'Yticklabels',testName,'YTickLabelRotation',0)
set(gca,'units','centimeters','position',[3,2,3,6],'YAxisLocation','right')
set(gca,'ylim',[0.5,numel(betas)+0.5]);
xlabel('\bf \beta')
fp.FormatAxes(gca);

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','behavioralbetas',savedir,1); close all


%% plot the ROCs
X = cellfun(@(x) x.X(:,1), allstats,'UniformOutput',0);
X = cat(2,X{:});
Y = cellfun(@(x) x.Y(:,1), allstats,'UniformOutput',0);
Y = cat(2,Y{:});
AUC = cellfun(@(x) x.AUC, allstats,'UniformOutput',0);
AUC = cat(2,AUC{:});

figure; hold on; 
plot(X,Y,'color',[0.5 0.5 0.5,0.3],'linewidth',1.5)
plot(nanmean(X,2),nanmean(Y,2),'color','k','linewidth',2)
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title(sprintf('AUC = %0.2g',nanmean(AUC)),'FontWeight','normal');
set(gca,'units','centimeters','position',[2 2 4 4]);
fp.FormatAxes(gca)
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','behavioralAUC',savedir,1); close all

%% plot the correlation between the axes
temp = cat(1,y,y_social,y_motor,y_dev)';
[rho,p] = corr(temp);

figure; hold on; 
rho(triu(rho)>1)=0;
imagesc(rho,[-0.5 1])
set(gca,'ydir','reverse','xdir','normal')
cmap = cat(1,[1 1 1],parula);
colormap(cmap)
colorbar
set(gca,'xtick',1:4,'xticklabel',{'G','S','M','D'},'XAxisLocation','top')
set(gca,'ytick',1:4,'yticklabel',{'G','S','M','D'})
set(gca,'ylim',[0.5,4.5],'xlim',[0.5,4.5]);
box on
for i = 1:size(rho,1)
    for j= 1:size(rho,2)        
        text(j,i,sprintf('%0.2g',rho(i,j)),'HorizontalAlignment','center');
    end
end
set(gca,'units','centimeters','position',[2 2 3 3])
fp.FormatAxes(gca);


figure; hold on; 
rho(triu(rho)>1)=0;
imagesc(rho,[-0.5 1])
set(gca,'ydir','reverse','xdir','normal')
cmap = cat(1,[1 1 1],parula);
colormap(cmap)
colorbar
set(gca,'xtick',1:4,'xticklabel',{'G','S','M','D'},'XAxisLocation','top')
set(gca,'ytick',1:4,'yticklabel',{'G','S','M','D'})
set(gca,'ylim',[0.5,4.5],'xlim',[0.5,4.5]);
box on
for i = 1:size(rho,1)
    for j= 1:size(rho,2)        
        text(j,i,sprintf('%0.2g',p(i,j)),'HorizontalAlignment','center');
    end
end
set(gca,'units','centimeters','position',[2 2 3 3])
fp.FormatAxes(gca);

%%
handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','behavioralaxiscorr',savedir,1); close all


% set(gca,'xtick',[min(get(gca,'xtick')),max(get(gca,'xtick'))],'xticklabel',{'SAL-like','VPA-like'})
%%

% %% Compare to PCA
% 
% [coeff, score, latent, tsquared, explained, mu] = pca(zscoreData(:,2:end-1)); 
% figure; hold on 
% plot(explained);
% xlabel('PC')
% ylabel('PEV')
% fp.SetTitle(gca,'Variance Explained By Each PC');
% fp.FormatAxes(gca);
% 
% % grp = zscoreData(:,end)+1;
% figure; hold on; 
% plot(score(grp==1,1),score(grp==1,2),'.','markersize',15,'color',fp.c_sal)
% plot(score(grp==2,1),score(grp==2,2),'.','markersize',15,'color',fp.c_vpa)
% fp.SetTitle(gca,'Behavioral Data in 2D PCA Space');
% fp.FormatAxes(gca);
% 
% figure; hold on; 
% plot3(score(grp==1,1),score(grp==1,2),score(grp==1,3),'.','markersize',15,'color',fp.c_sal)
% plot3(score(grp==2,1),score(grp==2,2),score(grp==2,3),'.','markersize',15,'color',fp.c_vpa)
% fp.SetTitle(gca,'Behavioral Data in 3D PCA Space');
% fp.FormatAxes(gca);
% 
% 
% figure; hold on; 
% [rho, p] = corr(y',score);
% plot(rho,'.','markersize',15,'color',[0.5 0.5 0.5])
% for i = 1:numel(rho)
%     text(i,rho(i)+0.1, sprintf('p=%0.2g',p(i)),'FontSize',12,'Rotation',90)
% end
% fp.SetTitle(gca,'Correlation Between Phenotype Axis and Each PC');
% fp.FormatAxes(gca);


   

   
   
   
   
   
   
   
   
   
   
   
   
   
   
   