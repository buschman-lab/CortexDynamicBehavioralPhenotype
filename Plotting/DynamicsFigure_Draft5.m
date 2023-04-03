%DynamicsFigure
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
fp = fig_params_vpa; 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\Dynamics';

% Load refitting data
file_list = GrabFiles(['\w*','chunk','\w*'],0,{data_dir});
% file_list = GrabFiles(['\w*','dynamics3','\w*'],0,{data_dir});
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
[nt, Y_null, betas, betas_null, allstats,bias] = DistanceToHyperplane(data_mouse,grp'+1,5,0,0.4);
% [nt, Y_null, betas, betas_null, allstats,bias] = DistanceToHyperplane(data_mouse,grp'+1,20,0,0.4);
[betas, idx] = sort(betas,'ascend');
data_mouse_sorted = data_mouse(:,idx);

%load the basis motifs
temp = load([data_dir filesep 'basismotifs.mat'],'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W = temp.W_basis;
W(:,temp.noise_clusters,:)=[];

[P,M,T] = size(W);
W = W(:,idx,:);

%Full Axis
[pt,~,~,~,~,~,~,betas_pt,~,~,bias_pt] = Phenotype(data_dir,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Social Axis
[pt_social,~,~,~,~,~,~,betas_social,~,~,bias_social] = Phenotype(data_dir,'behavior_col',2:5,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Developmental Axis
[pt_dev,~,~,~,~,~,~,betas_dev,~,~,bias_dev] = Phenotype(data_dir,'behavior_col',6:8,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Motor Phenotype Axis
[pt_motor,~,~,~,~,~,~,betas_motor,~,~,bias_motor] = Phenotype(data_dir,'behavior_col',9:13,'verbose',0,'num_classifiers',5,'holdout',0.5);

fp = fig_params_vpa;

%preload behavioral data
[data_behavioral_all,~,testName]= GetBehavioralData();

%just get the behavioral data for our mice of interest
imaged_mice = unique(mouse);
data_behavioral = NaN(numel(imaged_mice),size(data_behavioral_all,2)-1);
for i = 1:numel(imaged_mice)
    data_behavioral(i,:) = data_behavioral_all(data_behavioral_all(:,1)==imaged_mice(i),2:end);
end
data_behavioral(:,end) = [];

%% Ridge regression 
close all; rng('default')
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%first sort by idx
data_mouse_sorted = data_mouse(:,idx);
data_mouse_sorted = zscore(data_mouse_sorted,[],1);
X = data_mouse_sorted;
Y = zscore(pt_motor)';
k=logspace(-3,0,100);
nFold = 4;
%social, cross validate to find appropriate values
cvp = cvpartition(size(X,1),'KFold',nFold);
RMSE = ones(nFold,numel(k));
for i=1:nFold
   %train
   B = ridge(Y(cvp.training(i)),X(cvp.training(i),:),k,0);
   %get rmse for each k 
   for j = 1:numel(k)
       Yhat = B(1,j) + X(cvp.test(i),:)*B(2:end,j);
       RMSE(i,j) = sqrt(mean((Y(cvp.test(i))-Yhat).^2));
   end   
end
figure; hold on; 
plot(k,nanmean(RMSE))
set(gca,'xscale','linear')
figure
plot(k,B,'LineWidth',2)
grid on 
xlabel('Ridge Parameter') 
ylabel('Standardized Coefficient') 
title('Ridge Trace') 

%% use k = 0.35
close all; 
rng('default')
lambda=0.4;
Y = zscore(pt_motor)';
B = ridge(Y,X,lambda,0);
%predict
Yhat = B(1) + X*B(2:end);
Rsq = 1-sum((Y-Yhat).^2)/((length(Y)-1) * var(Y));
MSE = mean((Y-Yhat).^2);
rho = corr(Y,Yhat);

%compare with permuted
Rsq_perm = ones(1,1000);
MSE_perm = ones(1,1000);
rho_perm = ones(1,1000);
B_perm = ones(size(B,1),1000);
for i = 1:1000
   Y_perm = Y(randperm(size(Y,1),size(Y,1)),:);
   B_perm(:,i) = ridge(Y_perm,X,lambda,0);
   %predict
   Yhat = B_perm(1,i) + X*B_perm(2:end,i);
%    Rsq_perm(i) = 1-sum((Y-Yhat).^2)/((length(Y)-1) * var(Y)); 
%    MSE_perm(i) = mean((Y-Yhat).^2);
%    rho_perm(i) = corr(Y,Yhat);
   Rsq_perm(i) = 1-sum((Y_perm-Yhat).^2)/((length(Y_perm)-1) * var(Y_perm)); 
   MSE_perm(i) = mean((Y_perm-Yhat).^2);
   rho_perm(i) = corr(Y_perm,Yhat);   
end
% 
% figure; hold on; 
% histogram(Rsq_perm,'BinWidth',0.05)
% yval = get(gca,'ylim');
% line([Rsq,Rsq],yval,'color','r','linestyle','--')
% title(sprintf('%0.2g',sum(Rsq_perm>=Rsq)/numel(Rsq_perm)))
% 
% figure; hold on; 
% histogram(MSE_perm,'BinWidth',0.05)
% yval = get(gca,'ylim');
% line([MSE,MSE],yval,'color','r','linestyle','--')
% title(sprintf('%0.2g',sum(MSE_perm<=MSE)/numel(MSE_perm)))
% 
% figure; hold on; 
% histogram(rho_perm,'BinWidth',0.05)
% yval = get(gca,'ylim');
% line([rho,rho],yval,'color','r','linestyle','--')
% title(sprintf('%0.2g',sum(rho_perm>=rho)/numel(rho_perm)))

%confirm significant Betas.
pbeta = ones(1,size(B,1)-1);
for i = 1:size(B,1)-1
    temp = abs((B(1)+B(i+1)));
    temp_perm = abs([B_perm(1,:)+B_perm(i+1,:),temp]);
    pbeta(i)=sum(temp_perm>=temp)/numel(temp_perm);
end
bar(B(2:end));
set(gca,'xtick',1:size(B,1)-1,'xticklabels',pbeta)
set(gca,'xtick',1:size(B,1)-1,'xticklabels',pbeta,'xticklabelrotation',90)

%% 
rng('default')
%get the contribution of each variable to the ridge
delta_Rsq = ones(100,size(X,2));
for i = 1:100
   for j = 1:size(X,2)
       X_perm = X;
       X_perm(:,j) = X(randperm(size(X,1),size(X,1)),j);       
       B_perm = ridge(Y,X_perm,lambda,0);
       %predict
       Yhat = B_perm(1) + X_perm*B_perm(2:end);
%        Yhat = X_perm*B_perm;

       delta_Rsq(i,j) = 1-sum((Y-Yhat).^2)/((length(Y)-1) * var(Y)); 
   end
end
delta_Rsq = (Rsq-nanmean(delta_Rsq));
bar(B(2:end));
set(gca,'xtick',1:size(B,1)-1,'xticklabels',delta_Rsq,'xticklabelrotation',45)

%%



%% lasso
rng('default')
[B,FitInfo] = lasso(Y,pt_social','CV',10,'Alpha',0.5);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend
figure;
bar(B(:,idxLambdaMinMSE))

[B,FitInfo] = lasso(Y,pt_motor','CV',10,'Alpha',0.5);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend
figure;
bar(B(:,idxLambdaMinMSE))

[B,FitInfo] = lasso(Y,pt_dev','CV',10,'Alpha',0.5);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
lassoPlot(B,FitInfo,'PlotType','CV');
legend('show') % Show legend
figure;
bar(B(:,idxLambdaMinMSE))
%%
%fitlim
mdl = fitlm(Y,pt_social')
stemp = table2array(mdl.Coefficients(:,1));
stemp = stemp(2:end);
figure; bar(stemp);

mdl = fitlm(Y,pt_motor')
mtemp = table2array(mdl.Coefficients(:,1));
mtemp = mtemp(2:end);
figure; bar(mtemp);

mdl = fitlm(Y,pt_dev')
dtemp = table2array(mdl.Coefficients(:,1));
dtemp = dtemp(2:end);
figure; bar(dtemp);


%%
% X = score(:,1:5);
[A,B,r,U,V, stats] = canoncorr(X,Y);
nPerm = 1000;
wilks=stats.Wilks;

%just correlate
r=corr(Y,X)

%%
%try out pls
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,Y,2,'cv',5);

%permutation test
mse_perm = ones(2,size(MSE,2),nPerm);
for i = 1:nPerm
   Y_perm = Y(randperm(size(Y,1),size(Y,1)),:);
   [~,~,~,~,~,~,MSE_perm,~] = plsregress(X,Y_perm,2,'cv',4);
%    [~,~,r_perm(i,:),~,~,stats_perm] = canoncorr(X,Y_perm);
   mse_perm(:,:,i) = MSE_perm;
end
prctile(mse_perm,5,3)

%% kfold cross validation for strength and significance
rng('default')
nfold = 10;
A_cv = NaN(size(A,1),size(A,2),nfold);
B_cv = NaN(size(B,1),size(B,2),nfold);
cv_rho = ones(nfold,size(A,2));
nSamp = floor(0.4*20);
for i = 1:nfold
   temp = ones(size(data_mouse_sorted,1),1);
   temp(randperm(numel(temp),nSamp))=0;
   %compute cca without that animal
   [A_cv(:,:,i),B_cv(:,:,i),~,~,~,~] = canoncorr(X(temp==1,:),Y(temp==1,:));
   %project onto test data
   U_temp = (X(temp==0,:)-mean(X))*A_cv(:,:,i);
   V_temp = (Y(temp==0,:)-mean(Y))*B_cv(:,:,i);
   %get the correlation
   for j = 1:size(A,2)
      cv_rho(i,j) = corr(U_temp(:,j),V_temp(:,j));
   end
end
nanmean(cv_rho.^2)

%% now get the bootstrapped distribution of correlation levels for significance
rng('default');
cv_rho_perm = ones(1000,size(A,2));
for cur_perm = 1:1000
    Y_perm = Y(randperm(size(Y,1),size(Y,1)),:);
    cv_rho_temp = ones(nfold,size(A,2));
    for i = 1:nfold
       temp = ones(size(data_mouse_sorted,1),1);
       temp(randperm(numel(temp),nSamp))=0;
       %compute cca without that animal
       [A_loocv_perm,B_loocv_perm,~,~,~,~] = canoncorr(X(temp==1,:),Y_perm(temp==1,:));
       %project onto test data
       U_temp_perm = (X(temp==0,:)-mean(X))*A_loocv_perm;
       V_temp_perm = (Y(temp==0,:)-mean(Y))*B_loocv_perm;
       %get the correlation       
       for j = 1:size(A,2)
          cv_rho_temp(i,j) = corr(U_temp_perm(:,j),V_temp_perm(:,j));
       end       
    end    
    %get the strength
    cv_rho_perm(cur_perm,:) = nanmean(cv_rho_temp.^2);
end

prctile(cv_rho_perm,95,1)


%% assess stability with subtract 1 jack-knife
%LOOCV
U_temp = U;
V_temp = V;
A_loocv = NaN(size(A,1),size(A,2),size(data_mouse_sorted,1));
B_loocv = NaN(size(B,1),size(B,2),size(data_mouse_sorted,1));
for i = 1:size(data_mouse_sorted,1)
   temp = ones(size(data_mouse_sorted,1),1);
   temp(i)=0;
   %compute cca without that animal
   [A_loocv(:,:,i),B_loocv(:,:,i),~,~,~,~] = canoncorr(X(temp==1,:),Y(temp==1,:));
   U_temp(i,:) = (X(temp==0,:)-mean(X))*A_loocv(:,:,i);
   V_temp(i,:) = (Y(temp==0,:)-mean(Y))*B_loocv(:,:,i);   
end

%get the absolute strength of the loadings
A_loocv = nanvar(abs(A_loocv),[],3);
B_loocv = nanvar(abs(B_loocv),[],3);

%compare to the distribution of average strengths of permuted data. could
%also do this with variance
rng('default')
A_loocv_perm  = NaN(size(A_loocv,1),size(A_loocv,2),nPerm);
B_loocv_perm  = NaN(size(B_loocv,1),size(B_loocv,2),nPerm);
for j = 1:nPerm
    U_perm = U;
    V_perm = V;
    Y_perm = Y(randperm(size(Y,1),size(Y,1)),:);
    A_perm = NaN(size(A,1),size(A,2),size(data_mouse_sorted,1));
    B_perm = NaN(size(B,1),size(B,2),size(data_mouse_sorted,1));
    for i = 1:size(data_mouse_sorted,1)
       temp = ones(size(data_mouse_sorted,1),1);
       temp(i)=0;       
       %compute cca without that animal
       [A_perm(:,:,i),B_perm(:,:,i),~,~,~,~] = canoncorr(X(temp==1,:),Y_perm(temp==1,:));  
    end    
    A_loocv_perm(:,:,j) = nanvar(abs(A_perm),[],3);
    B_loocv_perm(:,:,j) = nanvar(abs(B_perm),[],3);
end

%get 95% of the strengths of the loadings
A_loocv_perm = prctile(A_loocv_perm,95,3);
B_loocv_perm = prctile(B_loocv_perm,95,3);

%get canonical variates where all loadings are more stable than permuted distribution
stable_idx = sum((A_loocv<A_loocv_perm)==0)+sum((B_loocv<B_loocv_perm)==0)==0;

%










p_perm_loocv = NaN(1,numel(r));
r_perm_mean_loocv = NaN(1,numel(r));
r_ci_mean_loocv_low = NaN(1,numel(r));
r_ci_mean_loocv_high = NaN(1,numel(r));
r_loocv_diff =  NaN(1,numel(r));
for i = 1:numel(r)
   p_perm_loocv(i) = sum(cat(1,r_loocv_perm(:,i),r_loocv(i))>=r_loocv(i))/(nPerm+1);
   r_perm_mean_loocv(i) = fisherInverse(nanmean(fisherZ(r_loocv_perm(:,i))));
   temp = bootci(1000,@nanmean,fisherZ(r_loocv_perm(:,i)));
   r_ci_mean_loocv_low(i) = fisherInverse(temp(1));
   r_ci_mean_loocv_high(i) = fisherInverse(temp(2)); 
   r_loocv_diff(i) = r_loocv(i)-r_perm_mean_loocv(i);
end


%find ones where both are significant
A_validated = A(:,1:2);
B_validated = B(:,1:2);
U_validated = U(:,1:2);
V_validated = V(:,1:2);

A_loocv_validated = nanmean(A_loocv(:,1:2,:),3);
B_loocv_validated = B_loocv(:,1:2,:);

cca_weights_loocv = nanmean(B_loocv_validated,3)'*coef(:,1:num_pcs)'+repmat(mu,size(B_loocv_validated,2),1);
cca_weights_behav_loocv = A_loocv_validated';

%Project B back into full space
cca_weights = B_validated'*coef(:,1:num_pcs)'+repmat(mu,size(B_validated,2),1); %the mu is ~0
cca_weights_behav = A_validated';

%Project A back into full space. Flip the social and developmental weights
%because their encoding betas as assigned to group ID (saline = negative).
%By Flipping them you match the direction of the original behavioral tests
%(e.g. increase social approach, increased social novelty. 
cca_weights_behav_full = NaN(2,numel(testName));
cca_weights_behav_full(:,1:4) = -1*A_validated(1,:)'*betas_social'+bias_social;
cca_weights_behav_full(:,8:12) = A_validated(2,:)'*betas_motor'+bias_motor;
cca_weights_behav_full(:,5:7) = -1*A_validated(3,:)'*betas_dev'+bias_dev;

%do the same thing for the cross validated cca weights
cca_weights_behav_full_loocv = NaN(2,numel(testName));
cca_weights_behav_full_loocv(:,1:4) = -1*(A_loocv_validated(1,:)'*betas_social'+bias_social);
cca_weights_behav_full_loocv(:,8:12) = A_loocv_validated(2,:)'*betas_motor'+bias_motor;
cca_weights_behav_full_loocv(:,5:7) = -1*(A_loocv_validated(3,:)'*betas_dev'+bias_dev);

%% Plot the neural data across animals
figure; hold on;
h = 2.4;
[n,m] = size(data_mouse_sorted);
imagesc(data_mouse_sorted');
caxis([-1.5,1.5]);
colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));
colorbar
set(gca,'xlim',[0.5, n+0.5],'ylim',[0.5 m+0.5],'units','centimeters',...
    'position',[3 3 5.5 h],'xtick',[1,n],'ytick',[1 m]);
xlabel('Animals')
ylabel('Motifs')
fp.FormatAxes(gca)  
box on

figure; hold on;
[n,m] = size(data_behavioral);
imagesc(data_behavioral');
caxis([-1.5,1.5]);
colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));
colorbar
set(gca,'xlim',[0.5, n+0.5],'ylim',[0.5 m+0.5],'units','centimeters',...
    'position',[3 3 5.5 (h/4)*3],'xtick',[1,n],'ytick',[1 m]);
xlabel('Animals')
ylabel('Assays')
fp.FormatAxes(gca)  
box on

%PCs
figure; hold on;
h = 2.4;
[n,m] = size(B');
imagesc(B);
caxis([-1.5,1.5]);
colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));
colorbar
set(gca,'xlim',[0.5, n+0.5],'ylim',[0.5 m+0.5],'units','centimeters',...
    'position',[3 3 5.5 h/2],'xtick',[1,n],'ytick',[1 m]);
xlabel('Mice')
ylabel('PCs')
fp.FormatAxes(gca)  
box on

%Weightings Matrices
figure; hold on;
h = 2.4;
[n,m] = size(data_mouse_sorted);
imagesc(data_mouse_sorted');
caxis([-1.5,1.5]);
colormap(gca,flipud(redgreencmap(256,'Interpolation','linear')));
colorbar
set(gca,'xlim',[0.5, n+0.5],'ylim',[0.5 m+0.5],'units','centimeters',...
    'position',[3 3 5.5 h],'xtick',[1,n],'ytick',[1 m]);
xlabel('Animals')
ylabel('Motifs')
fp.FormatAxes(gca)  
box on

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','zscoredmotifandbehavioraldata',savedir,1); close all


%% Plot the Choice of PC
figure; hold on; 
y = cumsum(explained(1:15));
x = 1:numel(y);
plot(x,y,'color','k','linewidth',2);
plot(x(num_pcs),y(num_pcs),'color','r','linestyle','none','Marker','o','markersize',10,'linewidth',2)
ylabel({'Explained','Variance (%)'})
xlabel('PCs')
set(gca,'units','centimeters','position',[3 2 3.5 3.5]);
fp.FormatAxes(gca)      


%% Plot the Beta Weights of the Significant Cannonical Variable Pairs
coef_lim_behav = {[-0.25,1],[-1,1]}; %[floor(min(cca_weights_behav(:))),ceil(max(cca_weights_behav(:)))];
coef_lim_behav_full = [-0.6,1.25]; 
coef_lim_motifs = [-0.6,0.6];
for i = 1:size(cca_weights,1)
    %motif weightings
    figure; hold on; 
    barh(cca_weights(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6])
    line([0,0],[-1,numel(cca_weights(i,:))+1],'color','k','linewidth',1.5)
    set(gca,'ytick',[1:numel(cca_weights(i,:))],'ydir','reverse')
    set(gca,'units','centimeters','position',[3,4,2,6],'YAxisLocation','right')
    set(gca,'ylim',[0.5,numel(cca_weights(i,:))+0.5]);
    set(gca,'xlim',[coef_lim_motifs]);
    xlabel({'Canonical','Coefficient'})
    ylabel('Motifs');
    fp.FormatAxes(gca)      
    
    %Behavioral weightings    
    figure; hold on; 
    bar(cca_weights_behav(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6],'FaceAlpha',0.5)     
    line([-1,numel(cca_weights_behav(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',1:numel(cca_weights_behav(i,:)),'xticklabels',{'Social','Motor','Dev.'},'xticklabelrotation',90)
    set(gca,'units','centimeters','position',[3,4,6,1.75],'YAxisLocation','left','XAxisLocation','top')
    set(gca,'xlim',[0.5,numel(cca_weights_behav(i,:))+0.5]);
    ylabel({'Canonical','Coefficient'})
    set(gca,'ylim',[coef_lim_behav{i}]);
    fp.FormatAxes(gca) 

    %Full Behavioral weightings    
    figure; hold on; 
    bar(cca_weights_behav_full(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6],'FaceAlpha',0.5)     
    line([-1,numel(cca_weights_behav_full(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',1:numel(cca_weights_behav_full(i,:)),'xticklabels',testName,'xticklabelrotation',90)
    set(gca,'units','centimeters','position',[3,2,6,1.75],'YAxisLocation','left','XAxisLocation','top')
    set(gca,'xlim',[0.5,numel(cca_weights_behav_full(i,:))+0.5]);
    ylabel({'Canonical','Coefficient'})
    set(gca,'ylim',[coef_lim_behav_full]);
    fp.FormatAxes(gca) 
        
    %Cannonical Correlation
    figure; hold on; 
    lm = fitlm(U_validated(:,i),V_validated(:,i)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(U_validated(grp==0,i),V_validated(grp==0,i),'.','markersize',25,'color',fp.c_sal)
    plot(U_validated(grp==1,i),V_validated(grp==1,i),'.','markersize',25,'color',fp.c_vpa)
    legend off
    set(gca,'xlim',[round(min(U_validated(:,i)),1)-0.1,round(max(U_validated(:,i)),1)+0.1])
    set(gca,'ylim',[round(min(V_validated(:,i)),1)-0.1,round(max(V_validated(:,i)),1)+0.1])
    ylabel('V')
    xlabel('U')
    [rho,~] = corr(U_validated(:,i),V_validated(:,i)); %correlation
    set(gca,'units','centimeters','position',[3 2 6 6]);
    fp.SetTitle(gca,{sprintf('rho=%0.2g',rho)});
    fp.FormatAxes(gca)      
    
     
    handles = get(groot, 'Children');
    fp.SaveFigs(handles,'-svg',sprintf('CCA%d',i),savedir,1); close all
end

%% Plot the Beta Weights of the Significant Cannonical Variable Pairs
coef_lim_behav = {[-0.25,1],[-1,1]}; %[floor(min(cca_weights_behav(:))),ceil(max(cca_weights_behav(:)))];
coef_lim_behav_full = [-1,1]; 
coef_lim_motifs = [-0.6,0.6];
for i = 1:size(cca_weights,1)
    %motif weightings
    figure; hold on; 
    yyaxis left
    bar(cca_weights(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.2 0.2 0.2],'FaceAlpha',0.25)
    line([-1,numel(cca_weights(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',[1:numel(cca_weights(i,:))],'xdir','normal')
    set(gca,'units','centimeters','position',[3,4,6,2])
    set(gca,'xlim',[0.5,numel(cca_weights(i,:))+0.5]);
    set(gca,'ylim',[coef_lim_motifs]);
    
    yyaxis right
    bar(cca_weights_loocv(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0 0 1],'FaceAlpha',0.15)
    line([-1,numel(cca_weights_loocv(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',[1:numel(cca_weights_loocv(i,:))],'xdir','normal')
    set(gca,'units','centimeters','position',[3,4,6,2])
    set(gca,'xlim',[0.5,numel(cca_weights_loocv(i,:))+0.5]);
    set(gca,'ylim',[-.1 .1]);
    
    xlabel({'Canonical','Coefficient'})
    ylabel('Motifs');
    fp.FormatAxes(gca)      
    
    %Behavioral weightings    
    figure; hold on; 
    bar(cca_weights_behav(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6],'FaceAlpha',0.5)     
    line([-1,numel(cca_weights_behav(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',1:numel(cca_weights_behav(i,:)),'xticklabels',{'Social','Motor','Dev.'},'xticklabelrotation',90)
    set(gca,'units','centimeters','position',[3,4,6,1.75],'YAxisLocation','left','XAxisLocation','top')
    set(gca,'xlim',[0.5,numel(cca_weights_behav(i,:))+0.5]);
    ylabel({'Canonical','Coefficient'})
    set(gca,'ylim',[coef_lim_behav{i}]);
    fp.FormatAxes(gca) 

    %Full Behavioral weightings    
    figure; hold on; 
    bar(cca_weights_behav_full(i,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6],'FaceAlpha',0.5)     
    line([-1,numel(cca_weights_behav_full(i,:))+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',1:numel(cca_weights_behav_full(i,:)),'xticklabels',testName,'xticklabelrotation',90)
    set(gca,'units','centimeters','position',[3,2,6,1.75],'YAxisLocation','left','XAxisLocation','top')
    set(gca,'xlim',[0.5,numel(cca_weights_behav_full(i,:))+0.5]);
    ylabel({'Canonical','Coefficient'})
    set(gca,'ylim',[coef_lim_behav_full]);
    fp.FormatAxes(gca) 
        
    %Cannonical Correlation
    figure; hold on; 
    lm = fitlm(U_validated(:,i),V_validated(:,i)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(U_validated(grp==0,i),V_validated(grp==0,i),'.','markersize',25,'color',fp.c_sal)
    plot(U_validated(grp==1,i),V_validated(grp==1,i),'.','markersize',25,'color',fp.c_vpa)
    legend off
    set(gca,'xlim',[round(min(U_validated(:,i)),1)-0.1,round(max(U_validated(:,i)),1)+0.1])
    set(gca,'ylim',[round(min(V_validated(:,i)),1)-0.1,round(max(V_validated(:,i)),1)+0.1])
    ylabel('V')
    xlabel('U')
    [rho,~] = corr(U_validated(:,i),V_validated(:,i)); %correlation
    set(gca,'units','centimeters','position',[3 2 6 6]);
    fp.SetTitle(gca,{sprintf('rho=%0.2g',rho)});
    fp.FormatAxes(gca)      
    
     
    handles = get(groot, 'Children');
    fp.SaveFigs(handles,'-svg',sprintf('CCAloocv%d',i),savedir,1); close all
end

%% Plot the motif beta weights on top of each other
coef_lim_motifs = [-0.6,0.6];
figure; hold on; 
barh(cca_weights(1,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0 0 1],'FaceAlpha',0.4)
barh(cca_weights(2,:),'FaceColor','flat','EdgeAlpha',0,'CData',[1 0 0],'FaceAlpha',0.4)
line([0,0],[-1,numel(cca_weights(1,:))+1],'color','k','linewidth',1.5)
set(gca,'ytick',[1:numel(cca_weights(1,:))],'ydir','reverse')
set(gca,'units','centimeters','position',[3,4,2,6],'YAxisLocation','right')
set(gca,'ylim',[0.5,numel(cca_weights(1,:))+0.5]);
set(gca,'xlim',coef_lim_motifs);
rho = corr(cca_weights(1,:)',cca_weights(2,:)','type','Spearman');
xlabel(sprintf('rho %0.2g',rho))
[~, temp] = maxk(abs(cca_weights(1,:)-cca_weights(2,:)),6);
title(sprintf('%d,',temp));
ylabel('Motifs');
fp.FormatAxes(gca)  

%same thing just for the loocv weights
coef_lim_motifs = [-0.25,0.25];
figure; hold on; 
barh(cca_weights_loocv(1,:),'FaceColor','flat','EdgeAlpha',0,'CData',[0 0 1],'FaceAlpha',0.4)
barh(cca_weights_loocv(2,:),'FaceColor','flat','EdgeAlpha',0,'CData',[1 0 0],'FaceAlpha',0.4)
line([0,0],[-1,numel(cca_weights_loocv(1,:))+1],'color','k','linewidth',1.5)
set(gca,'ytick',[1:numel(cca_weights_loocv(1,:))],'ydir','reverse')
set(gca,'units','centimeters','position',[3,4,2,6],'YAxisLocation','right')
set(gca,'ylim',[0.5,numel(cca_weights_loocv(1,:))+0.5]);
set(gca,'xlim',coef_lim_motifs);
rho = corr(cca_weights_loocv(1,:)',cca_weights_loocv(2,:)','type','Spearman');
xlabel(sprintf('rho %0.2g',rho))
[~, temp] = maxk(abs(cca_weights_loocv(1,:)-cca_weights_loocv(2,:)),6);
title(sprintf('%d,',temp));
ylabel('Motifs');
fp.FormatAxes(gca)  

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','comparedmotifcoefs',savedir,1); close all

%% Validate the subspaces on totally different animals MOTOR
rec_names = {'431-10-17','432-10-17','432-10-18'};
high_threshold = [0.7,0.2, 0.2]; 
low_threshold = [0.4, 0.05, 0.05];
dur = 15*13;
AUC = cell(1,numel(rec_names));
betas_behavior_all = cell(1,numel(rec_names));
behavioral_trace = cell(1,numel(rec_names));

for cur_rec = 1:numel(rec_names)
    if ~exist([rec_names{cur_rec},'behavioraldata.mat'],'file')
        [H_long,behavior] = CompileLongerRecordings(rec_names{cur_rec});
    else
        temp = load([rec_names{cur_rec},'behavioraldata.mat'],'data','features_downsampled','features_downsampled_nologzscore');
        H_long = temp.data;
        motorbehavior_raw = temp.features_downsampled_nologzscore(:,3);
    end

    %order, detrend and smooth
    H_long_avg = H_long(idx,:); 
    H_long_sem = H_long_avg;
    for j = 1:size(H_long_avg,1) 
       temp = detrend(tensor_convolve(nanmean(W(:,j,:),1),H_long_avg(j,:)));
       H_long_avg(j,:) = movmean(temp,dur); 
       H_long_sem(j,:) = movstd(temp,dur)/sqrt(dur); 
    end

    %same for behavior
    motorbehavior_raw = movmean(motorbehavior_raw,dur);

    %save for later
    behavioral_trace{cur_rec} = motorbehavior_raw;    

    rng('default');
    %binarize motor behavior
    low = find(motorbehavior_raw<low_threshold(cur_rec));
    high = find(motorbehavior_raw>=high_threshold(cur_rec));

    %match sample numbers (fewer motor epochs)                
    num_init = 1;
    betas_behavior = zeros(16,num_init);
    %preallocate since DistanceToHyperplane resets the rng
    data_fold = cell(1,num_init);
    for i = 1:num_init
        len_min = min(numel(low),numel(high));
        idx_low = low(randperm(numel(low),len_min));
        idx_high = high(randperm(numel(high),len_min));
        data_fold{i} = cat(1,H_long_avg(:,idx_low)',H_long_avg(:,idx_high)');
    end
    allstats = {};
    for i = 1:num_init
        [~, ~, betas_behavior(:,i), ~, temp,~] = DistanceToHyperplane(data_fold{i},cat(1,ones(len_min,1),ones(len_min,1)*2),5,0,0.4);  
        allstats = cat(2,allstats,temp);
    end
    betas_behavior = nanmean(betas_behavior,2);   
    betas_behavior_all{cur_rec} = betas_behavior;

    %Plot the classifier AUCs
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
    
    %plot the correlation between the betas of motifs
    figure; hold on; 
    plot(cca_weights(1,:)',betas_behavior,'.','markersize',20,'color',[0.5 0.5 0.5])
    lm = fitlm(cca_weights(1,:)',betas_behavior);
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2);
%     [rho,~] = corr(betas_behavior,cca_weights(1,:)','type','Pearson');
    [rho,~] = corr(betas_behavior,cca_weights_loocv(1,:)','type','Pearson');
    y = betas_behavior;
    % Compare to permuted brain-behavior axes

    Y_temp = (score(:,1:num_pcs));
    X_temp = zscore([pt_social',pt_motor',pt_dev']);
    %Test Significance with bootstrapped distribution
    rng('default')
    nPerm=1000;
    r_loocv_perm  = NaN(nPerm,2);
    for j = 1:nPerm
        Y_perm = Y_temp(randperm(size(Y_temp,1),size(Y_temp,1)),:);
        for i = 1:size(data_mouse_sorted,1)
           temp = ones(size(data_mouse_sorted,1),1);
           temp(i)=0;       
           %compute cca without that animal
           [~,B_perm(:,:,i),~,~,~,~] = canoncorr(X_temp(temp==1,:),Y_perm(temp==1,:));
        end

        %Compute the average motor and social correlation with beta weights
        temp_weight = nanmean(B_perm(:,1:2,:),3)'*coef(:,1:num_pcs)'+repmat(mu,size(B_perm(:,1:2,:),2),1);
        r_loocv_perm(j,1) = corr(temp_weight(1,:)',y,'type','Pearson');
        r_loocv_perm(j,2) = corr(temp_weight(2,:)',y,'type','Pearson');    
    end

    if rho >=0
        p = sum([r_loocv_perm(:,1),rho]>=rho)/(nPerm+1);
    else
        p = sum([r_loocv_perm(:,1),rho]<=rho)/(nPerm+1);
    end

    set(gca,'units','centimeters','position',[3 2 3 3]);
    xlabel('Motor Subspace Weights');
    ylabel({'Motor Activity','Decoder Betas'});
    fp.SetTitle(gca,{sprintf('Pearson rho=%0.2g p=%0.2g',rho,p)});
    fp.FormatAxes(gca)      
    legend off

    figure; hold on; 
    plot(cca_weights(2,:)',betas_behavior,'.','markersize',20,'color',[0.5 0.5 0.5])
    lm = fitlm(cca_weights(2,:)',betas_behavior);
    p = plot(lm,'markersize',0.1,'color','k','linewidth',2);
%     [rho,~] = corr(betas_behavior,cca_weights(2,:)','type','Pearson');
    [rho,~] = corr(betas_behavior,cca_weights_loocv(2,:)','type','Pearson');
    set(gca,'units','centimeters','position',[3 2 3 3]);
    xlabel('Social Subspace Weights');
    ylabel({'Motor Activity','Decoder Betas'});
   
    if rho >=0
        p = sum([r_loocv_perm(:,2),rho]>=rho)/(nPerm+1);
    else
        p = sum([r_loocv_perm(:,2),rho]<=rho)/(nPerm+1);
    end

    fp.SetTitle(gca,{sprintf('Pearson rho=%0.2g p=%0.2g',rho,p)});
    fp.FormatAxes(gca)      
    legend off   

    figure('position',[680   480   560   498]); hold on; 
    b = bar(betas_behavior,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
    line([-1,numel(betas_behavior)+1],[0,0],'color','k','linewidth',1.5)
    set(gca,'xtick',[1:1:16],'XTickLabelRotation',0,'xdir','normal')
    set(gca,'units','centimeters','position',[3,2,4,2],'YAxisLocation','left')
    set(gca,'xlim',[0.5,numel(betas_behavior)+0.5]);
    set(gca,'ylim',[round(min(betas_behavior),1)-0.1,round(max(betas_behavior),1)+0.1],'ytick',[round(min(betas_behavior),1)-0.1,round(max(betas_behavior),1)+0.1])
    ylabel('\bf \beta')
    fp.FormatAxes(gca);
 
    
    figure; hold on; 
    T = numel(motorbehavior_raw);
    plot(motorbehavior_raw,'color',[0.5 0.5 0.5],'linewidth',1.25);
    ylabel('Limb Speed'); 
    xlabel('Time (min)')
    line([0,T],[high_threshold(cur_rec),high_threshold(cur_rec)],'color',[0.75 0 0],'linewidth',2,'linestyle','--')
    line([0,T],[low_threshold(cur_rec),low_threshold(cur_rec)],'color',[0 0 0.75],'linewidth',2,'linestyle',':')
    set(gca,'xtick',linspace(0,T,5),'xticklabel',linspace(0,T/13/60,5),'xlim',[0, T])
    set(gca,'units','centimeters','position',[2 2 4 2]);
    fp.FormatAxes(gca)    
    
   
    handles = get(groot, 'Children');
    fp.SaveFigs(handles,'-svg',sprintf('SubspaceVsMotorActivity_%s',rec_names{cur_rec}),savedir,1); close all     
end
% close all

%% Now normalize the beta weights and combine
Y = cellfun(@(x) (x-nanmean(x))/std(x), betas_behavior_all,'UniformOutput',0);
rng('default');
N = size(cca_weights,2);
Y = cat(1,Y{:});
% X = repmat(cca_weights(1,:),1,3)';
X = repmat(cca_weights_loocv(1,:),1,3)';

figure; hold on; 
plot(X(1:N),Y(1:N),'.','markersize',fp.m_markersize*1.5,'color',[0 0 0])
plot(X(N+1:N*2),Y(N+1:N*2),'.','markersize',fp.m_markersize*1.5,'color',[0.3 0.3 0.3])
plot(X(N*2+1:end),Y(N*2+1:end),'.','markersize',fp.m_markersize*1.5,'color',[0.6 0.6 0.6])
lm = fitlm(X,Y);
p = plot(lm,'markersize',0.1,'color','k','linewidth',2);
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper

rho = corr(X,Y,'type','Pearson');

% Compare to permuted brain-behavior axes
Y_temp = (score(:,1:num_pcs));
X_temp = zscore([pt_social',pt_motor',pt_dev']);
%Test Significance with bootstrapped distribution
rng('default')
nPerm=1000;
r_loocv_perm  = NaN(nPerm,2);
for j = 1:nPerm
    Y_perm = Y_temp(randperm(size(Y_temp,1),size(Y_temp,1)),:);
    for i = 1:size(data_mouse_sorted,1)
       temp = ones(size(data_mouse_sorted,1),1);
       temp(i)=0;       
       %compute cca without that animal
       [~,B_perm(:,:,i),~,~,~,~] = canoncorr(X_temp(temp==1,:),Y_perm(temp==1,:));
    end
    
    %Compute the average motor and social correlation with beta weights
    temp_weight = nanmean(B_perm(:,1:2,:),3)'*coef(:,1:num_pcs)'+repmat(mu,size(B_perm(:,1:2,:),2),1);
    r_loocv_perm(j,1) = corr(repmat(temp_weight(1,:),1,3)',Y,'type','Pearson');
    r_loocv_perm(j,2) = corr(repmat(temp_weight(2,:),1,3)',Y,'type','Pearson');    
end

if rho >=0
    p = sum([r_loocv_perm(:,2),rho]>=rho)/(nPerm+1);
else
    p = sum([r_loocv_perm(:,2),rho]<=rho)/(nPerm+1);
end

set(gca,'units','centimeters','position',[3 2 5 5]);
xlabel('Motor Subspace Weights');
ylabel({'Motor Activity','Decoder Betas'});
fp.SetTitle(gca,{sprintf('Pearson rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)      
legend off

% X = repmat(cca_weights(2,:),1,3)';
X = repmat(cca_weights_loocv(2,:),1,3)';
figure; hold on; 
plot(X(1:N),Y(1:N),'.','markersize',fp.m_markersize*1.5,'color',[0 0 0])
plot(X(N+1:N*2),Y(N+1:N*2),'.','markersize',fp.m_markersize*1.5,'color',[0.3 0.3 0.3])
plot(X(N*2+1:end),Y(N*2+1:end),'.','markersize',fp.m_markersize*1.5,'color',[0.6 0.6 0.6])
lm = fitlm(X,Y);
p = plot(lm,'markersize',0.1,'color','k','linewidth',2);
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %bounds lower
p(4).Color = [0 0 0]; %bound upper
rho = corr(X,Y,'type','Pearson');

if rho >=0
    p = sum([r_loocv_perm(:,2),rho]>=rho)/(nPerm+1);
else
    p = sum([r_loocv_perm(:,2),rho]<=rho)/(nPerm+1);
end

set(gca,'units','centimeters','position',[3 2 5 5]);
xlabel('Social Subspace Weights');
ylabel({'Motor Activity','Decoder Betas'});
fp.SetTitle(gca,{sprintf('Pearson rho=%0.2g p=%0.2g',rho,p)});
fp.FormatAxes(gca)      
legend off

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','SubspaceVsMotorActivity_all',savedir,1); close all 


%% Bar plot of the average Beta Weights
Y = cellfun(@(x) (x-nanmean(x))/std(x), betas_behavior_all,'UniformOutput',0);
Y = cat(2,Y{:});
Y_mean = nanmean(Y,2);
% Y_mean = cca_weights(1,:);
figure('position',[680   480   560   498]); hold on; 
b = bar(Y_mean,'FaceColor','flat','EdgeAlpha',0,'CData',[0.6 0.6 0.6]);
line([-1,numel(Y_mean)+1],[0,0],'color','k','linewidth',1.5)
set(gca,'xtick',[1:1:16],'XTickLabelRotation',0,'xdir','normal')
set(gca,'units','centimeters','position',[3,2,4,2],'YAxisLocation','left')
set(gca,'xlim',[0.5,numel(Y_mean)+0.5]);
set(gca,'ylim',[round(min(Y_mean),1)-0.1,round(max(Y_mean),1)+0.1],'ytick',[round(min(Y_mean),1)-0.1,round(max(Y_mean),1)+0.1])
ylabel('\bf \beta')
fp.FormatAxes(gca);


handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','SubspaceVsMotorActivity_allbetas',savedir,1); close all 



