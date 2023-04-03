%EntropyFigure
data_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL';
fp = fig_params_vpa; 
savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\Figures\SoloVSPaired';
if ~exist(savedir)
    mkdir(savedir);
end

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
[betas, idx] = sort(betas,'ascend');
% data_mouse_sorted = data_mouse(:,idx);

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
[pt_dev,~,~,~,~,~,~,betas_dev,~,~,bias_dev]= Phenotype(data_dir,'behavior_col',6:8,'verbose',0,'num_classifiers',5,'holdout',0.5);
%Motor Phenotype Axis
[pt_motor,~,~,~,~,~,~,betas_motor,~,~,bias_motor]= Phenotype(data_dir,'behavior_col',9:13,'verbose',0,'num_classifiers',5,'holdout',0.5);

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


%% CCA Analysis Full Space
%Find the cannoncial variables using the full model; 
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});    
data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
data_mouse = cat(1,data_mouse{:});
%first sort by idx
data_mouse_sorted = data_mouse(:,idx);
data_mouse_sorted = zscore(data_mouse_sorted,[],1);

rng('default')
[coef, score, ~, ~, explained, mu] = pca(data_mouse_sorted);
num_pcs = 8; %find(cumsum(explained)>85); %Elbow_pt(cat(1,0,cumsum(explained)));
Y = (score(:,1:num_pcs));
X = zscore([pt_social',pt_motor',pt_dev']);
[A,B,r,U,V, stats] = canoncorr(X,Y);
nPerm = 1000;

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

r_loocv  = NaN(1,numel(r));
for i = 1:numel(r)
   r_loocv(i) = corr(U_temp(:,i),V_temp(:,i));   
end

%Test Significance with bootstrapped distribution
rng('default')
r_loocv_perm  = NaN(nPerm,numel(r));
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
       U_perm(i,:) = (X(temp==0,:)-mean(X(temp==1,:)))*A_perm(:,:,i);
       V_perm(i,:) = (Y(temp==0,:)-mean(Y(temp==1,:)))*B_perm(:,:,i);   
    end    
    
    for i = 1:numel(r)
       r_loocv_perm(j,i) = corr(U_perm(:,i),V_perm(:,i));   
    end
end

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

%% PLOTTING 
%get the total weightings
solo_activity = data_mouse(:,idx)'*100;

% Load the Social Recording data
file_list = GrabFiles(['\w*','refit','\w*'],0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social'});
data_social = cellfun(@(x) load(x,'loadings','mouseID','stats_refit'),file_list,'UniformOutput',0);
mouse_social_id = cellfun(@(x) x.mouseID, data_social,'UniformOutput',0);
mouse_social_id = [mouse_social_id{:}];

%number of vpa and saline animals used for social recs. 


%average explained variance across 2 min fits
temp = ones(1,numel(data_social));
for i = 1:numel(data_social)
   temp(i) = nanmean(cellfun(@(x) x.pev, data_social{1,i}.stats_refit));
end
%average per mouse
temp = arrayfun(@(x) nanmean(temp(mouse_social_id==x))*100, unique(mouse_social_id));
soc_pev = nanmean(temp);
soc_pev_ci = bootci(1000,@nanmean,temp);


%average per mouse
social_activity = cellfun(@(x)  nanmean(x.loadings(idx,:),2)*100,data_social,'UniformOutput',0);
social_activity = cat(3,social_activity{:});
social_activity = arrayfun(@(x) nanmean(social_activity(:,:,mouse_social_id==x),3), unique(mouse_social_id), 'UniformOutput', 0);
%combine animals
social_activity = cat(2,social_activity{:});

delta_activity = social_activity - solo_activity;

%% Generate the reference distribution
target_axis = [pt_motor',pt_social']; %reflip social so 
Y_temp = (score(:,1:num_pcs));
X_temp = zscore([pt_social',pt_motor',pt_dev']);
%Test Significance with reference distribution
rng('default')
nPerm=1000;
perm_stat_ref = NaN(nPerm,2);
for j = 1:nPerm
    Y_perm = Y_temp(randperm(size(Y_temp,1),size(Y_temp,1)),:);
    for i = 1:size(data_mouse_sorted,1)
       temp = ones(size(data_mouse_sorted,1),1);
       temp(i)=0;       
       %compute cca without that animal
       [~,B_perm(:,:,i),~,~,~,~] = canoncorr(X_temp(temp==1,:),Y_perm(temp==1,:));
    end

    %Compute the average motor and social correlation with beta weights
    temp_weight = nanmean(B_perm(:,1:2,:),3)'*coef(:,1:num_pcs)';%+repmat(mu,size(B_perm(:,1:2,:),2),1);
    temp_weight = (temp_weight*delta_activity)';
    for q = 1:2
        perm_stat_ref(j,q) = corr(target_axis(:,q),temp_weight(:,q));    
    end
end

delta_activity_axis = (cca_weights_loocv*delta_activity)';
label = [{'motor'},{'social'}];
for i = 1:numel(label)
    figure; hold on; 
    lm = fitlm(target_axis(:,i),delta_activity_axis(:,i)); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
    p(1).Color = [1 1 1]; %poitns
    p(2).Color = [0 0 0]; %fit
    p(3).Color = [0 0 0]; %bounds lower
    p(4).Color = [0 0 0]; %bound upper
    p(2).LineWidth = 2;
    p(3).LineWidth = 1.25;
    p(4).LineWidth = 1.25;
    plot(target_axis(grp==0,i),delta_activity_axis(grp==0,i),'.','markersize',25,'color',fp.c_sal)
    plot(target_axis(grp==1,i),delta_activity_axis(grp==1,i),'.','markersize',25,'color',fp.c_vpa)
%     plot(target_axis(:,i),temp(:,i),'.','markersize',25,'color',[0.25 0.25 0.25])
    legend off
    set(gca,'xlim',[round(min(target_axis(:,i)),1)-0.1,round(max(target_axis(:,i)),1)+0.1])
    set(gca,'ylim',[min(delta_activity_axis(:,i))-0.01, max(delta_activity_axis(:,i))+0.01])
    set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'Negative','Positive'});
    ylabel('\Delta Activity','Interpreter','tex')
    xlabel('Phenotype')
    rho = corr(target_axis(:,i),delta_activity_axis(:,i),'type','pearson'); %correlation
    set(gca,'units','centimeters','position',[3 2 5,5]);       
    
    if rho<0
        pref = sum(cat(1,perm_stat_ref(:,i),rho)<=rho)/(1000+1);
    else
        pref = sum(cat(1,perm_stat_ref(:,i),rho)>=rho)/(1000+1);
    end       
            
    %reference distribution
    rhoref = fisherInverse(nanmean(fisherZ(perm_stat_ref(:,i))));
    ciref = fisherInverse(bootci(1000,@nanmean,fisherZ(perm_stat_ref(:,i))));
    diffref = rho-nanmean(perm_stat_ref(:,i));
      
    %also do with within group shuffle
    perm_stat = NaN(1,nPerm);
    for j = 1:nPerm              
       within = delta_activity_axis(:,i)';
       temp = within(grp==0);
       within(grp==0) = temp(randperm(numel(temp),numel(temp))); 
       temp = within(grp==1);
       within(grp==1) = temp(randperm(numel(temp),numel(temp))); 
       perm_stat(j) = corr(within',target_axis(:,i));
    end
   
    if rho<0
        pval_withinall = sum(cat(2,perm_stat,rho)<=rho)/(1000+1);
    else
        pval_withinall = sum(cat(2,perm_stat,rho)>=rho)/(1000+1);
    end

    rho_within = fisherInverse(nanmean(fisherZ(perm_stat)));
    ciwithin = fisherInverse(bootci(1000,@nanmean,fisherZ(perm_stat)));
    diffwithin = rho-nanmean(perm_stat);
    
    
    fp.SetTitle(gca,{sprintf('%s rho=%0.2g pref=%0.2g \n rhoref=%0.2g, diffref=%0.2g \n CIref=%0.2g %0.2g \n rhowin = %0.2g diffwin=%0.2g pwin=%0.2g \n CIwin=%0.2g %0.2g',...
        label{i}, rho, pref, rhoref, diffref, ciref(1), ciref(2), rho_within, diffwithin, pval_withinall, ciwithin(1), ciwithin(2))})  
    
    fp.FormatAxes(gca)   
end
%%

%% Plot groupwise comparison 

N = size(delta_activity_axis,1);
pos = rand(N,1)/2-0.25;
for j = 1:2
    figure; hold on; 
    sal = delta_activity_axis(grp==0,j);
    vpa = delta_activity_axis(grp==1,j);
    b = bar(1,nanmean(sal),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
    b.CData(1,:) = fp.c_sal;
    errorbar(1,nanmean(sal),sem(sal),'LineWidth',1.25,'Color','k');
%     plot(pos(grp==0)+1,sal,'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
    b = bar(2,nanmean(vpa),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
    b.CData(1,:) = fp.c_vpa;
    errorbar(2,nanmean(vpa),sem(vpa),'LineWidth',1.25,'Color','k');
%     plot(pos(grp==1)+2,vpa,'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)    
    if nanmean(sal)-nanmean(vpa)>0
        p = ranksum(sal,vpa,'tail','right');
    else
        p = ranksum(sal,vpa,'tail','left');
    end
    text(1.5,max(get(gca,'ylim')),sprintf('diff=%0.2g p=%0.2g',nanmean(sal)-nanmean(vpa),p),'VerticalAlignment','bottom','HorizontalAlignment','center')
    
    if nanmean(sal)>0
        psal = signrank(sal,0,'tail','right');
    else
        psal = signrank(sal,0,'tail','left');
    end   
    if nanmean(vpa)>0
        pvpa = signrank(vpa,0,'tail','right');
    else
        pvpa = signrank(vpa,0,'tail','left');
    end    
    
    
    title(sprintf('sal=%0.2g,sem%0.2g,p=%0.2g \n vpa=%0.2g,sem%0.2g,p=%0.2g',nanmean(sal),sem(sal),psal,nanmean(vpa),sem(vpa),pvpa),'fontweight','normal','fontsize',6,'Position',[2.5 2.5]);
    %format
    set(gca,'xlim',[0.5, 2.5],'units','centimeters','position',[4,3,1.25,2])
    fp.FormatAxes(gca)
end

%% Plot the change between the environment as a whole

%Compare
N = size(delta_activity_axis,1);
figure; hold on; 
xloc = [1,2];
pos = repmat((rand(N,1)/2-0.25),1,size(delta_activity_axis,2)) + xloc;
col = [0 0 1; 1 0 0];
for i = 1:size(delta_activity_axis,2)    
    b = bar(xloc(i),nanmean(delta_activity_axis(:,i)),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]); 
    b.CData(1,:) = col(i,:);
    plot(pos(:,i),delta_activity_axis(:,i),'.','markersize',fp.m_markersize*1.25,'color',col(i,:))
    plot(pos(grp==0,i),delta_activity_axis(grp==0,i),'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
    plot(pos(grp==1,i),delta_activity_axis(grp==1,i),'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)
end
%link mice
for i = 1:N
%     plot(pos(i,:),delta_activity(i,:),'linewidth',1.25,'color',[0.5 0.5 0.5 0.4])
    if grp(i)==0
        plot(pos(i,:),delta_activity_axis(i,:),'linewidth',1.5,'color',[fp.c_sal 0.4])
    else
        plot(pos(i,:),delta_activity_axis(i,:),'linewidth',1.5,'color',[fp.c_vpa 0.4])
    end
end
%paired permutation testing
yloc = repmat(ceil(max(delta_activity_axis(:))),1,2);
% start + vs -
p = signrank(delta_activity_axis(:,1),delta_activity_axis(:,2),'tail','right');
line([1,2],yloc,'linewidth',2,'color','k')
text(1.5,yloc(1),sprintf('p=%0.2g',p),'VerticalAlignment','bottom','HorizontalAlignment','center')

%format
set(gca,'ylim',[floor(min(delta_activity_axis(:))-0.1),ceil(yloc(1)+0.1)],'xlim',[0.5 2.5])
ylabel({'Weighted Percent','Explained Variance'})
set(gca,'units','centimeters','position',[4,3,2.5,5])
fp.FormatAxes(gca)

%Additional significance testing
p_vs_zero_motor = signrank(delta_activity_axis(:,1),0);
p_vs_zero_social = signrank(delta_activity_axis(:,2),0);


%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','ComparingSoloToPaired',savedir,1); close all
%%


%%%% Now Investigate the Dynamics
%% Load the data by chunks and concatenate
shift = 13*30;
dur = 13*60*2;
cd('Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution')
if ~exist('ConcatenatedMotifData_2Min_5holdout4.mat','file')
    file_list_ordered = cell(1,numel(file_list));
    file_list_ordered(1:2:end) = file_list(2:2:end);
    file_list_ordered(2:2:end) = file_list(1:2:end);

    %get the original data as well
    original_data_dir = {'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution'};

    %look through and combine each into the single whole recording
    mice_list = zeros(1,numel(file_list_ordered)/6);
    loadings = cell(1,numel(file_list_ordered)/6);
    H_raw_all = cell(1,numel(file_list_ordered)/6);
    H_convolve_all =cell(1,numel(file_list_ordered)/6); 
    for i = 1:(numel(file_list_ordered)/6)
       fprintf('\t\nWorking on Recording %d of %d',i,(numel(file_list_ordered)/6))
       [mice_list(i), ~] = MouseNumFromFileName(file_list_ordered((i-1)*6+1));

       H = cell(1,6);
       W = cell(1,6);
       X = cell(1,6);
       for j = 1:6 %6 is the number of blocks per recording
           temp = load(file_list_ordered{ (((i-1)*6)+j)},'w','H');
           H{j} = temp.H(idx,:);
           W{j} = temp.w(:,idx,:);    

           %get recording name to link with original data
           [~, name] = fileparts(file_list_ordered{(i-1)*6+j});
           if ~isempty(regexp(name,'train','ONCE')) %if training then add that
               name = erase(name,'train');
               temp = GrabFiles([name, '.mat'],0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution'});
               temp = load(temp{1},'data_train');
               X{j} = temp.data_train;

           else %add testing
               name = erase(name,'test');
               temp = GrabFiles([name, '.mat'],0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution'});
               temp = load(temp{1},'data_test');
               X{j} = temp.data_test;
           end
       end       
       X = cat(2,X{:});
       %loop through each motif and compute it's pev sliding window style
       pev = cell(1,numel(idx));
       H_convolve = zeros(numel(idx),size(X,2));
       for j = 1:numel(idx)
          Xhat = cell(1,6);  
          temp_H_convolve = cell(1,6); 
          for cur_chunk = 1:6
              Xhat{cur_chunk} = tensor_convolve(W{cur_chunk}(:,j,:),H{cur_chunk}(j,:));
              temp_H_convolve{cur_chunk} = tensor_convolve(nanmean(W{cur_chunk}(:,j,:),1),H{cur_chunk}(j,:));
          end
          Xhat = cat(2,Xhat{:});
          H_convolve(j,:) = cat(2,temp_H_convolve{:});
          pev{j} = SlidingWindowPEV(X,Xhat,dur,shift);
       end
       pev = cat(1,pev{:});
       loadings{i} = pev./sum(pev,1);%save off loadings
       H_raw_all{i} = cat(2,H{:}); %save of raw H's
       H_convolve_all{i} = H_convolve; %save of convolved H's
    end
    %save off
    save('ConcatenatedMotifData_2Min_5holdout4.mat','loadings','mice_list','idx','cca_weights','H_raw_all','H_convolve_all','dur','shift','-v7.3');
else
    load('ConcatenatedMotifData_2Min_5holdout4.mat','loadings','mice_list');
end

%% Plot the activity of the subspaces over time
subspace_activity = cellfun(@(x) (cca_weights_loocv*x*100), loadings, 'UniformOutput', 0);
%average per mice first
subspace_activity = cat(3,subspace_activity{:});
subspace_activity = arrayfun(@(x) nanmean(subspace_activity(:,:,mice_list==x),3), unique(mice_list), 'UniformOutput', 0);
subspace_activity = (cat(1,subspace_activity{:}));
% subspace_activity = zscore(cat(1,subspace_activity{:}),[],2);
motor_sub = subspace_activity(1:2:end,:)';
social_sub = subspace_activity(2:2:end,:)';

social_sub = social_sub-social_sub(1,:);
motor_sub = motor_sub-motor_sub(1,:);

ylimits = [-1.5,2.5];
T = size(social_sub,1);
f1 = figure; hold on; 
plot(social_sub,'color',[0.5 0.5 0.5,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(social_sub,1),nanmean(social_sub,2),sem(social_sub,2),'lineprops',{'color',[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(social_sub,2),'color',[0.5 0.5 0.5 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    temp = social_sub(i,:);
    if nanmean(temp(:))>=0
        [p(1,i)] = signrank(temp(:),0,'tail','right');
    else
        [p(1,i)] = signrank(temp(:),0,'tail','left');
    end
end
f2 = figure; hold on; 
line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
line([1,T+0.5],[-log10(0.05/21),-log10(0.05/21)],'linestyle',':','linewidth',1.5,'color','r');
plot(-log10(p(1,:)),'color','k','linewidth',1.5)
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[],'YAxisLocation','right','ytick',[])
set(gca,'ylim',[0,5])
% set(gca,'color',[0.95 0.95 0.95 1]);
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on
ylimits = [-0.6,0.3];
f3 = figure; hold on; 
plot(motor_sub,'color',[0.5 0.5 0.5,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(motor_sub,1),nanmean(motor_sub,2),sem(motor_sub,2),'lineprops',{'color',[0.5 0.5 0.5]},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(motor_sub,2),'color',[0.5 0.5 0.5 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    temp = motor_sub(i,:);
    if nanmean(temp(:))>=0
        [p(1,i)] = signrank(temp(:),0,'tail','right');
    else
        [p(1,i)] = signrank(temp(:),0,'tail','left');
    end
end
f4 = figure; hold on; 
line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
line([1,T+0.5],[-log10(0.05/21),-log10(0.05/21)],'linestyle',':','linewidth',1.5,'color','r');
plot(-log10(p(1,:)),'color','k','linewidth',1.5)
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[],'YAxisLocation','right','ytick',[])
set(gca,'ylim',[0,4])
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on
social_sub_solo = social_sub; 
motor_sub_solo = motor_sub;

% Now compare with social
subspace_activity = cellfun(@(x)  cca_weights_loocv*x.loadings(idx,:)*100,data_social,'UniformOutput',0);
subspace_activity = cat(3,subspace_activity{:});
subspace_activity = arrayfun(@(x) nanmean(subspace_activity(:,:,mouse_social_id==x),3), unique(mouse_social_id), 'UniformOutput', 0);
subspace_activity = (cat(1,subspace_activity{:}));

motor_sub = subspace_activity(1:2:end,:)';
social_sub = subspace_activity(2:2:end,:)';

social_sub = social_sub-social_sub(1,:);
motor_sub = motor_sub-motor_sub(1,:);

ylimits = [-1.5,2.5];
T = size(social_sub,1);
set(0, 'currentfigure', f1);
plot(social_sub,'color',[fp.c_nmf,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(social_sub,1),nanmean(social_sub,2),sem(social_sub,2),'lineprops',{'color',[fp.c_nmf]},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(social_sub,2),'color',[fp.c_nmf 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    temp = social_sub(i,:);
    if nanmean(temp(:))>=0
        [p(1,i)] = signrank(temp(:),0,'tail','right');
    else
        [p(1,i)] = signrank(temp(:),0,'tail','left');
    end
end
set(0, 'currentfigure', f2);
% line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
% line([1,T+0.5],[-log10(0.05/21),-log10(0.05/21)],'linestyle',':','linewidth',1.5,'color','r');
plot(-log10(p(1,:)),'color',fp.c_nmf,'linewidth',1.5)
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[],'YAxisLocation','right','ytick',[])
set(gca,'ylim',[0,5])
% set(gca,'color',[0.95 0.95 0.95 1]);
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on
ylimits = [-0.6,0.3];
set(0, 'currentfigure', f3);
plot(motor_sub,'color',[fp.c_nmf,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(motor_sub,1),nanmean(motor_sub,2),sem(motor_sub,2),'lineprops',{'color',fp.c_nmf},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(motor_sub,2),'color',[fp.c_nmf 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    temp = motor_sub(i,:);
    if nanmean(temp(:))>=0
        [p(1,i)] = signrank(temp(:),0,'tail','right');
    else
        [p(1,i)] = signrank(temp(:),0,'tail','left');
    end
end
set(0, 'currentfigure', f4);
line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
line([1,T+0.5],[-log10(0.05/21),-log10(0.05/21)],'linestyle',':','linewidth',1.5,'color','r');
plot(-log10(p(1,:)),'color',fp.c_nmf,'linewidth',1.5)
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[],'YAxisLocation','right','ytick',[])
set(gca,'ylim',[0,4])
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on

social_sub_paired = social_sub; 
motor_sub_paired = motor_sub;
%% Compare solo vs social temporal dynamics 
signrank(sum(social_sub_paired.^2)-sum(social_sub_solo.^2),0,'tail','left') 
nanmean(sum(social_sub_paired.^2))
nanmean(sum(social_sub_solo.^2))


%% Compare the change be the end of the recording in social and saline animals

ranksum(sum(social_sub_paired(:,grp==0).^2),0);
ranksum(sum(social_sub_paired(:,grp==1).^2),0)

temp = social_sub_paired(end,:);
nanmean(temp(grp==0))
signrank(temp(grp==0),0)
ci = bootci(100,@nanmean,temp(grp==0));
nanmean(temp(grp==1))
ci = bootci(100,@nanmean,temp(grp==1));
signrank(temp(grp==1),0)

perm_stat = NaN(1,nPerm);
for i = 1:nPerm
    temp_perm = temp(randperm(numel(temp),numel(temp)));
    perm_stat(i) = nanmean(temp_perm(grp==1))-nanmean(temp_perm(grp==0));
end
orig_val = nanmean(temp(grp==1))-nanmean(temp(grp==0));

pval = sum(perm_stat>=orig_val)/(nPerm+1);


%comapre change in social between groups
figure; hold on; 
pos = rand(20,1)/2-0.25;
sal = nanmean(social_sub(:,grp==0))';
vpa = nanmean(social_sub(:,grp==1))';
b = bar(1,nanmean(sal),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_sal;
errorbar(1,nanmean(sal),sem(sal),'LineWidth',1.25,'Color','k');
plot(pos(grp==0)+1,sal,'.','markersize',fp.m_markersize*1.25,'color',fp.c_sal)
b = bar(2,nanmean(vpa),'FaceColor','flat','FaceAlpha',0.25,'EdgeColor',[1 1 1]);
b.CData(1,:) = fp.c_vpa;
errorbar(2,nanmean(vpa),sem(vpa),'LineWidth',1.25,'Color','k');
plot(pos(grp==1)+2,vpa,'.','markersize',fp.m_markersize*1.25,'color',fp.c_vpa)    

%format
set(gca,'xlim',[0.5, 2.5],'units','centimeters','position',[4,3,2,5])
fp.FormatAxes(gca)


% handles = get(groot, 'Children');
% fp.SaveFigs(handles,'-svg','Timecourse2',savedir,1); close all

%% Plot the social activity over time
subspace_activity = cellfun(@(x)  cca_weights_loocv*x.loadings(idx,:)*100,data_social,'UniformOutput',0);
subspace_activity = cat(3,subspace_activity{:});
subspace_activity = arrayfun(@(x) nanmean(subspace_activity(:,:,mouse_social_id==x),3), unique(mouse_social_id), 'UniformOutput', 0);
subspace_activity = cat(1,subspace_activity{:});
motor_sub = subspace_activity(1:2:end,:)';
social_sub = subspace_activity(2:2:end,:)';

%center to the average of the saline animals
social_sub = social_sub-nanmean(social_sub(1,grp==0));
motor_sub = motor_sub-nanmean(motor_sub(1,grp==0));

ylimits = [-3,1.1];
T = size(social_sub,1);
figure; hold on; 
plot(social_sub(:,grp==0),'color',[fp.c_sal,0.2],'linewidth',1.25);
plot(social_sub(:,grp==1),'color',[fp.c_vpa,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(social_sub,1),nanmean(social_sub(:,grp==0),2),sem(social_sub(:,grp==0),2),'lineprops',{'color',[fp.c_sal]},'transparent',1,'patchSaturation',0.075);
shadedErrorBar(1:size(social_sub,1),nanmean(social_sub(:,grp==1),2),sem(social_sub(:,grp==1),2),'lineprops',{'color',[fp.c_vpa]},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(social_sub(:,grp==0),2),'color',[fp.c_sal 1],'linewidth',2);
p1 = plot(nanmean(social_sub(:,grp==1),2),'color',[fp.c_vpa 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    p(1,i) = nanmean(social_sub(i,grp==0))-nanmean(social_sub(i,grp==1));
    if nanmean(social_sub(i,grp==0))-nanmean(social_sub(i,grp==1))>=0
        [p(2,i)] = ranksum(social_sub(i,grp==0),social_sub(i,grp==1),'tail','right');
    else
        [p(2,i)] = ranksum(social_sub(i,grp==0),social_sub(i,grp==1),'tail','left');
    end  
end
figure; hold on; 
yyaxis left
plot(p(1,:),'color','k','linewidth',2);
set(gca,'ylim',[0.4,0.85])
yyaxis right
set(gca,'ylim',[0.75,2])
plot(-log10(p(2,:)),'color',[1 0 0 0.5],'linewidth',1.5)
line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[])
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on

ylimits = [-0.6,0.6];
figure; hold on; 
plot(motor_sub(:,grp==0),'color',[fp.c_sal,0.2],'linewidth',1.25);
plot(motor_sub(:,grp==1),'color',[fp.c_vpa,0.2],'linewidth',1.25);
line([0,size(motor_sub,1)+1],[0 0],'linestyle','--','color','k','linewidth',2)
shadedErrorBar(1:size(motor_sub,1),nanmean(motor_sub(:,grp==0),2),sem(motor_sub(:,grp==0),2),'lineprops',{'color',[fp.c_sal]},'transparent',1,'patchSaturation',0.075);
shadedErrorBar(1:size(motor_sub,1),nanmean(motor_sub(:,grp==1),2),sem(motor_sub(:,grp==1),2),'lineprops',{'color',[fp.c_vpa]},'transparent',1,'patchSaturation',0.075);
p1 = plot(nanmean(motor_sub(:,grp==0),2),'color',[fp.c_sal 1],'linewidth',2);
p1 = plot(nanmean(motor_sub(:,grp==1),2),'color',[fp.c_vpa 1],'linewidth',2);
set(gca,'xlim',[0.5 T+0.5])
set(gca,'xtick',[1,T],'xticklabel',[1,12])
set(gca,'ylim',ylimits);
xlabel('Time (min)')
ylabel('\Delta Activity','interpreter','tex')
set(gca,'units','centimeters','position',[4,3,3.5,5])
fp.FormatAxes(gca)
p = NaN(2,T);
for i = 1:T
    p(1,i) = nanmean(motor_sub(i,grp==0))-nanmean(motor_sub(i,grp==1));
    if nanmean(motor_sub(i,grp==0))-nanmean(motor_sub(i,grp==1))>=0
        [p(2,i)] = ranksum(motor_sub(i,grp==0),motor_sub(i,grp==1),'tail','right');
    else
        [p(2,i)] = ranksum(motor_sub(i,grp==0),motor_sub(i,grp==1),'tail','left');
    end  
end
figure; hold on; 
yyaxis left
plot(p(1,:),'color','k','linewidth',2);
set(gca,'ylim',[-0.2,-.1])
yyaxis right
set(gca,'ylim',[0.5,2])
plot(-log10(p(2,:)),'color',[1 0 0 0.5],'linewidth',1.5)
line([1,T+0.5],[-log10(0.05),-log10(0.05)],'linestyle',':','linewidth',1.5,'color','b');
set(gca,'units','centimeters','position',[4,3,3.5,1]);
set(gca,'xlim',[0.5 T+0.5],'xtick',[])
fp.FormatAxes(gca)
set(gca,'linewidth',1)
box on


%% Get the average motor and social activity without referencing

sal_motor = nanmean(motor_sub(:,grp==0));
sal_social = nanmean(social_sub(:,grp==0));
vpa_motor = nanmean(motor_sub(:,grp==1));
vpa_social = nanmean(social_sub(:,grp==1));
    
nanmean(sal_social);
ci = bootci(100,@nanmean,sal_social);
nanmean(vpa_social);
ci = bootci(100,@nanmean,vpa_social);

temp = nanmean(social_sub); 
perm_stat = NaN(1,nPerm);
for i = 1:nPerm
    temp_perm = temp(randperm(numel(temp),numel(temp)));
    perm_stat(i) = nanmean(temp_perm(grp==1))-nanmean(temp_perm(grp==0));
end
orig_val = nanmean(temp(grp==1))-nanmean(temp(grp==0));

sum([perm_stat,orig_val]<=orig_val)/(nPerm+1)


nanmean(sal_motor);
ci = bootci(100,@nanmean,sal_motor);
nanmean(vpa_motor);
ci = bootci(100,@nanmean,vpa_motor);

temp = nanmean(motor_sub); 
perm_stat = NaN(1,nPerm);
for i = 1:nPerm
    temp_perm = temp(randperm(numel(temp),numel(temp)));
    perm_stat(i) = nanmean(temp_perm(grp==1))-nanmean(temp_perm(grp==0));
end
orig_val = nanmean(temp(grp==1))-nanmean(temp(grp==0));

sum([perm_stat,orig_val]>=orig_val)/(nPerm+1)


%% Comparing within group 
temp = nanmean(social_sub_paired([end-3:end],:))-nanmean(social_sub_paired([1:3],:)); %similar results are see with block size ranging from 1-5. 
rho_both = corr(pt_social',temp');
perm_stat = PermutationStatistic(pt_social',temp', @(x,y) corr(x,y),1000);
if rho<0
    pval = sum([perm_stat,rho]<=rho)/(1000+1);
else
    pval = sum([perm_stat,rho]>=rho)/(1000+1);
end

%also do the within group together
nt_centered = temp;
nt_centered(grp==0) = nt_centered(grp==0)-nanmean(nt_centered(grp==0)); %0.08
nt_centered(grp==1) = nt_centered(grp==1)-nanmean(nt_centered(grp==1));

pt_centered = pt_social; 
pt_centered(grp==0) = pt_centered(grp==0)-nanmean(pt_centered(grp==0));
pt_centered(grp==1) = pt_centered(grp==1)-nanmean(pt_centered(grp==1));

rng('default');
nt_centered_perm = nt_centered;
nPerm = 1000;
perm_stat = NaN(1,nPerm);
for j = 1:nPerm
    temp_poop = nt_centered(grp==0);
    N = numel(temp_poop);
    nt_centered_perm(grp==0) = temp_poop(randperm(N,N)); 
    temp_poop = nt_centered(grp==1);
    N =  numel(temp_poop);
    nt_centered_perm(grp==1) = temp_poop(randperm(N,N));     
    perm_stat(j) = corr(nt_centered_perm',pt_centered');
end

rho_both = corr(pt_centered',nt_centered');

if rho_both<0
    pval_both = sum([perm_stat,rho_both]<=rho_both)/(1000+1);
else
    pval_both = sum([perm_stat,rho_both]>=rho_both)/(1000+1);
end

figure; hold on; 
lm = fitlm(pt_social,temp); p = plot(lm,'markersize',0.1,'color','k','linewidth',2); 
p(1).Color = [1 1 1]; %poitns
p(2).Color = [0 0 0]; %fit
p(3).Color = [0 0 0]; %fit
p(4).Color = [0 0 0]; %fit
p(2).LineWidth = 2;
p(3).LineWidth = 1.25;
p(4).LineWidth = 1.25;
plot(pt_social(grp==0),temp(grp==0),'.','markersize',25,'color',fp.c_sal)
plot(pt_social(grp==1),temp(grp==1),'.','markersize',25,'color',fp.c_vpa)
legend off
set(gca,'xlim',[round(min(pt_social),1)-0.1,round(max(pt_social),1)+0.1])
set(gca,'ylim',[min(temp)-0.01, max(temp)+0.01])
set(gca,'xtick',get(gca,'xlim'),'xticklabel',{'Negative','Positive'});
ylabel('\Delta Activity','Interpreter','tex')
title(sprintf('change in paired rho %0.2g p= %0.2g',rho_both,pval_both));
xlabel('Social Phenotype')
set(gca,'units','centimeters','position',[4,3,5,5]);
fp.FormatAxes(gca)


%%

handles = get(groot, 'Children');
fp.SaveFigs(handles,'-svg','Timecourse',savedir,1); close all












