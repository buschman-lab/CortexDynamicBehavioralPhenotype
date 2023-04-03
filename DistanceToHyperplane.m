function [Y, Y_null, betas_avg, betas_null, allstats,bias_avg] = DistanceToHyperplane(data,grp,num_classifiers,verbose,holdout, resetshf)
%Camden MacDowell - timeless
%computes a real and null hyperplane using features in data for 2 groups in
%grp 
fp = fig_params_vpa; 

if nargin <6; resetshf=1; end
%loop through multiple classifiers
if resetshf
    rng('default');
end
auc_train = zeros(1,num_classifiers); %for gutcheck. e.g. this should be approx >0.9 
auc = zeros(1,num_classifiers); 
auc_null = zeros(1,num_classifiers);  %null = classifier trained on shuffled data
auc_shuf = zeros(1,num_classifiers); %shuf = reg class evaluated on shuffled data
betas = zeros(size(data,2),num_classifiers);
bias = zeros(1,num_classifiers);
betas_null = zeros(size(data,2),num_classifiers);
bias_null = zeros(1,num_classifiers);
allstats = cell(1,num_classifiers);
for i = 1:num_classifiers
    if num_classifiers ==1 %train on the whole thing
        partition = [];
        partition.training = true(numel(grp),1); %train on entire dataset
        partition.test = true(numel(grp),1); %doesn't matter as you aren't testing this classifier
        [~, observed, shuffled, trainauc] = SVMClassifier_Binary([data,grp],partition,'optimize',0,'kernel','linear','featureselect','none','nshuf',1,'verbose',0); 
    else    
        [~, observed, shuffled, trainauc] = SVMClassifier_Binary([data,grp],[],'holdout',holdout,'optimize',0,'kernel','linear','featureselect','none','nshuf',1,'verbose',0); 
    end
    %Train a totally shuffled classifier (this is different then applying the other classifier to shuffled data because it allows us to also capture null biases and betas
    [~, null, ~, ~] = SVMClassifier_Binary([data,grp(randperm(numel(grp),numel(grp)))],[],'holdout',holdout,'optimize',0,'kernel','linear','featureselect','none','nshuf',0,'verbose',0); 
    auc_train(i) = nanmean([trainauc(:).AUC]);
    auc(i) = observed.AUC;
    auc_shuf(i) = shuffled.AUC;   
    auc_null(i) = null.AUC;
    betas(:,i) = observed.Classifier.Beta;
    betas_null(:,i) = null.Classifier.Beta;
    bias(i) = observed.Classifier.Bias;
    bias_null(i) = null.Classifier.Bias;
    allstats{i} = observed;
end
    

%generate an average classifier and average null classifier
bias_avg = mean(bias);
bias_null = mean(bias_null);
betas_avg = mean(betas,2);
betas_null = mean(betas_null,2);

%compute distance to the hyperplane
DistToHyperplane= @(x) arrayfun(@(n) dot(x(n,:),betas_avg)+bias_avg,...
    (1:1:size(x,1)),'UniformOutput',1);
Y = DistToHyperplane(data);

%get the AUC of the average classifier on the whole dataset
validationResponse = Y; 
validationResponse(validationResponse<0)=0;
validationResponse(validationResponse>0)=1;
validationResponse = validationResponse+1;
[~,~,~,AUC_general] = perfcurve(grp,validationResponse(:),2);
allstats{end+1} = AUC_general;

%distance using null classifier
DistToHyperplane= @(x) arrayfun(@(n) dot(x(n,:),betas_null)+bias_null,...
    (1:1:size(x,1)),'UniformOutput',1);
Y_null = DistToHyperplane(data);

%compare the AUC of observed, shuffled, and training classifiers
if verbose
    figure; hold on; 
    plot(rand(num_classifiers,1)/4,auc_train,'.','color',[0.5 0.5 0.5],'markersize',15)
    plot(rand(num_classifiers,1)/4+1,auc,'.','color',[0.5 0.5 0.5],'markersize',15)
    plot(rand(num_classifiers,1)/4+2,auc_shuf,'.','color',[0.5 0.5 0.5],'markersize',15)
    plot(rand(num_classifiers,1)/4+3,auc_null,'.','color',[0.5 0.5 0.5],'markersize',15)
    ylim([0 1]); 
    [~,p] = ttest2(auc,auc_shuf,'Tail','right'); line([1,2.5],[0.925,0.925],'linewidth',2,'color','r'); text(1.5,.9,sprintf('p=%0.2g',p),'color','red'); %significance testing
    [~,p] = ttest2(auc,auc_null,'Tail','right'); line([1,3.5],[0.975,0.975],'linewidth',2,'color','b'); text(2.5,.95,sprintf('p=%0.2g',p),'color','blue'); %significance testing
    set(gca,'xtick',[0:3],'xticklabel',{'Train','Real','Shuf','Null'})
    xlabel('classifier');
    ylabel('AUC');
    fp.SetTitle(gca,'Classifier > Shuffled/Null '); 
    fp.FormatAxes(gca);

    %check if any correlation between the null and real 
    figure; hold on; 
    plot(Y,Y_null,'.','markersize',15,'color','k'); axis equal;
    xlabel('Real Neurotype Axis');
    ylabel('Null Neurotype Axis');
    [rho,p] = corr(Y',Y_null');
    fp.SetTitle(gca,sprintf('rho=%0.2g p=%0.2g',rho,p));
end %verbosity 

end