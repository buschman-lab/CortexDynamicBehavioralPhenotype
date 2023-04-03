function MotifDecorder_VPA(data_dir,varargin)
%Evaluating the ability to use motif loadings to decode between VPA and Saline
opts.data_name = 'stats_refit';
opts.file_string = 'chunk';
opts.type = 'box';
opts.norm_method = 'zscore'; 
opts.balanced = 1; %boolean that will use random subset of the data to balance between groups. 
opts = ParseOptionalInputs(opts,varargin);

%load data
file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});
[mouse, ~] = MouseNumFromFileName(file_list);
grp = isVPA(mouse)'+1;
data_all = cellfun(@(x) load(x,opts.data_name),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.(opts.data_name),data_all,'UniformOutput',0);
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});

%normalize
switch opts.norm_method
    case 'zscore'
        data = zscore(data,[],1);
    case 'none'
    otherwise 
        error('Unknown normalization method');
end

%balance 
if opts.balanced
   temp = min(sum(grp==1),sum(grp==2));
   vpa = find(grp==2);
   vpa = vpa(randperm(numel(vpa),temp));
   sal = find(grp==1);
   sal = sal(randperm(numel(sal),temp));   
   data = data(cat(1,vpa,sal),:);
   grp = grp(cat(1,vpa,sal));
end
    
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

end %function end


