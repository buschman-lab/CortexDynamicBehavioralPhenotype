function [Y,Y_null,Y_full,Y_null_full,grp,Imaged_only,Imaged_only_null,betas,betas_imaged,allstats,bias_avg,bias_avg_imaged] = Phenotype(data_dir,varargin)
%Classifier-based dimensionality reduction shift 
%output is either just the imaged mice (Y/Y_null) or all mice (full suffix)

opts.file_string = 'chunk';
opts.num_classifiers = 25;
opts.verbose = 1; %boolean, 1 if make figures.
opts.behavior_col = 2:13;  %list of behavioral columns to use, default = [2:13]
opts.holdout = 0.5;
opts.reversegroups = 0;
opts = ParseOptionalInputs(opts,varargin);

%load imaging mouse numbers
file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});
% file_list = substring_filename(file_list,'nodynamics',opts.nodynamics);
[imaged_mice, ~] = unique(MouseNumFromFileName(file_list));

%load data
data = GetBehavioralData();
mouse = data(:,1);
grp = data(:,end)+1;
if opts.reversegroups
    temp=grp;
    grp(temp==1)=2;
    grp(temp==2)=1;
end
    
data = data(:,opts.behavior_col);

%compute distance
[Y_full, Y_null_full, betas,~,allstats,bias_avg] = DistanceToHyperplane(data,grp,opts.num_classifiers,opts.verbose,opts.holdout);

%loop through each so that you make sure it's in the order of imaged_mice
Y = zeros(1,numel(imaged_mice));
Y_null = zeros(1,numel(imaged_mice));
idx = zeros(1,numel(imaged_mice));
for i = 1:numel(Y_null)   
   Y(i) = Y_full(mouse==imaged_mice(i));
   Y_null(i) = Y_null_full(mouse==imaged_mice(i)); 
   idx(i) = find(mouse==imaged_mice(i)==1); %for plotting below
end

%only using imaged mice
[~, data] = GetBehavioralData();
imaged_mouse_data = zeros(numel(imaged_mice),size(data,2));
imaged_grp = zeros(numel(imaged_mice),1);
for i = 1:numel(imaged_mice)
    imaged_mouse_data(i,:) = data(mouse==imaged_mice(i),:);
    imaged_grp(i,:) = grp(mouse==imaged_mice(i));
end
%zscore
imaged_mouse_data = zscore(imaged_mouse_data,0,1);
[Imaged_only, Imaged_only_null, betas_imaged,bias_avg_imaged] = DistanceToHyperplane(imaged_mouse_data,imaged_grp,opts.num_classifiers,opts.verbose,opts.holdout);


if opts.verbose 
   fp = fig_params_vpa;
   figure('position',[432   480   797   227]); hold on;   
   temp = rand(numel(Y_full),1);
   plot(Y_full(grp==1),temp(grp==1),'.','markersize',15,'color',fp.c_sal);    
   plot(Y_full(grp==2),temp(grp==2),'.','markersize',15,'color',fp.c_vpa);    
   plot(Y,temp(idx),'o','markersize',10,'color','k');    
   set(gca,'ytick',[],'xtick',get(gca,'xlim'),'xticklabel',{'SAL-like','VPA-like'});
   xlabel('Phenotype')
   legend('SAL','VPA','Imaged-Mice','Location','eastoutside')
   fp.FormatAxes(gca)
end

end %function end


