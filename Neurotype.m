function [Y,Y_null, betas, betas_null,allstats] = Neurotype(data_dir,varargin)
%Classifier-based dimensionality reduction shift 

opts.data_name = 'stats_refit';
opts.file_string = 'chunk';
opts.norm_method = 'zscore'; 
opts.avg_method = 'mean';
opts.num_classifiers = 5;
opts.neuro_col = 1:16;
opts.holdout = 0.4;
opts.verbose = 0; %boolean, 1 if make figures. 

opts = ParseOptionalInputs(opts,varargin);

fp = fig_params_vpa; 

%load data
file_list = GrabFiles(['\w*',opts.file_string,'\w*'],0,{data_dir});

[mouse, ~] = MouseNumFromFileName(file_list);
data_all = cellfun(@(x) load(x,opts.data_name),file_list,'UniformOutput',0);
data_all = cellfun(@(x) x.(opts.data_name),data_all,'UniformOutput',0);
data = cellfun(@(x) x.loadings{1}, data_all,'UniformOutput',0);
data = cat(1,data{:});
data = data(:,opts.neuro_col);

%average by mouse
switch opts.avg_method
    case 'mean'        
        data_mouse = arrayfun(@(x) nanmean(data(mouse==x,:)), unique(mouse),'UniformOutput',0);
    case 'median'
        data_mouse = arrayfun(@(x) nanmedian(data(mouse==x,:)), unique(mouse),'UniformOutput',0);        
    otherwise 
        error('Unknown averaging method'); 
end
data_mouse = cat(1,data_mouse{:});

%normalize
switch opts.norm_method
    case 'zscore'
        data_mouse = zscore(data_mouse,[],1);
    case 'none'
    otherwise 
        error('Unknown normalization method');
end

grp = isVPA(unique(mouse))'+1;        

[Y, Y_null, betas, betas_null, allstats] = DistanceToHyperplane(data_mouse,grp,opts.num_classifiers,opts.verbose,opts.holdout);

end %function end


