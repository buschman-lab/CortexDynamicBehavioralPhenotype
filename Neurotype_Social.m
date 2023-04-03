function [Y,Y_null, betas, betas_null,allstats,mice] = Neurotype_Social(varargin)
%Classifier-based dimensionality reduction shift 

opts.data_name = 'stats_refit';
opts.file_string = 'refit';
opts.norm_method = 'zscore'; 
opts.avg_method = 'mean';
opts.num_classifiers = 5;
opts.neuro_col = 1:16;
opts.holdout = 0.4;
opts.verbose = 0; %boolean, 1 if make figures. 

opts = ParseOptionalInputs(opts,varargin);

fp = fig_params_vpa; 

%load data
file_list = GrabFiles(['\w*','refit','\w*'],0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social'});
%only use files where same groups paired
data_social = cellfun(@(x) load(x,'loadings','mouseID','stats_refit'),file_list,'UniformOutput',0);
mouse = cellfun(@(x) x.mouseID, data_social,'UniformOutput',0);
mouse = [mouse{:}];
%average per mouse
switch opts.avg_method
    case 'mean' 
        data = cellfun(@(x)  nanmean(x.loadings,2)*100,data_social,'UniformOutput',0);
        data = cat(2,data{:});
        data_mouse = arrayfun(@(x) nanmean(data(:,mouse==x),2), unique(mouse), 'UniformOutput', 0);
    case 'median'
        data = cellfun(@(x)  nanmedian(x.loadings,2)*100,data_social,'UniformOutput',0);
        data = cat(2,data{:});
        data_mouse = arrayfun(@(x) nanmedian(data(:,mouse==x),2), unique(mouse), 'UniformOutput', 0);   
    otherwise 
        error('Unknown averaging method'); 
end
%combine animals 
data_mouse = cat(2,data_mouse{:})'; %(data mouse should be mice x motifs)
data_mouse = data_mouse(:,opts.neuro_col);

%normalize
switch opts.norm_method
    case 'zscore'
        data_mouse = zscore(data_mouse,[],1);
    case 'none'
    otherwise 
        error('Unknown normalization method');
end

grp = isVPA(unique(mouse))'+1;      

mice = unique(mouse);

[Y, Y_null, betas, betas_null, allstats] = DistanceToHyperplane(data_mouse,grp,opts.num_classifiers,opts.verbose,opts.holdout);

end %function end


