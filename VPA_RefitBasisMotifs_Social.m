function [w,H,stats_refit,bad_pxl] = VPA_RefitBasisMotifs_Social(data_train,nanpxs, basis_dir,parameter_class)

%Add paths
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis'));
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis'));
end

gp = loadobj(feval(parameter_class));

%load the basis motifs
temp = load([basis_dir filesep 'basismotifs.mat'],'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W = temp.W_basis;
W(:,temp.noise_clusters,:)=[];

%recondition, smooth, and flatten
if numel(gp.smt_kernel)>2 %autodetermine appropriate smoothing value
    fprintf('\n\t Autofitting Smoothing Value');
    gp.smt_kernel = AutoFitSmoothingLevel(data_train,nanpxs,W,gp);
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));    
elseif numel(gp.smt_kernel)==2
    fprintf('\n\t Set Smoothing Value');
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
else
    fprintf('\n\t No Smoothing');
    data_train = reshape(conditionDffMat(data_train',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
end

%Only work on shared pixels
[~, bad_pxl] = SharedPixels(W,data_train);
data_train(bad_pxl,:) = [];
W(bad_pxl,:,:) = []; 
    
%optional add of refitting the H sparsity term
[w,H] = fpCNMF(data_train,'non_penalized_iter',...
    gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
    'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
    'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
stats_refit = CNMF_Stats(w,H,data_train,0);
stats_refit.smoothingkernel = gp.smt_kernel;


end
