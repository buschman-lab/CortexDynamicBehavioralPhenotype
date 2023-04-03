function VPA_RefitBasisMotifs(fn, basis_dir, dynamics)

if nargin <3; dynamics = 1; end

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

gp = general_params_vpa;

%load the basis motifs
temp = load([basis_dir filesep 'basismotifs.mat'],'W_basis','noise_clusters');
%ignore noise clusters when refitting (you could also refit these and then remove them and their representation in the data
W = temp.W_basis;
W(:,temp.noise_clusters,:)=[];

%optionally remove the dynamics
if dynamics == 0 %flatten and 
   W = repmat(nanmean(W,3),[1,1,size(W,3)]);
elseif dynamics == 2 %reverse
   for i = 1:size(W,2)   
       W(:,i,:) = fliplr(squeeze(W(:,i,:)));
   end
elseif dynamics ==3 %single frame
    W = nanmean(W,3);
end

%load the data_train and the data_test data
fprintf('\n\tLoadings Data\n'); 
temp = load(fn,'data_test','data_train','nanpxs','stats_train');
data_train = temp.data_train;
data_test = temp.data_test;
nanpxs = temp.nanpxs; 
[~,name] = fileparts(fn);
if dynamics ~= 1
    name = [name,sprintf('nodynamics%d',dynamics)];
end

%recondition, smooth, and flatten
if numel(gp.smt_kernel)>2 %autodetermine appropriate smoothing value
    fprintf('\n\t Autofitting Smoothing Value');
    gp.smt_kernel = AutoFitSmoothingLevel(cat(2,data_train,data_test),nanpxs,W,gp);
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
elseif numel(gp.smt_kernel)==2
    fprintf('\n\t Set Smoothing Value');
    data_train = reshape(SpatialGaussian(conditionDffMat(data_train',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(SpatialGaussian(conditionDffMat(data_test',nanpxs),gp.smt_kernel),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
else
    fprintf('\n\t No Smoothing');
    data_train = reshape(conditionDffMat(data_train',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_train,2));
    data_test = reshape(conditionDffMat(data_test',nanpxs),gp.pixel_dim(1)*gp.pixel_dim(2),size(data_test,2));
end

%Only work on shared pixels
[~, bad_pxl] = SharedPixels(W,cat(2,data_test,data_train));
data_train(bad_pxl,:) = [];
data_test(bad_pxl,:) = [];
W(bad_pxl,:,:) = []; 
    
%optional add of refitting the H sparsity term
[w,H] = fpCNMF(data_train,'non_penalized_iter',...
    gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
    'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
    'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
stats_refit = CNMF_Stats(w,H,data_train,0);
stats_refit.smoothingkernel = gp.smt_kernel;

%save off
save([basis_dir filesep name 'train.mat'],'w','H','stats_refit','bad_pxl');

[w,H] = fpCNMF(data_test,'non_penalized_iter',...
    gp.non_penalized_iter,'penalized_iter',gp.penalized_iter_refit,...
    'speed','fast','verbose',0,'lambda',0,'w_update_iter',0,...
    'ortho_H',gp.ortho_H,'W',W,'sparse_H',0);
stats_refit = CNMF_Stats(w,H,data_test,0);

save([basis_dir filesep name 'test.mat'],'w','H','stats_refit','bad_pxl');


%save off
