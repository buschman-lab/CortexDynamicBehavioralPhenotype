function VPA_FitMotifs(cur_file,save_dir)

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

%Load Preprocessed Data
temp = load('AllOriginalDataFileList.mat');
if ispc
   file_list = temp.file_list;
else
   file_list = temp.bucket_list;
end

%load parameters
gp = general_params_vpa;

%Load the data
temp = load(file_list{cur_file});
%spatially bin to 68x68
stack = SpatialBin(temp.dff,2,[],1);
stack(repmat(nanvar(stack,[],3)<=eps,[1,1,size(stack,3)])==1)=NaN;

%skip if too short (<12min) recording
if size(stack,3) < gp.fps*11*60; return; end

%mask with conservative mask
temp = load('FigureMask.mat');
mask = repmat(temp.mask,[1,1,size(stack,3)]);
stack(mask==0)=NaN;

%threshold and split
[~, nanpxs, data_train_all, data_test_all] = ProcessAndSplitData(stack,[],gp);

%Loop through the chunks and fit/Xval motifs; 
for i = 1:size(data_train_all,2)
    fprintf('\nWorking on chunk %d\n',i);
    if gp.reverse_fit %flip onto testing data
        [W, H, stats_train, stats_test] = FitandValidateMotifs(squeeze(data_test_all(:,:,i)),squeeze(data_train_all(:,:,i)),gp,0);            
    else
        [W, H, stats_train, stats_test] = FitandValidateMotifs(squeeze(data_train_all(:,:,i)),squeeze(data_test_all(:,:,i)),gp,0);            
    end  
    [~,fn] = fileparts(file_list{cur_file});
    data_train = squeeze(data_train_all(:,:,i));
    data_test = squeeze(data_test_all(:,:,i));
    if gp.reverse_fit
        save([save_dir, filesep, sprintf('%s_reverse_chunk%d.mat',fn,i)],'W','H','stats_train','stats_test','nanpxs','data_train','data_test','-v7.3')
    else
        save([save_dir, filesep, sprintf('%s_chunk%d.mat',fn,i)],'W','H','stats_train','stats_test','nanpxs','data_train','data_test','-v7.3')
    end    
end

fprintf('Done');

end
    







