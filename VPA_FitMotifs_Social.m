function VPA_FitMotifs_Social(cur_file,save_dir,basis_dir,parameter_class)

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
temp = load('Social_AllOriginalDataFileList.mat');
if ispc
   file_list = temp.file_list;
else
   file_list = temp.bucket_list;
end

%remove the 'Final'
file_list = cellfun(@(x) erase(x,'_Final'), file_list,'UniformOutput',0);
[~,fn] = fileparts(file_list{cur_file});
%load parameters
gp = loadobj(feval(parameter_class)); 

%Load the data
temp = load(file_list{cur_file});
%spatially bin to 68x68
stack = SpatialBin(temp.dff,2,[],1);
stack(repmat(nanvar(stack,[],3)<=eps,[1,1,size(stack,3)])==1)=NaN;

%skip if too short (<12min) recording
if size(stack,3) < gp.fps*11*60; return; end

%get the mouse ID
if ~isempty(regexp(fn,'mouse1_','ONCE'))
    mouseID = str2double(fn(13:16));
elseif ~isempty(regexp(fn,'mouse2_','ONCE'))
    mouseID = str2double(fn(18:21));
else
    return
end   

%mask with conservative mask
temp = load('FigureMask.mat');
mask = repmat(temp.mask,[1,1,size(stack,3)]);
stack(mask==0)=NaN;

%threshold and split
[~, nanpxs, data_train_all, data_test_all] = ProcessAndSplitData(stack,[],parameter_class);

%combine 
data = NaN(size(data_test_all,1),size(data_test_all,2),size(data_test_all,3)*2);
data(:,:,1:2:end)= data_train_all;
data(:,:,2:2:end)= data_test_all;

%Loop through the chunks and refit basis motifs
w = cell(1,size(data,3));
H = cell(1,size(data,3));
stats_refit = cell(1,size(data,3));
bad_pxl = cell(1,size(data,3));
X = cell(1,size(data,3));
for i = 1:size(data,3)
    [w{i},H{i},stats_refit{i},bad_pxl{i}] = VPA_RefitBasisMotifs_Social(data(:,:,i), nanpxs, basis_dir,parameter_class);  
    X{i} = data(:,:,i);
end

%combine and get sliding window loadings
X = cat(2,X{:});

N = size(H{1},1); %number of basis motifs
%loop through each motif and compute it's pev sliding window style
dur = 13*2*60;
shift = 13*30;
pev = cell(1,N);
for j = 1:N
  Xhat = cell(1,6);  
  for cur_chunk = 1:6
      Xhat{cur_chunk} = tensor_convolve(w{cur_chunk}(:,j,:),H{cur_chunk}(j,:));
  end
  Xhat = cat(2,Xhat{:});
  pev{j} = SlidingWindowPEV(X,Xhat,dur,shift);
end
pev = cat(1,pev{:});
loadings = pev./sum(pev,1);%save off loadings

save([save_dir, filesep, sprintf('%s_refit.mat',fn)],'w','H','stats_refit','nanpxs','loadings','data','bad_pxl','dur','shift','mouseID','-v7.3')

fprintf('Done');

end
    







