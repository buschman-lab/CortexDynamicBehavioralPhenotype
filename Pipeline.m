%VPA_Pipeline
% Open ssh connection
username = input(' Spock Username: ', 's');
password = passcode();
s_conn = ssh2_config('spock.princeton.edu',username,password);

%Add paths
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\Widefield_Imaging_Analysis'));
addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\VPA_Mesoscale_Analysis'));

%Load data location
temp = load('AllOriginalDataFileList.mat');
file_list = temp.file_list;

%save location 
save_dir ='Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end


%% Post processing and fit motifs
job_id = cell(1,numel(file_list));
for cur_file = 1:numel(file_list)      
    %deconvolve and split the data
    script_name = WriteBashScript(sprintf('%d',cur_file),'VPA_FitMotifs',{cur_file,...
        ConvertToBucketPath(save_dir)},...
        {'%d',"'%s'"},...
        'sbatch_time',600,'sbatch_memory',12,...
        'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis/");
    %Run job
    response = ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory of the bash script
        sprintf('sbatch %s',script_name)]);         
    %get job id
    if cur_file~=numel(file_list)
        job_id{cur_file} = [erase(response.command_result{1},'Submitted batch job '),','];
    else
        job_id{cur_file} = erase(response.command_result{1},'Submitted batch job ');
    end
end


%% Cluster motifs as a whole, as VPA, and as SAL
% basis_save_dir_base ='Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_nosmooth_90\';
basis_save_dir_base ='Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\';
% group = [{1},{0},{2}];
group = [{2}];
dynamics = 1; %0 = flattened motifs, 1 = normal, 2 = reversed motifs, 3=single frame average
job = [];
for i = 1:numel(group)   
    %split by the fit type
    cur_group = group{i};
    if cur_group == 1
        save_fn = 'VPA';
    elseif cur_group ==0
        save_fn = 'SAL'; 
    else
        save_fn = 'ALL_sparse';
    end        
    basis_save_dir = [basis_save_dir_base save_fn];
    if ~exist(basis_save_dir,'dir')
        mkdir(basis_save_dir);
    end
%     
%     %Cluster Motifs
%     script_name = WriteBashScript(sprintf('%d',i),'VPA_ClusterMotifs',{ConvertToBucketPath(save_dir),...
%         ConvertToBucketPath(basis_save_dir),cur_group},...
%         {"'%s'","'%s'",'%d'},...
%         'sbatch_time',1440,'sbatch_memory',32,...
%         'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis/");
%     
%     %Run job
%     if exist('job_id','var')
%         response = ssh2_command(s_conn,...
%             ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory of the bash script
%             sprintf('sbatch --dependency=afterany:%s %s',[job_id{:}],script_name)]);  
%     else
%         response = ssh2_command(s_conn,...
%             ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory of the bash script
%             sprintf('sbatch %s',script_name)]);          
%     end
%     
%     %get job id
%     job = erase(response.command_result{1},'Submitted batch job ');
    
%     Refit basis motifs to the entire data set
    VPA_RefitBasisMotifs_Swarm(save_dir,basis_save_dir,job,s_conn,dynamics)
    
end

ssh2_close(s_conn);
clear username password sconn







