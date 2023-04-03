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
temp = load('Social_AllOriginalDataFileList.mat');
file_list = temp.file_list;

%save location 
save_dir ='Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

parameter_class = 'general_params_vpa';
gp = loadobj(feval(parameter_class)); %this is not needed here, but demonstrates how I load this class in other functions by just passing the string. 

basis_dir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution\ALL_sparse';
%% Post processing and fit old motifs to this data
job_id = cell(1,numel(file_list));
for cur_file = 1:numel(file_list)      
    %deconvolve and split the data
    script_name = WriteBashScript(sprintf('%d',cur_file),'VPA_FitMotifs_Social',{cur_file,...
        ConvertToBucketPath(save_dir),ConvertToBucketPath(basis_dir),parameter_class},...
        {'%d',"'%s'","'%s'","'%s'"},...
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


%%
ssh2_close(s_conn);
clear username password sconn







