function VPA_RefitBasisMotifs_Swarm(motif_dir,basis_dir,job_id,s_conn,dynamics)
%Camden MacDowell

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

file_list = GrabFiles('\w*chunk\w*.mat', 0, {motif_dir});


%generate swarm
for cur_file = 1:numel(file_list) 
    script_name = WriteBashScript(sprintf('%d',1),'VPA_RefitBasisMotifs',{ConvertToBucketPath(file_list{cur_file}),ConvertToBucketPath(basis_dir),dynamics},...
    {"'%s'","'%s'",'%d'},...
    'sbatch_time',59,'sbatch_memory',10,...
    'sbatch_path',"/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA_Mesoscale_Analysis/");
    
    if ~isempty(job_id)
        ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch --dependency=afterok:%s %s',job_id,script_name)]); 
    else
        ssh2_command(s_conn,...
        ['cd /jukebox/buschman/Rodent\ Data/Wide\ Field\ Microscopy/Widefield_Imaging_Analysis/Spock/DynamicScripts/ ;',... %cd to directory
        sprintf('sbatch %s',script_name)]);  
    end
end

