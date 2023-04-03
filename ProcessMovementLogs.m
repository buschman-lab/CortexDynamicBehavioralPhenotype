function ProcessMovementLogs(fn)
%Camden MacDowell - timeless
%fn is path to the opts file for each recording
%see "ProcessLogs_old.m" for reference
if nargin <1
    fn = GrabFiles('\w*OptsFile.mat');
end

%loop through each opts file
for cur_file = 1:numel(fn)
    try
    opts = load(fn{cur_file},'opts'); opts = opts.opts;
    [fn_path, fn_base] = fileparts(fn{cur_file});
    
    log_fn = erase(fn_base,'-OptsFile');
    log_fn = [fn_path filesep log_fn '-RestingState-log.mat'];
    %load the data logger. Number of log inputs to read (e.g. # analog inputs + 1 (the timing vector))
    numIn = length(opts.AnalogInputMaps)+1;
    %Logs are stored as binary data .mat files so use fopen/fread
    %row 1 is timing data
    [logdata, ~] = fread(fopen(log_fn,'r'),[numIn,inf],'double');    
    
    %match to frames
    camsig = logdata([0 strcmp('ai2 - Frame readout',opts.AnalogInputMaps)]==1,:); %(preceeds exposure)
    motsig = abs(logdata([0 strcmp('ai0 - Piezo',opts.AnalogInputMaps)]==1,:));
    
    %average change in motor signal during each camera frame
    camsig(camsig<1)=0; camsig(camsig>=1)=1;
    idx = find(diff(camsig)==1);
    mot_ds = NaN(1,numel(idx));
    for i = 1:numel(idx) 
       if i == numel(idx)
          temp_idx = [idx(i),idx(i)+opts.exposuretime];
       else
          temp_idx = [idx(i),idx(i+1)];
       end
       mot_ds(i) = nanmean(abs(diff(motsig(temp_idx(1):temp_idx(2)))));
    end
    
    %threshold at 1std and combine short bouts and remove tiny bouts
    mot_thresh = mot_ds>std(mot_ds);    
    se = strel('square',50);
    mot_thresh = imclose(mot_thresh,se);
    %remove motion epochs shorter than user defined duration
    se = strel('square',20);
    mot_thresh = imdilate(mot_thresh,se);    
    
    %figure
    x = (1:numel(mot_thresh))*opts.exposuretime/1000/60;
    figure; hold on; 
    plot(x,mot_ds,'color',[0.6 0.6 0.6 0.6],'linewidth',1);
    ylabel('raw activity');
    yyaxis right    
    plot(x,mot_thresh,'color','k','linewidth',1.5);
    set(gca,'ycolor','k')
    set(gca,'ytick',[0 1],'YTickLabel',{'Low','High'});
    ylabel('arousal');
    xlabel('time (min)');
    
    
    log_fn = erase(fn_base,'-OptsFile');
    mouse = str2num(log_fn(1:4));
    
    %save off the data     
    savedir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Processed_MovementLogs';
    save([savedir filesep erase(fn_base,'-OptsFile') '.mat'],'mot_thresh','opts','mouse','mot_ds');    
    saveCurFigs(gcf,{'-dpng','-dsvg'},erase(fn_base,'-OptsFile'),savedir,0); close all
    catch
    end
end











