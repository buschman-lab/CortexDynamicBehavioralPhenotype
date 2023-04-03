%%%% Processing Pipeline for Resting State Widefield Imaging Logs
%Camden MacDowell 2018 (camdenm@princeton.edu Princeton: Buschman Lab)

%TL;DR : User selects log files and save directory. Code saves timepoints
%and imaging frames to exclude due to animal movement. 

%SUMMARY: Takes log files (either individually or as batch) and
%outputs epochs of animal motion (e.g. raw timestamps and raw framestamps when an
%animal was moving +/- a user defined window). Also outputs quantifiable
%stats about the movement (e.g. total movement epochs). Combines this
%infromation with the recording opts files and saves as RecordingLogData in
%a user defined folder (or in the original directory of the image if no
%target directory). As a gut check this also checks if there were any
%dropped frames in the imaging sequence. This works for both individual
%animals resing state paired animals resting state (e.g. two piezos). 

%% Choose Logs and initialize opts
if nargin<2, logopts = setLogOpts; end

COUNT = 0; 
for cur_folder = 1:length(folder_list)
    cd(folder_list{cur_folder});
    try
        if logopts.verbose
            fprintf('\n\tNow processing log %d, out of %d...\n',cur_folder,length(folder_list))
        end
        %load recording options file (first find the one that matches date/time
        %of log file
        [path, base] = fileparts(folder_list{cur_folder});
        base = erase(base,'_Processing'); %Name of single resting state log file
        base = erase(base,'-RestingState-log'); %Name of single resting state log file
        base = erase(base,'-PairedSocial-log'); %Name of the paired resting state log file
        opts_fn = [path filesep sprintf('%s-OptsFile.mat',base)];
        if ~exist(opts_fn,'file')        
            error('No options file for file %s \n double check the file name',folder_list{cur_folder})
        else
            load(opts_fn)
        end
        %Number of log inputs to read (e.g. # analog inputs + 1 (the timing vector))
        numIn = length(opts.AnalogInputMaps)+1;
        %Logs are stored as binary data .mat files so use fopen/fread
        %row 1 is timing data
        [logdata, ~] = fread(fopen(folder_list{cur_folder},'r'),[numIn,inf],'double');
        %downsample
        logdata = downsample(logdata',logopts.dsFactor)';
        %Compute window (s) given samplign rate and downsampling
        window =(logopts.VarWindow*opts.AnalogInRate)/logopts.dsFactor; 

        %Animal loop to check motion if there are multiple animals
        MotionTrace = zeros(length(logopts.mChan),length(logdata(1,:)));
        for cur_animal = 1:length(logopts.mChan) %animal loop
            %Compute moving variance of signal (motion = High variance)
            mData = movvar(logdata(logopts.mChan(cur_animal)+1,:),window); 
            tData = logdata(1,:);
            %Mean and std of variance for motion detection
            threshold = median(mData)+logopts.threshold*std(mData);
            motion = zeros(1,length(mData));
            motion(mData>=threshold)=1;
            %combined motion epochs seperated by < user defined duration by 'closing' the 
            %motion logical (as you would an image to remove spurius noise). 
            %So you use a structuring element to remove all motion epochs shorter than
            %user defined duration
            try
                se = strel('square',(logopts.removestilldur*opts.AnalogInRate)/logopts.dsFactor);
                motion = imclose(motion,se);
                %remove motion epochs shorter than user defined duration
                se = strel('square',(logopts.motiondur*opts.AnalogInRate)/logopts.dsFactor);
                motion = imdilate(motion,se);
                se = strel('square',(logopts.excludewindow*opts.AnalogInRate)/logopts.dsFactor);
                motion = imdilate(motion,se);
            catch
                warning('memory error - you need to set logopts.dsFactor to at least 10')
            end
            MotionTrace(cur_animal,:) = motion; 
        end
        %Combined motion traces if mutliple animals
        MotionTrace = sum(MotionTrace,1);
        MotionTrace(MotionTrace>=1) = 1; %combined for when both animals are active at same time

        %Add summary stats to logopts
        logopts.MotionTrace = MotionTrace;
        logopts.numEpochs = (sum(diff(MotionTrace)~=0))/2;
        logopts.totalMovementDurationInSecs = (length(MotionTrace(MotionTrace==1))...
            *logopts.dsFactor)/opts.AnalogInRate;
        logopts.PercentTimeActive = (logopts.totalMovementDurationInSecs*opts.AnalogInRate)...
            /(logopts.dsFactor*length(MotionTrace))*100;
        
        %To combine everything while testing (comment out later)
        COUNT = COUNT+1;
        AllAnimalMovementData(COUNT).AnimalNum = str2num(base(1:4));
        AllAnimalMovementData(COUNT).Date = str2num(base(6:7));
        AllAnimalMovementData(COUNT).numEpochs = logopts.numEpochs;
        AllAnimalMovementData(COUNT).PercTime = logopts.PercentTimeActive;
        AllAnimalMovementData(COUNT).totalSecs = logopts.totalMovementDurationInSecs;
        
        %%DETERMINE WHICH FRAMES OCCUR DURING MOVEMENT 
        %Note for recs spring-fall 2018: Unfortunately the data log is not capturing every frame due to too low sampling frequency in spring 2018 experiments, so we
        %must figure it out it post-hoc. This should be improved in future.

        trigData = logdata(logopts.trigChan+1,:)>=1; %make logical of trigger channel
        startTime = tData(diff(trigData)~=0); %first trigger is start time (s)
        exposec = opts.exposuretime/1000; %exposure time in sec 

        %If >1 rec per session,  then there will be a save delay between recs, need to
        %include that to get accurate frame time points
        if opts.numRec > 1 %if more than one recording during session
            %preallocate
            finalframetime = NaN(1,opts.numRec);
            FrameList = NaN(opts.numFrames,opts.numRec);
            FrameTimes = NaN(opts.numFrames,opts.numRec);  
            for rec = 1:opts.numRec
                if rec == 1 %first recorinding in session
                    finalframetime(:,rec) = (exposec*opts.numFrames)+startTime(1)-exposec;
                    FrameList(:,rec) = (1:1:opts.numFrames); %frame number
                    FrameTimes(:,rec) = (startTime(1):exposec:finalframetime(rec));    
                else %all subsquent recordings in session
                    finalframetime(:,rec) =  (exposec*opts.numFrames)+finalframetime(rec-1)+opts.saveDelay/1000-exposec; 
                    tempstarttime = finalframetime(rec-1)+opts.saveDelay/1000;
                    FrameList(:,rec) = (FrameList(end,(rec-1))+1:1:opts.numFrames+FrameList(end,(rec-1))); %frame numbers starting at previous number
                    FrameTimes(:,rec) = (tempstarttime:exposec:finalframetime(rec));  %frame times starting at previous time
                end
            end
            FrameList = FrameList(:)';
            FrameTimes = FrameTimes(:)';
        else %if only one rec per session
            finalframetime = (exposec*opts.numFrames)+startTime(1)-exposec;
            FrameList = (1:1:opts.numFrames); %frame number
            FrameTimes = (startTime(1):exposec:f);    
        end
        FrameTimesComplete = FrameTimes; %non-rounded version for plotting below
        FrameTimes = round(FrameTimes,0); %round heavily to not miss any frames when checking if in the downsample motion trace 

        % Find the frames times that occur during motion
        tMotion = tData(MotionTrace==1); 
        indx = ismember(FrameTimes,tMotion);
        logopts.MotionFrameIndx = FrameList(indx==1);

        if logopts.verbose
            MotionTraceFig = figure('name','Piezo Trace','units','normalized','outerposition',[0 0 1 1]);    
            plot(tData,logdata(logopts.mChan+1,:),'k'); hold on; 
            title(sprintf('PiezoTrace for rec %d',cur_folder));
            plot(tData,ones(1,length(mData))*median(mData),'g');
            plot(tData,ones(1,length(mData))*threshold,'r');    
            plot(tData,mData,'c');
            plot(tData,MotionTrace,'y','linewidth',3);
            scatter(FrameTimesComplete(indx==1),indx(indx==1)*3,'.','r');
            legend('Raw','Median Variance','Threshold Variance','Variance','Motion Epoch','Imaging Frames During Motion')
            pause(1);
            logopts.MotionCaptureTrace = MotionTraceFig;
            close
        end

        RecordingOpts = opts;  %To keep opts variables straight moving forward
        save_fn = [SaveDir filesep sprintf('Mouse%s_Processing',base) filesep sprintf('%s_RecLogAndMovementInfo.mat',base)];
        if ~exist([SaveDir filesep sprintf('Mouse%s_Processing',base)],'dir')
            mkdir([SaveDir filesep sprintf('Mouse%s_Processing',base)]);
        end
        save(save_fn,'RecordingOpts','logopts'); 
        if logopts.verbose
            fprintf('\n\tDone processing log %d, out of %d...\n',cur_folder,length(folder_list))
        end
        
    catch
        warning('There was an issue processing movement log (probably because on live mode, not triggered) %s',folder_list{cur_folder})
    end
 end %Gather Rec opts loop
 
%  save('C:\Users\Camden\Desktop\AllAnimalMovementData.mat','AllAnimalMovementData', '-v7.3');















