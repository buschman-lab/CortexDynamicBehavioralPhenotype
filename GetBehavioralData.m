function [zscoreData, rawData, testName] = GetBehavioralData()
%Compile data for each animal in each test:
warning off
% Load data
if ispc
    %VR change path 
%     vrpath = '\\cup\buschman\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\';
%     load(append(vrpath,'SummaryUSV.mat'));
%     load(append(vrpath,'SummaryMSP.mat'));
%     load(append(vrpath,'SummaryMarble.mat'));
%     load(append(vrpath,'SummaryGroom.mat')); 
%     load(append(vrpath,'SACombined.mat'));
%     load(append(vrpath,'SNCombined.mat'));
%     load(append(vrpath,'SummarySRR.mat'));
%     load(append(vrpath,'SummaryEYE.mat'));
%     load(append(vrpath,'SummaryW.mat')); %change from SummaryryW?
%     load(append(vrpath,'HABCombined.mat'));
    
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryUSV.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryMSP.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryMarble.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryGroom.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SACombined.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SNCombined.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummarySRR.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryEYE.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\SummaryW.mat');
    load('Z:\Users\Camden\Projects\VPA Model\Behavioral and Developmental Assessments\Analysis Code\Classifier Analysis\HABCombined.mat');
else
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryUSV.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryMSP.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryGroom.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryMarble.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SACombined.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SNCombined.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummarySRR.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryEYE.mat');
    load('/Volumes/buschman/Users/Camden/Projects/VPA Model/Behavioral and Developmental Assessments/Analysis Code/NeuroDevelopmental Analsis/SummaryData_8-12-2018/SummaryW.mat');
end
%%
clear zscoreData
% Preallocate matrices
% VPA (2)  or SAL (1)  group
mice = [SummaryUSV(1).Mouse'; SummaryUSV(2).Mouse'];   
tot_test = 10;
zscoreData = NaN(length(mice),tot_test); 
rawData= NaN(length(mice),tot_test); 
testName = cell(1,tot_test);
opts.randshuffle  = 0; 
opts.removenan = 0;
% Add mouse numbers
for cur_mouse = 1:length(mice)
    zscoreData(cur_mouse,1) = mice(cur_mouse); %mouse num is first column
    rawData(cur_mouse,1) = mice(cur_mouse);
end
%%%%%%%%%%%%%%%%%%%SOCIAL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adult Social Approach # is in Column 6;
SAIndex = NaN(1,length(SACombined));
SA_mice = NaN(1,length(SACombined));
total = 1; 
testName{1} = 'Social Approach';
if total
    for m = 1:length(SACombined)
        SAIndex(m) = length(SACombined(m).SAComp)/(length(SACombined(m).SAComp)+length(SACombined(m).EmptyComp));
        %The following subtracts the side preference from habituation from each animal 
        SAIndex1 = length(SACombined(m).SAComp)./(length(SACombined(m).EmptyComp)+length(SACombined(cur_mouse).SAComp));    
        SideIndx = length(HABCombined(m).SAComp)./(length(HABCombined(m).EmptyComp)+length(HABCombined(cur_mouse).SAComp));
        SAIndex(m)= SAIndex1-(SideIndx-0.5);
        SA_mice(m) = SACombined(m).MouseNum; 
        
    end
    for cur_mouse = 1:length(SA_mice) 
        indx = find(zscoreData==SA_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
        %calc zscoreData
        zscoreData(indx,2) = (SAIndex(cur_mouse)-nanmean(SAIndex))/nanstd(SAIndex);
        rawData(indx,2) = SAIndex(cur_mouse);
    end
else
    secs = 120; %calc index during first x seconds of the sa test.
    for m = 1:length(SACombined)
        SA_mice(m) = SACombined(m).MouseNum; 
        dur_social(m) = length(find(SACombined(m).SAComp < secs & SACombined(m).SAComp >= 1));
        dur_empty(m) = length(find(SACombined(m).EmptyComp < secs & SACombined(m).EmptyComp >=1));
        SAIndex(m)= dur_social/(dur_social+dur_empty);
    end
    for cur_mouse = 1:length(SA_mice) 
        indx = find(zscoreData==SA_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
        %calc zscoreData
        zscoreData(indx,2) = (SAIndex(cur_mouse)-nanmean(SAIndex))/nanstd(SAIndex);
        rawData(indx,2) = SAIndex(cur_mouse);
    end
end


% Adult Social Novelty # is in Column 7;
testName{2} = 'Social Novelty';
SNIndex = NaN(1,length(SNCombined));
SN_mice = NaN(1,length(SNCombined));
for m = 1:length(SNCombined)
    SNIndex(m) = length(SNCombined(m).StrComp)/(length(SNCombined(m).SAComp)+length(SNCombined(m).StrComp));
    SN_mice(m) = SNCombined(m).MouseNum; 
end
for cur_mouse = 1:length(SN_mice) 
    indx = find(zscoreData==SN_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,3) = (SNIndex(cur_mouse)-nanmean(SNIndex))/nanstd(SNIndex);
    rawData(indx,3) = SNIndex(cur_mouse);
end

% MaternalPreference Index for 1 min test is in Column 3; 
testName{3} = 'Maternal Preference';
MSP = [SummaryMSP(1).Indx1; SummaryMSP(2).Indx1];  %This has same order as the USV_mice numbers fyi
MSP_mice = [SummaryMSP(1).MouseNum; SummaryMSP(2).MouseNum];

for cur_mouse = 1:length(MSP_mice) 
    indx = find(zscoreData==MSP_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,4) = (MSP(cur_mouse)-nanmean(MSP))/nanstd(MSP);
    rawData(indx,4) = MSP(cur_mouse);
end

%%%%%%%%%%%%%%%%%%%%%%Development Tests%%%%%%%%%%%%%%%%%%%%%%%%%%

% USVs Are in Column 2; 
%(e.g. lower z score = 'more autistic phenotype')
testName{4} = 'USVs';
USVs = [SummaryUSV(1).CallRatePerMinPND5'; SummaryUSV(2).CallRatePerMinPND5'];  %This has same order as the USV_mice numbers fyi
USV_mice = [SummaryUSV(1).Mouse'; SummaryUSV(2).Mouse'];
for cur_mouse = 1:length(USV_mice)
    indx = find(zscoreData==USV_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,5) = (USVs(cur_mouse)-nanmean(USVs))/nanstd(USVs); 
    rawData(indx,5) = USVs(cur_mouse);
end


% EYE is in Column 9
testName{5} = 'Eye Opening';
EYE = [SummaryEYE(1).FirstDay; SummaryEYE(2).FirstDay];  %This has same order as the USV_mice numbers fyi
EYE_mice = [SummaryEYE(1).MouseNum; SummaryEYE(2).MouseNum];
for cur_mouse = 1:length(EYE_mice)
    indx = find(zscoreData==EYE_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,6) = (EYE(cur_mouse)-nanmean(EYE))/nanstd(EYE);
    rawData(indx,6) = EYE(cur_mouse);
end

% Weaning weight is in Column 10
testName{6} = 'Weaning Weight';
Weight = [SummaryW(1).WeanWeight; SummaryW(2).WeanWeight];  %This has same order as the USV_mice numbers fyi
Weight_mice = [SummaryW(1).Mouse; SummaryW(2).Mouse];
for cur_mouse = 1:length(Weight_mice)
    indx = find(zscoreData==Weight_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,7) = (Weight(cur_mouse)-nanmean(Weight))/nanstd(Weight);
    rawData(indx,7) = Weight(cur_mouse);
end

% SRR is in Column 8
testName{7} = 'Self Righting Reflex';
SRR = [SummarySRR(1).Data; SummarySRR(2).Data];  %This has same order as the USV_mice numbers fyi
SRR_mice = [SummarySRR(1).MouseNum; SummarySRR(2).MouseNum];
for cur_mouse = 1:length(SRR_mice)
    indx = find(zscoreData==SRR_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,8) = (SRR(cur_mouse)-nanmean(SRR))/nanstd(SRR); 
    rawData(indx,8) = SRR(cur_mouse);
end


%%%%%%%%%%%%%%%MOTOR%%%%%%%%%%%%%%%%%%%
%Marbles in Column 4
testName{8} = 'Marble Burying';
Marbles = ([SummaryMarble(1).Number'; SummaryMarble(2).Number']);  %Convert to % marble's buried
Marbles_mice = [SummaryMarble(1).Mouse'; SummaryMarble(2).Mouse'];

for cur_mouse = 1:length(Marbles_mice) 
    indx = find(zscoreData==Marbles_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    zscoreData(indx,9) = (Marbles(cur_mouse)-nanmean(Marbles))/nanstd(Marbles);
    rawData(indx,9) = Marbles(cur_mouse);
end

% Grooming Epochs # is Column 5; 
Epoch = [SummaryGroom(1).NumEpoch';SummaryGroom(2).NumEpoch'];   
Epoch_mice = [SummaryGroom(1).Mouse'; SummaryGroom(2).Mouse'];
testName{9} = 'Grooming Epochs';
for cur_mouse = 1:length(Epoch_mice)
    indx = find(zscoreData==Epoch_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    zscoreData(indx,10) = (Epoch(cur_mouse)-nanmean(Epoch))/nanstd(Epoch); 
    rawData(indx,10) = Epoch(cur_mouse);
end

% Mouse Grooming Duration
testName{10} = 'Grooming Duration';
GroomDuration = [SummaryGroom(1).Duration';SummaryGroom(2).Duration'];
GroomDuration_mice = [SummaryGroom(1).Mouse'; SummaryGroom(2).Mouse'];

for cur_mouse = 1:length(GroomDuration_mice)
    indx = find(zscoreData==GroomDuration_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    zscoreData(indx,11) = (GroomDuration(cur_mouse)-nanmean(GroomDuration))/nanstd(GroomDuration);
    rawData(indx,11) = GroomDuration(cur_mouse);
end

% Mouse Number of transitions (total)
testName{11} = 'Chamber Transitions';
TotTransitions = NaN(1,length(SNCombined));
TotTransitions_mice = NaN(1,length(SNCombined));
for m = 1:length(SNCombined) 
    TotTransitions(m) = SNCombined(m).numTransitions+SACombined(m).numTransitions+HABCombined(m).numTransitions; 
    TotTransitions_mice(m) = SNCombined(m).MouseNum;    
end
for cur_mouse = 1:length(TotTransitions_mice) 
    indx = find(zscoreData==TotTransitions_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,12) = (TotTransitions(cur_mouse)-nanmean(TotTransitions))/nanstd(TotTransitions);
    rawData(indx,12) = TotTransitions(cur_mouse);
end


% Mouse Avg Velocity
testName{12} = 'Exploration Velocity';
AvgVelocity = NaN(1,length(SACombined));
AvgVelocity_mice = NaN(1,length(SACombined));
for m = 1:length(HABCombined) 
    AvgVelocity(m) = nanmean([HABCombined(m).velocity]); 
    AvgVelocity_mice(m) = HABCombined(m).MouseNum; 
end
for cur_mouse = 1:length(AvgVelocity_mice) 
    indx = find(zscoreData==AvgVelocity_mice(cur_mouse),1); %it's set this way so that you can randomize the order and still assign the correct value
    %calc zscoreData
    zscoreData(indx,13) = (AvgVelocity(cur_mouse)-nanmean(AvgVelocity))/nanstd(AvgVelocity);
    rawData(indx,13) = AvgVelocity(cur_mouse);
end


%Add group labels 
for i = 1:length(zscoreData)
    zscoreData(i,14) = isVPA(zscoreData(i,1));
    rawData(i,14) = isVPA(zscoreData(i,1));
end

if opts.randshuffle 
    temp = zscoreData(:,14);
    indx = randperm(length(temp))';
    zscoreData(:,14) = temp(indx);
    rawData(:,14) = rawData(:,indx);
end

warning on
%NOTE: Your removing any mice that do not have complete data sets.
if opts.removenan
    warning(spintf('you are removing %d animals with nan in their zscores',sum(any(isnan(zscoreData,2)))));
    zscoreData(any(isnan(zscoreData),2),:) =[];
    rawData(any(isnan(rawData),2),:) =[];
end


end








