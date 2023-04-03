function MouseGroup = isVPA(mouse_num,varargin)
%%%%%%%%%%%%%%%%%%%%%
% For Brandy: be sure that you add the path of the folder as opts.BaseDir
% 
% @SYNOPSIS: isVPA is a function that will take a list of mouse numbers 
% and output the group that those mice belong to. 1 = VPA 0 = SAL. 
% For examples how to call and use the function see below. 
%
%
% @INPUTS: 
% @mouse_num (required): 
% Either a single mouse number (e.g. 9021) or a list of mouse numbers (e.g.
% [9021 9022]); 
%
% @varagins (optional):
% List of optional inputs included in define options section below. Primarily are
% directory locations. 
%
% @OUTPUTS:
% @MouseGroup: List of group ID (0 or 1) equivalent length to number of
% input mouse numbers. 
% 0 = SAL
% 1 = VPA
%
%Example 1: 
%mouse_num = [9021 9022];
%MouseGroup = isVPA(mouse_num);
%....
%>> MouseGroup = [2 1]
%
%Example 2: 
%mouse_num = [9021];
%MouseGroup = isVPA(mouse_num);
%....
%>> MouseGroup = [1]
%
%%%%%%%%%%%%%%%%%%%%%

%Define options
opts.TableFile = 'Cohort2_Animals.mat';
opts.LitterFile = 'Litters.xlsx';
opts.MouseNumColumn = 1; %column of mouse numbers on excel sheet
opts.LitterColumn = 2; %column on the excel sheet
opts.TableName = 'keyTable';
opts.MotherNumColumn = 'MouseNumber';
opts.DoseColumn = 'Instructions';
opts.VPA_String = {'Give VPA'};
if ispc
    opts.BaseDir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO';
    %'\\cup\buschman\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO';
    % VR temp change path 11/1/22
else
    opts.BaseDir = '/jukebox/buschman/Rodent Data/Wide Field Microscopy/VPA Experiments_Spring2018/Spock_Code_Repository/MOUSEGROUPINFO';
end
%%
%Process optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin)
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end

%% Load MouseNum --> LitterNum file and find each passed mouse litter number. 

%load litter ket data from the excel file 
[lkdata] = xlsread([opts.BaseDir filesep opts.LitterFile]);

%Preallocate
litter_num = nan(1,length(mouse_num));

for cur_mouse = 1:length(mouse_num)
    %find current mouse row
    row = find(lkdata(:,opts.MouseNumColumn) == mouse_num(cur_mouse),1);
    if isempty(row)
        error('Mouse %d is not in excel sheet. Please double check number %d in input list',mouse_num(cur_mouse),cur_mouse)
    end 
    litter_num(cur_mouse) = lkdata(row,opts.LitterColumn);
end

%% Load file and test whether each passed mouse is in the VPA group

%Load from file
temp = load([opts.BaseDir filesep opts.TableFile], opts.TableName);
keyTable = temp.(opts.TableName);

%Which mice recieved VPA?
mother_vpa = false(size(keyTable, 1), 1);
for cur_str = 1:length(opts.VPA_String)
    mother_vpa = mother_vpa | strcmpi(keyTable.(opts.DoseColumn), opts.VPA_String{cur_str});
end

%Check whether our passed list is in this list
is_vpa = double(ismember(litter_num, keyTable.(opts.MotherNumColumn)(mother_vpa)));

%Any numbers not in our list should be set as NaN
is_vpa(~ismember(litter_num, keyTable.(opts.MotherNumColumn))) = NaN;

%% Parse input mouse numbers into groups 1:4 depending on ket and vpa status
MouseGroup = is_vpa;
end

%%









