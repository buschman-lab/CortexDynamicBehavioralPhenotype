function sex = whatSex(mouse_num,varargin)
%%%%%%%%%%%%%%%%%%%%%
% For Brandy: be sure that you add the path of the folder as opts.BaseDir
% 
% @SYNOPSIS: whatSex is a function that will take a list of mouse numbers 
% and output the sex of those mice. 1 = Female 0 = Male. 
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
% @sex: List of sexes (0 or 1) equivalent length to number of
% input mouse numbers. 
% 0 = Male
% 1 = Female
%
%Example 1: 
%mouse_num = [9021 9022];
%sex = whatSex(mouse_num);
%....
%>> sex = [1 1]
%
%Example 2: 
%mouse_num = [9021];
%sex = whatSex(mouse_num);
%....
%>> sex = [1]
%
%%%%%%%%%%%%%%%%%%%%%

%Define options
opts.TableFile = 'Cohort2_Animals.mat';
opts.LitterFile = 'Litters.xlsx';
opts.MouseNumColumn = 1; %column of mouse numbers on excel sheet
opts.SexColumn = 3; %column on the excel sheet
opts.TableName = 'keyTable';
opts.MotherNumColumn = 'MouseNumber';
opts.DoseColumn = 'Instructions';
opts.VPA_String = {'Give VPA'};
if ispc
    opts.BaseDir = 'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\Spock_Code_Repository\MOUSEGROUPINFO';
else
    opts.BaseDir = '/Volumes/Users/macdo/Desktop/Gould_Brandy';
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

%load litter data from the excel file 
[lkdata] = xlsread([opts.BaseDir filesep opts.LitterFile]);

%Preallocate
sex = nan(1,length(mouse_num));

for cur_mouse = 1:length(mouse_num)
    %find current mouse row
    row = find(lkdata(:,opts.MouseNumColumn) == mouse_num(cur_mouse),1);
    if isempty(row)
        error('Mouse %d is not in excel sheet. Please double check number %d in input list',mouse_num(cur_mouse),cur_mouse)
    end 
    sex(cur_mouse) = lkdata(row,opts.SexColumn);
end

%1 = female, 0 = male
end











