function [pairedID, mouseID, grp] = GetPairedMouseID(file_names)

file_list = GrabFiles(['\w*','refit','\w*'],0,{'Z:\Rodent Data\Wide Field Microscopy\VPA Experiments_Spring2018\VPA_Mesomapping\FitMotifs_deconvolution_Social'});

pairedID = NaN(numel(file_list),1);
mouseID = NaN(numel(file_list),1);
for i = 1:numel(file_list)
    temp = load(file_list{i},'mouseID');
    mouseID(i) = temp.mouseID;
    [~,name] = fileparts(file_list{i});    
    if ~isempty(regexp(name,'mouse1'))
        name = erase(name,sprintf('mouse1_Mouse%d-',mouseID(i)));
        pairedID(i) = str2num(name(1:4));
    else
        name = erase(name,sprintf('mouse2_Mouse'));
        pairedID(i) = str2num(name(1:4));    
    end
end

grp = sum([isVPA(pairedID)', isVPA(mouseID)'],2); %zero = both saline, 1 = overlap, 2= both VPA

end





        
    