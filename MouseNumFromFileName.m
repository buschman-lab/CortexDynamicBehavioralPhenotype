function [mouse, group] = MouseNumFromFileName(file_list)
[~, fn] = cellfun(@(x) fileparts(x), file_list, 'UniformOutput',0); 
mouse = cellfun(@(x) str2double(x(6:9)), fn,'UniformOutput',0);
mouse = [mouse{:}];
group = isVPA(mouse);
end