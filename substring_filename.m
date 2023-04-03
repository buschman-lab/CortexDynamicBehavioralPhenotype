function file_list = substring_filename(file_list,name,flag)

%flag ==1 = remove the substring. Otherwise remove files with it
temp = cellfun(@(x) ~isempty(regexp(x,name)), file_list,'UniformOutput',1);

file_list = file_list(temp==flag);


end