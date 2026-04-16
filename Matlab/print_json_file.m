function print_json_file(data,file_name)
str = jsonencode(data, "PrettyPrint", true);
% jsonencode uses 2-space indentation by default. Replace leading two-space groups
% with four spaces for each indentation level.
str = regexprep(str, '(?m)^( +)', '${repmat('' '', 1, length($1)*2)}');
fid = fopen(file_name,"w");
fwrite(fid,str);
fclose(fid);
end