function [d]=import_output_file(file_name)
    d=import_json_file(file_name);
    bin_filename=replace(file_name,'.json','.bin');
    fileID=fopen(bin_filename);
    var_list=listStructFields(d);
    for i=1:numel(var_list)
        %fprintf(1,'%s\n',var_list{i});
        if endsWith(var_list{i}, '.format') 
            var_name=extractBetween(var_list{i},1,length(var_list{i})-length('.format'));
            var=eval(var_name{1});
            if (strcmp(var.format,'binary')==1)
                set_binary_data(var_name{1},var.info);
            end    
        end    
    end    
    fclose(fileID);

    function y=get_binary_data(offset,count)
        fseek(fileID,offset,'bof');
        y=fread(fileID,count,d.datatype);
    end

    function set_binary_data(varname,metadata)
        %fprintf(1,'%s\n',varname);
        expression=sprintf('%s=get_binary_data(%i,%i);',varname,metadata(1),metadata(2));
        evalin('caller',expression);
    end

    % function names = listStructFields(s, prefix)
    %     if nargin < 2
    %         prefix = inputname(1);
    %         if isempty(prefix)
    %             prefix = 'root';
    %         end
    %     end
    %     names = {};
    %     if ~isstruct(s)
    %         names = {prefix};
    %         return
    %     end
    %     fields = fieldnames(s);
    %     for i = 1:numel(fields)
    %         fname = fields{i}; 
    %         %fprintf(1,'%s\n',fname);
    %         for k = 1:numel(s)   
    %             if numel(s) > 1
    %                 new_prefix = sprintf('%s(%d).%s', prefix, k, fname);
    %             else
    %                 new_prefix = sprintf('%s.%s', prefix, fname);
    %             end
    %             value = s(k).(fname);
    %             if isstruct(value)
    %                 subnames = listStructFields(value, new_prefix);
    %                 names = [names; subnames]; 
    %             else
    %                 names = [names; {new_prefix}]; 
    %             end
    %         end
    %     end
    % end

    function names = listStructFields(s, prefix)
        if nargin < 2
            prefix = inputname(1);
            if isempty(prefix)
                prefix = 'root';
            end
        end
    
        names = {};
    
        % if not a struct, it's a leaf
        if ~isstruct(s)
            names = {prefix};
            return
        end
    
        fields = fieldnames(s);
        for i = 1:numel(fields)
            fname = fields{i};
            for k = 1:numel(s)
                % build prefix for this field (handles struct arrays)
                if numel(s) > 1
                    new_prefix = sprintf('%s(%d).%s', prefix, k, fname);
                else
                    new_prefix = sprintf('%s.%s', prefix, fname);
                end
    
                value = s(k).(fname);
    
                % Case A: direct struct or struct array -> recurse
                if isstruct(value)
                    subnames = listStructFields(value, new_prefix);
                    names = [names; subnames]; %#ok<AGROW>
    
                % Case B: cell (may contain structs or other types)
                elseif iscell(value)
                    % iterate each cell element
                    for ci = 1:numel(value)
                        elem = value{ci};
                        cell_prefix = sprintf('%s{%d}', new_prefix, ci);
    
                        if isstruct(elem)
                            % if element is struct (or struct array), recurse
                            subnames = listStructFields(elem, cell_prefix);
                            names = [names; subnames]; %#ok<AGROW>
                        else
                            % leaf cell element (non-struct)
                            names = [names; {cell_prefix}]; %#ok<AGROW>
                        end
                    end
    
                % Case C: everything else -> leaf
                else
                    names = [names; {new_prefix}]; %#ok<AGROW>
                end
            end
        end
    end
end