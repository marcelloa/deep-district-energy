function data=import_json_file(file_name)
fid = fopen(file_name); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
data = jsondecode(str);
forceArrayFields = ["buildings","units","thermal_grids", "double_connections","substations","parameters","energy_centres"];
data = forceArraysByFieldName(data, forceArrayFields);

function x = forceArraysByFieldName(x, forceFields)
% forceFields: string array of field names that must be encoded as JSON arrays

    if iscell(x)
        for i = 1:numel(x)
            x{i} = forceArraysByFieldName(x{i}, forceFields);
        end
        return;
    end

    if isstruct(x)
        % if struct array, recurse into each element
        if ~isscalar(x)
            for i = 1:numel(x)
                x(i) = forceArraysByFieldName(x(i), forceFields);
            end
            return;
        end

        fn = fieldnames(x);
        for k = 1:numel(fn)
            f = fn{k};
            v = x.(f);

            % If this field should be an array, but MATLAB collapsed it to scalar struct,
            % wrap it in a cell array so jsonencode outputs [...]
            if any(forceFields == string(f)) && isstruct(v) && isscalar(v)
                x.(f) = {v};
                v = x.(f); % update for recursion
            end

            % Recurse into the value
            x.(f) = forceArraysByFieldName(v, forceFields);
        end
    end
end
end