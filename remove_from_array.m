function [ result ] = remove_from_array( init_array, removed_elelements )
result = init_array;
result(ismember(result, removed_elelements)) = [];
end

