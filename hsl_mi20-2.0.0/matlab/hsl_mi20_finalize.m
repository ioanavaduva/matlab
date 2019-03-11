function hsl_mi20_finalize(i)
% mi20_finalize removes the GLOBAL variables

global mi20_handle
global MI20STATUS

if nargin == 0
    i = 1;
    destroy_all = 1;
else
    destroy_all = 0;
end
nmax = length(MI20STATUS);

if i > nmax
    error(['The index sent to the function does not correspond ' ...
               'to an instance of hsl_mi20']);
end

fullCells = ~cellfun(@isempty,mi20_handle);

if destroy_all
    for j = 1:nmax
        if fullCells(j), hsl_mi20('destroy', mi20_handle{j}); end
    end
    clearvars -global mi20_handle MI20STATUS;
else
    % Check if mi20_handle{i} exists
    if fullCells(i)
        hsl_mi20('destroy',mi20_handle{i});
        MI20STATUS(i) = 0;
    
        remaining_instances = sum(MI20STATUS);
        if ~remaining_instances
            clearvars -global mi20_handle MI20STATUS;
        end
    else
        error(['The index sent to the function does not correspond ' ...
               'to an instance of hsl_mi20']);
    end
end

