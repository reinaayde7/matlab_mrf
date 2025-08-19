function [idxs,values, fp] = dict_entry_lookup(dict, r, t1, t2, b0)

    for pair = 1:length(t1)
        t_diff = abs(r(:,1)-t1(pair)) + abs(r(:,2)-t2(pair));
        ii = find(t_diff == min(t_diff));
    
        id = find(r(ii,3) == b0(pair));
        idxs(pair) = ii(id(1)); %retain one solution only
        values(pair) = strcat(string(t1(pair)), '/', string(t2(pair)), '/', string(b0(pair)));
        fp(:,pair) = dict(:,1,idxs(pair));
    end


end