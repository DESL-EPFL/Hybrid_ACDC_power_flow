function index_order = get_index(Bus_list, index_unordered)
    index_order = [];
    for i = 1:length(index_unordered) 
        index_order = [index_order ; find(Bus_list == index_unordered(i))];
    end

end

