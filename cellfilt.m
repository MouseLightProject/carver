function result = cellfilt(predicate, cell_array)
    is_predicate_true = cellfun(predicate, cell_array) ;
    result = cell_array(is_predicate_true) ;    
end