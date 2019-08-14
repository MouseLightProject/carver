function progress_bar(varargin)
    persistent percent_to_two_decimals_last
    if nargin<2 ,
        % this means reset
        percent_to_two_decimals_last = [] ;
        return
    end
    i = varargin{1} ;
    n = varargin{2} ;
    percent = 100*(i/n) ;    
    percent_to_two_decimals = round(percent*100)/100 ;    
    if ~isequal(percent_to_two_decimals, percent_to_two_decimals_last) ,
        bar = repmat('*', [round(percent/2) 1]) ;
        fprintf('[%-50s]: %6.2f%%\n', bar, percent_to_two_decimals) ;
    end
    percent_to_two_decimals_last = percent_to_two_decimals ;
end
