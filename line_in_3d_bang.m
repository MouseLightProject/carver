function result = line_in_3d_bang(ax, xyzs, varargin)
    result = line(ax, 'XData', xyzs(:,1), 'YData', xyzs(:,2), 'ZData', xyzs(:,3), varargin{:}) ;
end
