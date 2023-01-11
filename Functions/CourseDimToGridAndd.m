function [Grid,d,Size,L_upd] = CourseDimToGridAndd(Size,d_set,Grid)

L_upd = [];

if isempty(Size)
    Size = d_set .* (Grid-1);
    d = d_set;
elseif isempty(Grid)
    
    % Create design space of grid dim. integers
    Gx_vec = 2:1:300;
    Gy_vec = 2:1:300;
    
    % Calculate the necessary discretizations
    d1 = Size(1)./(Gx_vec-1);
    d2 = Size(2)./(Gy_vec-1);
    
    % Find common discretizations
    %[d_common,idx_Gx,idx_Gy] = intersect(d1,d2);
    
    [d1_idx,d2_idx] = ismembertol(d1,d2);
    
    idx_Gx = find(d1_idx);
    idx_Gy = d2_idx(d1_idx);
    d_common = d1(idx_Gx);
    
    % Find discretization closest to set value
    [~,idx] = min(abs(d_common-d_set));
    
    if ~isempty(d_common)
        d = d_common(idx);
        Grid(1) = Gx_vec(idx_Gx(idx));
        Grid(2) = Gy_vec(idx_Gy(idx));
    else
        % Find closest discretization for the width
        [~,d1_idx] = min(abs(d1-d_set));
        
        % Pick the found discertization and Gx
        d = d1(d1_idx);
        Grid(1) = Gx_vec(d1_idx);
        
        % Calculate Gy with that discertization and round up
        Grid(2) = ceil(Size(2)/d + 1);
        
        % Check the new length
        L_upd = d * (Grid(2)-1);
    end
    
    if ~isempty(L_upd)
        fprintf('\n\nAdjusted L to %g\n\n',L_upd)
    end
    
else
    Size = d_set .* (Grid-1);
    d = d_set;
end
end