function AnalyzeCursorInfo(cursor_info)
% This function displays the distance between two points in mm chosen in a 
% figure window and which are saved to the workspace as the default name 
% 'cursor_info'

Dist = norm(cursor_info(2).Position-cursor_info(1).Position)*1e3;

fprintf('\n\nThe distance between the points is %g mm \n\n',Dist)

end