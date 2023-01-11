function DefineOtherScreen
% This function creates a figure window that the user can move to another
% screen/display and maximize such that that "OuterPosition" property can
% be recorded. In this way, plots can appear maximized on the other screen
% when created with the draping program. The data is stored in the
% Temporary files folder.

% Create the figure and the title with instructions
S.f = figure;
title('Please move this figure window to the other monitor and maximize')

% Create a pushbtton with a callback function
S.pb = uicontrol('style','push',...
                 'units','normalized',...
                 'position',[0.3 0.5 0.4 0.1],...
                 'fontsize',14,...
                 'string','Record this location',...
                 'callback',{@pb_call,S});          
end

function pb_call(varargin)
%Callback function for the pushbutton
% Get the structure
S = varargin{3};  
% Get the outerposition in normalized units
set(gcf, 'Units', 'Normalized')
FigPos = S.f.OuterPosition;
% Write to a text file and close the figure
if ~isdeployed
	writematrix(FigPos,'./Temporary files/ScreenConfig.txt');
else
	writematrix(FigPos,'C:/Drape_Output/Temporary files/ScreenConfig.txt');
end
close(S.f);
end