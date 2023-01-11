function AddLayerSliderToFigure(PlyPltObj)
%% Interactively display layers using a slider
% Modified from slider_figure by Asbjørn M. Olesen
% which was modified from https://stackoverflow.com/a/28256519
% (also consider: https://se.mathworks.com/matlabcentral/fileexchange/29544 )
nLayers = length(PlyPltObj); 

% Get figure and plotting area
hFig = gcf; 
handles.axes1 = gca; 

% Create slider and listener object for smooth visualization
handles.SliderFrame1 = uicontrol('Style','slider', ...
                                'Position',[60 20 400 30], ...
                                'Min',1, ...
                                'Max',nLayers, ...
                                'Value',1, ...
                                'SliderStep',[1/nLayers 2/nLayers], ...
                                'Callback',@XListenerCallBack);

handles.SliderFrame2 = uicontrol('Style','slider', ...
                                'Position',[650 20 400 30], ...
                                'Min',1, ...
                                'Max',nLayers, ...
                                'Value',nLayers, ...
                                'SliderStep',[1/nLayers 2/nLayers], ...
                                'Callback',@XListenerCallBack);

handles.SliderxListener1 = addlistener(handles.SliderFrame1,'Value',...
    'PostSet',@(s,e) XListenerCallBack);
handles.SliderxListener2 = addlistener(handles.SliderFrame2,'Value',...
    'PostSet',@(s,e) XListenerCallBack);

% Create text object to indicate active frame
handles.Text1 = uicontrol('Style','Text','Position',[500 20 60 30],...
    'String','Current layers');
handles.Edit1 = uicontrol('Style','Edit','Position',[550 20 50 30],...
    'String',['1:' num2str(nLayers)]);

% Use setappdata to store the image stack and in callbacks, use getappdata 
% to retrieve it and use it. Check the docs for the calling syntax.
setappdata(hFig,'AllLayers',PlyPltObj);

% IMPORTANT. Update handles structure.
guidata(hFig,handles);

% Listener callback, executed when you drag the slider.
    function XListenerCallBack(~,~)
        % Retrieve handles structure. Used to let MATLAB recognize the
        % edit box, slider and all UI components.
        handles = guidata(hFig);
       
        % Here retrieve PlyPltObj using getappdata.
        PlyPltObj = getappdata(hFig,'AllLayers');
        
        % Get current layers
        CurrentLayer1 = round(get(handles.SliderFrame1,'Value'));
        CurrentLayer2 = round(get(handles.SliderFrame2,'Value'));

        CurrentLayer1 = min(CurrentLayer1,CurrentLayer2);
        CurrentLayer2 = max(CurrentLayer2,CurrentLayer1);

        set(handles.Edit1,'String',[num2str(CurrentLayer1) ':' ...
            num2str(CurrentLayer2)]);
        
        for i = 1:nLayers
            if i >= CurrentLayer1 && i <= CurrentLayer2
                set(PlyPltObj{i},'visible','on')
            else
                set(PlyPltObj{i},'visible','off')
            end
        end
    end
end
