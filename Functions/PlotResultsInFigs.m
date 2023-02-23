function Plt = PlotResultsInFigs(Inp,Set,Plt,Dra_S,Mold,F,DT,z,nLayers,...
    nCourses_trim,OptRes_S)
% This function plots the analyzed or optimized draping model. If more than
% one layer is present, UI sliders to control the visible layers are added

PlyPlt_original = Plt.PlyPlt;
Plt.PlyPlt = [];
for PNo = 1:length(PlyPlt_original)
    Plt.PlyPlt{1} = PlyPlt_original{PNo};
    % Initialize
    PlyPltObj = cell(1,nLayers);
    % Loop over all layers
    for LNo = 1:nLayers
        % Loop over the courses in the current layer
        for CNo = 1:nCourses_trim(LNo)
            if CNo == 1 && LNo == 1 && Plt.Display
                % Plot the mold, points, lines and the draped ply
                Plt = DrapePlot(Plt,Set,Dra_S{LNo}(CNo),Mold{LNo},F{LNo},...
                    DT{LNo},z{LNo},Dra_S{LNo}(CNo).OrgAdj);
            elseif (CNo > 1 || LNo > 1) && Plt.Display
                % Only plot the current ply
                Plt.PlyColorNo = CNo;
                Plt = PlyPlot(Dra_S{LNo}(CNo),Plt,Dra_S{LNo}(CNo).OrgAdj,F{LNo});
            end

            % Store the patch plot objects (plies), scatter objects (org point)
            % and line of color m (nodes, if selected) for the slider figure
            PlyPltObj{LNo} = [findobj('type','patch') ; ...
                findobj('type','scatter') ; ...
                findobj('type','line','color','m') ; ...
                findobj('type','line','color','r')];

            if LNo > 1
                % Remove previously plotted patch objects from array
                PlyPltObj{LNo} = setdiff(PlyPltObj{LNo},...
                    vertcat(PlyPltObj{1:LNo-1}));
            end

        end
    end

    if nLayers > 1
        AddLayerSliderToFigure(PlyPltObj)
    end

    dcm = datacursormode(gcf);
    datacursormode on;
    set(dcm, 'updatefcn',@AdvancedDatacursor)

    set(gca,'clipping','off')

end
Plt.PlyPlt = PlyPlt_original;

% If non considering a single course, plot the stack
if ~strcmpi(Set.InpFileType,'course')
    if isempty(OptRes_S) && ~strcmpi(Set.InpFileType,'course')
        % If no optimization has been run, insert data into OptRes_S to it
        % can be used in Plot2DFlattenedStack
        OptRes_S.CourseWidths = [];
        OptRes_S.FirstTransOffsets = [];
        OptRes_S.TransOffsets = [];
        for LNo = 1:length(nCourses_trim)
            OptRes_S(LNo).nCourses_trim = nCourses_trim(LNo);
        end
    end
    % Plot 2D flattened stack (here for right-to-left draping order)
    Plt.StackFig = Plot2DFlattenedStack(Inp,Mold,Set,OptRes_S);
end

end