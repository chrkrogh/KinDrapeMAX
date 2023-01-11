function [Plt, FigHan, AxHan] = DrapePlot(Plt,Set,Dra,Mold,F,DT,z,Org)
% This function plots the output of a draping analysis or optimization. It
% creates the figure window(s), plots the mold and other requested objects.
% The plotting of the courses (patch plot) is done in a separate function
% 'PlyPlot', which is called at the end of this function. Thus 'DrapePlot' 
% is only called for the first course and for the subsequent courses 
% 'PlyPlot' is called.

%% Figure window size
if strcmpi(Plt.FigWindowSize,'Normal')
    ScrOutPos = [0.3495 0.5575 0.3000 0.4275];
elseif strcmpi(Plt.FigWindowSize,'Half')
    ScrOutPos = [0 0 0.42 1];
elseif strcmpi(Plt.FigWindowSize,'Full')
    ScrOutPos = [0 0 1 1];
elseif strcmpi(Plt.FigWindowSize,'Full-other')
    % Obtained from OuterPosition property of get(gcf) for full figure on
    % other monitor (first use set(gcf, 'Units', 'Normalized'))
	if ~isdeployed
		ScreenConfigFile = './Temporary files/ScreenConfig.txt';
	else
		ScreenConfigFile = 'C:/Drape_Output/Temporary files/ScreenConfig.txt';
	end
    if exist(ScreenConfigFile,'file')
        ScrOutPos = readmatrix(ScreenConfigFile);
    else
        ScrOutPos = [0.9958 0.0267 1.0083 0.9800];
    end
end

for k = 1:length(Plt.PlyPlt)
    LegCtr = 1;
    if strcmpi(Plt.PlyPlt{k},'Shear')
        Plt.ShearFig = figure('units','normalized','outerposition',ScrOutPos);
    elseif strcmpi(Plt.PlyPlt{k},'Flat')
        Plt.FlatFig = figure('units','normalized','outerposition',ScrOutPos);
    elseif strcmpi(Plt.PlyPlt{k},'FibDev')
        Plt.FibDevFig = figure('units','normalized','outerposition',ScrOutPos);
    end
    hold on
    
    %% Plot mold
    
    % Mold surface, patches and origin
    surf(Mold.X,Mold.Y,Mold.Z,'EdgeColor','none','FaceColor',[162 181 205]/255,...
        'FaceAlpha',0.5,'FaceLighting','gouraud');
    %surf(X,Y,Z,'EdgeColor','none')
    
    % Triangular mesh
    if Plt.TriMesh
        trimesh(DT.ConnectivityList,DT.Points(:,1),DT.Points(:,2),z,...
            'edgecolor','k','facecolor','none')
    end
    
    % Boundary of triangular mesh
    if Plt.MeshBoundary
        MeshBd = freeBoundary(DT);
        GrObj(LegCtr) = plot3(DT.Points(MeshBd,1),DT.Points(MeshBd,2),...
            F(DT.Points(MeshBd,1),DT.Points(MeshBd,2)),'b-');
        LegEntry{LegCtr} = 'Mold ext. boundary';
        LegCtr = LegCtr + 1;
    end
    
    % Net boundary
    if ~isempty(Mold.MoldEdge) && Plt.MoldNetBd
        GrObj(LegCtr) = plot3(Mold.MoldEdge{1}(:,1),Mold.MoldEdge{1}(:,2),...
            Mold.MoldEdge{1}(:,3),'k--');
        plot3(Mold.MoldEdge{2}(:,1),Mold.MoldEdge{2}(:,2),...
            Mold.MoldEdge{2}(:,3),'k--')
        plot3(Mold.MoldEdge{3}(:,1),Mold.MoldEdge{3}(:,2),...
            Mold.MoldEdge{3}(:,3),'k--')
        
        LegEntry{LegCtr} = 'Mold net boundary';
        LegCtr = LegCtr + 1;
    end
    
    % Plot mold elements above aspect ratio limit
    if Set.CheckAspect
        WarnElem = Mold.ElemAspect > Set.AspectLim;
        trimesh(DT.ConnectivityList(WarnElem,:),DT.Points(:,1),...
            DT.Points(:,2),z,...
            'edgecolor','y','facecolor','none','LineWidth',2)
        if any(WarnElem)
            trimesh(DT.ConnectivityList(Mold.MaxElemAspect,:),DT.Points(:,1),...
                DT.Points(:,2),z,...
                'edgecolor','r','facecolor','none','LineWidth',2)
        end
    end
    
    %% Plot ply and related things
    
    % Origin
    if ~isempty(Org)
        GrObj(LegCtr) = scatter3(Org(1),Org(2),F(Org(1),Org(2)),'kx',...
            'LineWidth',5);
        LegEntry{LegCtr} = 'Origin point';
        LegCtr = LegCtr + 1;
    end
    
    % Generator curves
    if isfield(Dra,'GenCurve') && Plt.GenCurve
        for i = find(~cellfun(@isempty,Dra.GenCurve))
            GrObj(LegCtr) = plot3(Dra.GenCurve{1,i}(:,1),...
                Dra.GenCurve{1,i}(:,2),...
                Dra.GenCurve{1,i}(:,3),'r--','LineWidth',3);
            LegEntry{LegCtr} = 'Generator curve';
            if (i == 1 || i == 3) && ~isempty(Dra.SteerPtsAbs{2})
                SteerPts_curr = Dra.SteerPtsAbs{2};
            elseif (i == 2 || i == 4) && ~isempty(Dra.SteerPtsAbs{1})
                SteerPts_curr = Dra.SteerPtsAbs{1};
            else
                continue
            end
            GrObj(LegCtr+1) = scatter3(SteerPts_curr(1,:),SteerPts_curr(2,:),...
                F(SteerPts_curr(1,:),SteerPts_curr(2,:)),100,'ro');
            LegEntry{LegCtr+1} = 'Steering points';
        end
    end
    LegCtr = length(LegEntry)+1;
    
    axis('equal')
    axis('tight');
    xlabel('x'); ylabel('y'); zlabel('z');
    colormap('jet');
    if strcmpi(Plt.PlyPlt{k},'Shear')
        cb = colorbar('Eastoutside');
        cb.Label.String = 'Shear Angle [deg]';
    elseif strcmpi(Plt.PlyPlt{k},'FibDev')
        cb = colorbar('Eastoutside');
        cb.Label.String = 'Fiber Angle Deviation [deg]';
    end
    view(3)
    
    if Plt.DispLegend
        legend(GrObj,LegEntry,'Location','WestOutside','AutoUpdate','off');
    end
    
    set(gcf,'color','w');
    
    FigHan = gcf;
    AxHan = gca;
    
end
% Patch plot of ply cells
% Get the min/max shear and fib. dev. values for the color limits
Plt.MinShear = min(min(Dra.P(:,:,4)));
Plt.MaxShear = max(max(Dra.P(:,:,4)));

Plt.MinFibDev = min(min(Dra.P(:,:,5)));
Plt.MaxFibDev = max(max(Dra.P(:,:,5)));

Plt.PlyColorNo = 1;

Plt = PlyPlot(Dra,Plt,[],F);

end