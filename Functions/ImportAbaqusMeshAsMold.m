function [X,Y,Z,MoldEdge,F,DT,z] = ...
    ImportAbaqusMeshAsMold(Filename,MoldEdgeDetec)
% This function imports an Abaqus FE mesh stored in a .inp file. The mesh
% must consist of triangular elements only, i.e. to work with the Delaunay
% triangularion. 
% For the mold edge detection to work, the mold must be approximately 
% rectangular with four edges and oriented with the y-axis in the 
% longitudinal direction and with the z axis pointing upwards.
% The Abaqus mesh is stored as nodes and their connectivity. This data is
% processed to mold data that can be used with KinDrape: a triangulation
% and an interpolation object is created, matrices for surface plotting are
% created.

MoldPlotRes = 0.020;

% Read the .inp file to get the nodes and their connectivity
[Node, Connectivity] = ReadAbaqusInp(Filename);

% Nodal coordinates
x = Node(:,1);
y = Node(:,2);
z = Node(:,3);

% Create the triangulation
DT = triangulation(double(Connectivity),x,y);

% Create interpolation object
F = scatteredInterpolant(x(:),y(:),z(:),'linear','linear');

% Create matrices for surface plotting
nPts_x = ceil((max(x)-min(x))/MoldPlotRes);
nPts_y = ceil((max(y)-min(y))/MoldPlotRes);
[X,Y] = meshgrid(linspace(min(x),max(x),nPts_x),...
    linspace(min(y),max(y),nPts_y));
Z = F(X,Y);


% Determine the mold edges
if ~MoldEdgeDetec
    MoldEdge = [];
else
    % Set angle tolerance between adjacent element edges on the boundary.
    % If the tolerance is exceeded, the node is a corner node
    AngleTol = 45;

    % Find the boundary nodes
    BdNo = DetermineMeshBoundary(Connectivity,3);
    
    % Sort to make the points in order
    xCenter = mean(x(BdNo));
    yCenter = mean(y(BdNo));
    % Calculate the angles from the center
    Angles = atan2d((y(BdNo)-yCenter) , (x(BdNo)-xCenter));
    % Then sort
    [~, SortIdx] = sort(Angles);
    BdNo_Sort = BdNo(SortIdx);
    
    % Divide into segments based on angle between adjacent segments
    Angle(1:length(BdNo_Sort)) = 0.0;
    BdNo_Sort_Aug = [BdNo_Sort(end) ; BdNo_Sort ; BdNo_Sort(1)];
    for i = 2:length(BdNo_Sort_Aug)-1
        
        Pt_curr = [x(BdNo_Sort_Aug(i)) y(BdNo_Sort_Aug(i)) z(BdNo_Sort_Aug(i))];
        Pt_m1 = [x(BdNo_Sort_Aug(i-1)) y(BdNo_Sort_Aug(i-1)) z(BdNo_Sort_Aug(i-1))];
        Pt_p1 = [x(BdNo_Sort_Aug(i+1)) y(BdNo_Sort_Aug(i+1)) z(BdNo_Sort_Aug(i+1))];
        
        v1 = Pt_p1 - Pt_curr;
        v2 = Pt_m1 - Pt_curr;
        
        v1(3) = 0.0; v2(3) = 0.0;
        
        Angle(i) = atan2d(vecnorm(cross(v1,v2),2,2),dot(v1,v2));
    end
    
    Angle(1) = [];
    Angle(end) = [];
    
    CornerIdx = find(abs(Angle-180)>AngleTol);
    
    if length(CornerIdx) ~= 4
        fprintf(2,'Mold edge detection failed during Abaqus mesh processing\n\n')
        keyboard
    end
    
    Seg{1} = [BdNo_Sort(CornerIdx(4):end) ; BdNo_Sort(1:CornerIdx(1))];
    Seg{2} = BdNo_Sort(CornerIdx(1):CornerIdx(2));
    Seg{3} = BdNo_Sort(CornerIdx(2):CornerIdx(3));
    Seg{4} = BdNo_Sort(CornerIdx(3):CornerIdx(4));
    
    SegCenter = [cellfun(@(c)mean(x(c)),Seg)' cellfun(@(c)mean(y(c)),Seg)'];

    
    [~,LeftSeg] = min(SegCenter(:,1));
    [~,RightSeg] = max(SegCenter(:,1));
    [~,BotSeg] = min(SegCenter(:,2));
    
    LeftEdge = [x(Seg{LeftSeg}) y(Seg{LeftSeg}) z(Seg{LeftSeg})];
    RightEdge = [x(Seg{RightSeg}) y(Seg{RightSeg}) z(Seg{RightSeg})];
    BottomEdge = [x(Seg{BotSeg}) y(Seg{BotSeg}) z(Seg{BotSeg})];
    
    % Check that the left and right edges are in ascending y and that the
    % bottom edge is in ascending x
    LeftEdge = sortrows(LeftEdge,2);
    RightEdge = sortrows(RightEdge,2);
    BottomEdge = sortrows(BottomEdge,1);
    
    % Store in cell array
    MoldEdge{1} = LeftEdge;
    MoldEdge{2} = RightEdge;
    MoldEdge{3} = BottomEdge;
end
end