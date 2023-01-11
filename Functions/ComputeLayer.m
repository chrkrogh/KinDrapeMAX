function [Dra_L, nCourses_trim] = ComputeLayer(Mold,Plt,Inp,Set,F,DT,z)

% Unpack variables from input struct
nCourses = Inp.nCourses;
if strcmpi(Set.InpFileType,'course')
    % Single course
    Grid_L{1} = Inp.Grid;
    d_L{1} = Inp.d;
    Org_L{1} = Inp.Org;
    Ang_L = Inp.Ang;
    OrgNode_L{1} = Inp.OrgNode;
    PreShear_L{1} = Inp.PreShear;
    DraContr_L = Inp.DraContr;
    SteerPts_L = Inp.SteerPts;
    SteerPtsRef = Inp.SteerPtsRef;
    LayerPropagation = 'N/A';
    TrimEdge = 1;
else
    % Layer
    LayerPropagation = Inp.LayerPropagation;
    Grid_L = Inp.Grid_L;
    d_L = Inp.d_L;
    Org_L = Inp.Org_L;
    Ang_L = []; % Not used for layers / stack
    OrgNode_L = Inp.OrgNode_L;
    PreShear_L = Inp.PreShear_L;
    DraContr_L = Inp.DraContr_L;
    SteerPts_L = Inp.SteerPts_L;
    SteerPtsRef_L = Inp.SteerPtsRef_L;
    SteerPtsRef = SteerPtsRef_L;
end

% Initialize
nCourses_trim = Inp.nCourses;
Dra_L = struct('Node',[],'GenCurve',{},'nIniNode',[],'SteerPtsAbs',{},...
    'OrgAdj',[],'P',[],'TrimArea',[]);

for i = 1:nCourses

    if Set.MeshOrCurveWarning
        fprintf(2,'\n|=== Course #%d ===|\n',i);
    end

    % Unpack cell
    d = d_L{i};
    Grid = Grid_L{i};
    Ang = Ang_L;
    Org = Org_L{i};
    PreShear = PreShear_L{i};
    DraContr = DraContr_L;
    SteerPts = SteerPts_L(:,i)';

    % Make sure that the OrgNode is at the course edge (left or right) for
    % courses # 2:end
    if i > 1 && strcmpi(LayerPropagation,'right-to-left')
        OrgNode = [Grid(1) OrgNode_L{i}(2)];
    elseif i > 1 && strcmpi(LayerPropagation,'left-to-right')
        OrgNode = [1 OrgNode_L{i}(2)];
    elseif i == 1
        OrgNode = OrgNode_L{1};
    end

    % Step 1: Generators (geodesic from unfolded tri or steering curve)
    Dra = Step1(d,Grid,Ang,Org,OrgNode,PreShear,F,DT,z,DraContr,...
        SteerPts,SteerPtsRef,Mold,Plt,Set);

    % Step 2: Constrained nodes
    Dra = Step2(Dra,Grid,d,OrgNode,F);

    % Trim ply to mold boundary
    if Set.TrimPlyToNetBd
        [Dra, Node_trim] = TrimPlyToBoundary(Mold,Dra,OrgNode,Grid,F,d);
    else
        Dra.TrimArea = 0;
        Node_trim = Dra.Node;
    end
    
    % Store data for the current course
    Dra_L(i) = Dra;

    % Check if the entire length of the current course is trimmed. If so,
    % stop draping more courses in the layer
    if strcmpi(LayerPropagation,'right-to-left')
        TrimEdge = squeeze(Node_trim(1,:,:));
    elseif strcmpi(LayerPropagation,'left-to-right')
        TrimEdge = squeeze(Node_trim(size(Node_trim,1),:,:));
    end
    if all(isnan(TrimEdge(:)))
        if Set.MeshOrCurveWarning && i < nCourses
            fprintf(2,['Mold coverage reached with course #%d / %d.' ...
                ' Aborting draping of layer\n'],i,nCourses);
        end
        nCourses_trim = i;
        break
    end

    if i <= nCourses-1
        % Get edge to use as the new steering curve
        if strcmpi(LayerPropagation,'right-to-left')
            SteerPtsRef{1} = 'Right-CourseEdge';
            Mold.LeftEdge = squeeze(Dra.Node(1,:,:));

            % Remove NaN's in the previous course if it has been cropped
            if all(isnan(Mold.LeftEdge(:)))
                fprintf(2,'Could not determine left edge of previous course\n');
                keyboard
            elseif any(isnan(Mold.LeftEdge(:)))
                Mold.LeftEdge(isnan(Mold.LeftEdge(:,1)),:) = [];
            end
            % Check the number of points extracted from prev. course edge
            % (there must be at least two for the spline interp. to work)
            if size(Mold.LeftEdge,1) < 2
                fprintf(2,'\nOnly one point identified in left course edge\n')
                nCourses_trim = i;
                break
            end

        elseif strcmpi(LayerPropagation,'left-to-right')
            SteerPtsRef{1} = 'Left-CourseEdge';
            Mold.RightEdge = squeeze(Dra.Node(size(Dra.Node,1),:,:));
            % Remove NaN's in the previous course has been cropped
            if all(isnan(Mold.RightEdge(:)))
                fprintf(2,'Could not determine right edge of previous course\n');
                keyboard
            elseif any(isnan(Mold.RightEdge(:)))
                Mold.RightEdge(isnan(Mold.RightEdge(:,1)),:) = [];
            end

            % Check the number of points extracted from prev. course edge
            % (there must be at least two for the spline interp. to work)
            if size(Mold.RightEdge,1) < 2
                fprintf(2,'\nOnly one point identified in right course edge\n')
                nCourses_trim = i;
                break
            end
        end
    end

end
end