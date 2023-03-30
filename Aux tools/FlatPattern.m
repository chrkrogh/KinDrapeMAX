function [FlatGrid,Boundary,FigHan] = FlatPattern(Dra_S,LNo,CNo,OutputName)
% BETA
% This function computes (estimates) the flat pattern of a draped course.
% If the course has been split up in multiple parts due to trimming, only
% one part is located and converted into a flat pattern.
% There is a known issue which means that in some cases, the flat pattern
% can not be computed, e.g. in some cases where the origin node is outside 
% of the net boundary and thereby is trimmed away. Maybe a different
% algorithm needs to be considered.
% The first input argument is Dra_S, i.e. a struct with draping information
% for the stack. The second and third argument are the Layer number and
% course number, respectively. The fourth argument is the name used for 
% storing the output files (if ommitted, no output is stored).

% The idea with the present algorithm is to go through all cells in the
% same way they were placed in Step 2 during the draping analysis. Then,
% for each cell. retrieve the updated / trimmed cell edge lengths from the 
% P array and place the nodes in a grid with 90 deg between all cell edges 
% using these edge lengths.

P = Dra_S{LNo}(CNo).P;
Grid = size(Dra_S{LNo}(CNo).Node,[1 2]);
nIniNode_curr = Dra_S{LNo}(CNo).nIniNode;
OrgNode =  1 + nIniNode_curr(3:4);
% Estimate d (discretization) from distance between Nodes in grid
Node = Dra_S{LNo}(CNo).Node;
d = mean(vecnorm(Node(1,2:Grid(2),:)-Node(1,1:Grid(2)-1,:),2,3));

FlatGrid = NaN([Grid 2]);

Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';

% figure
% hold on
% axis equal

HasFoundOrigin = false;
ErrorFlag = false;
% Number of constrained nodes in each quadrant
nConNode = nIniNode_curr([1 2 ; 3 2 ; 3 4 ; 1 4]);
% Loop over quadrant, row and column analogous to Step 2
for i = 1:4
    for j = OrgNode(1) + (0:nConNode(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = OrgNode(2) + (0:nConNode(i,2)-1)*(Dir1(i,2)+Dir2(i,2));
            % Def. single idx
            [Idx,CellNo] = LinearCellIdx(Grid,j,k,Dir1,Dir2,i);

            % Check if the current vertex #1 (OrgNode in 1st iteration) is
            % finite, i.e. not NaN, which means that it is inside the net
            % boundary. If not continue with iterations until a vertex #1
            % is found which is finite
            if HasFoundOrigin == false
                if ~all(isnan(squeeze(P(CellNo,1,1:3))))
                    % Use the current vertex #1 as the zero point
                    FlatGrid(j,k,1:2) = 0.0;
                    % Set the flag to true
                    HasFoundOrigin = true;
                else
                    % If the current vertex #1 is NaN, continue with the
                    % iterations
                    continue
                end
            end

            % Order of cell edges: 1-2, 2-3, 3-4, 1-4s
            CellEdgeLength = vecnorm(squeeze(P(CellNo,[1 2 3 4],1:3) - ...
                P(CellNo,[2 3 4 1],1:3)),2,2);

            % If any of the cell edges are NaN, i.e. the cell is trimmed
            if any(isnan(CellEdgeLength))
                continue
            end

            % If there is a difference in cell lengths in the first column,
            % define adjustment variable that makes sure that the
            % perimeter/boundary vertices gets moved rather than the
            % internal vertices
            % Only relevant if the origin node is on the boundary
            if (i == 2 && k == 1) || (i == 4 && k == OrgNode(2))
                Vert2_Dir2adj = CellEdgeLength(4) - CellEdgeLength(2);
            else
                Vert2_Dir2adj = 0;
            end

            if (i == 1 && k == 1) || (i == 3 && k == OrgNode(2))
                Vert4_Dir2adj = CellEdgeLength(1) - CellEdgeLength(3);
            else
                Vert4_Dir2adj = 0;
            end

            if any(isnan(FlatGrid(Idx(1,1:2))))
                fprintf(2,'\n\nCould not complete flat pattern\n\n')
                ErrorFlag = true;
                break
            end

            % Define the cell vertices
            % Order of vertices from 1: 1,2,3,4
            Vert1 = FlatGrid(Idx(1,1:2));
            Vert2 = Vert1 + Dir2(i,:)*CellEdgeLength(1) + Dir1(i,:)*Vert2_Dir2adj;
            Vert3 = Vert2 + Dir1(i,:)*CellEdgeLength(2);
            Vert4 = Vert1 + Dir1(i,:)*CellEdgeLength(4) + Dir2(i,:)*Vert4_Dir2adj;

            % Insert into FlatGrid
            FlatGrid(Idx(:,1:2)) = [Vert1 ; Vert2 ; Vert3 ; Vert4];

            % Incremental plot?
            %plot(FlatGrid(Idx(:,1)),FlatGrid(Idx(:,2)),'o')
        end
        if ErrorFlag == true
            break
        end
    end
    if ErrorFlag == true
        break
    end
end

if ErrorFlag == true
    FlatGrid = [];
    % Instead, estimate the flat pattern as an undeformed grid with all 
    % complete (i.e. not NaN) cells having the nominal discretization 
    % distance, d. Thus there is no effect of cell trimming
    P_reshape = reshape(mean(Dra_S{LNo}(CNo).P(:,:,1),2,"includenan"),Grid(1)-1,[]);
    FlatGrid = NaN([Grid 2]);
    for i = 1:size(P_reshape,1)
        for j = 1:size(P_reshape,2)
            if ~isnan(P_reshape(i,j))
                for k = 0:1
                    for l = 0:1
                        FlatGrid(i+k,j+l,1:2) = d * [i+k-1,j+l-1];
                    end
                end
            end
        end
    end
    fprintf(2,'Estimating the flat pattern as all cells that are not removed completely\n\n')
end


% Determine boundary
x_vec = reshape(FlatGrid(:,:,1),[],1);
x_vec(isnan(x_vec)) = [];
y_vec = reshape(FlatGrid(:,:,2),[],1);
y_vec(isnan(y_vec)) = [];

% Shrink factor. 0: convex hull, 1: compact boundary
ShrinkFactor = 0.5;
BdIdx = boundary(x_vec,y_vec,ShrinkFactor);
if ~isempty(BdIdx)
    Boundary(:,1) = x_vec(BdIdx);
    Boundary(:,2) = y_vec(BdIdx);
else
    Boundary = [NaN NaN];
end


FigHan = figure;
hold on
plot(x_vec,y_vec,'bo')
plot(Boundary(:,1),Boundary(:,2),'r-')
axis equal

%% Store output?

if nargin > 3 && ~ErrorFlag
    Header = {'x-boundary','y-boundary'};

    ExportCell = [Header ; num2cell(Boundary)];

    writecell(ExportCell,['./Output/' OutputName '_FlatPattern.csv']);

    print(FigHan,['./Output/' OutputName '_FlatPattern'],'-dpng','-r300')

    fprintf('\n\nSuccessfully saved .csv and figure to output folder\n\n')
end
