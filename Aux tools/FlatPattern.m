function [FlatGrid,Boundary,FigHan] = FlatPattern(Dra_S,LNo,CNo)
% This function computes (estimates) the flat pattern of a draped course.
% If the course has been split up in multiple parts due to trimming, only
% one part is located and converted into a flat pattern.
% The first input argument is Dra_S, i.e. a struct with draping information 
% for the stack. The second and third argument are the Layer number and 
% course number, respectively. 

P = Dra_S{LNo}(CNo).P;
Grid = size(Dra_S{LNo}(CNo).Node,[1 2]);
nIniNode_curr = Dra_S{LNo}(CNo).nIniNode;
OrgNode =  1 + nIniNode_curr(3:4);

FlatGrid = NaN([Grid 2]);

Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';

% figure 
% hold on
% axis equal

HasFoundOrigin = false;
% Number of constrained nodes in each quadrant
nConNode = nIniNode_curr([1 2 ; 3 2 ; 3 4 ; 1 4]);
% Loop over quadrant, row and column analogous to Step 2
for i = 1:4
    for j = OrgNode(1) + (0:nConNode(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = OrgNode(2) + (0:nConNode(i,2)-1)*(Dir1(i,2)+Dir2(i,2))
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
            
            % Define the cell vertices
            % Order of vertices from 1: 1,2,3,4
            Vert1 = FlatGrid(Idx(1,1:2));
            Vert2 = Vert1 + Dir2(i,:)*CellEdgeLength(1);
            Vert3 = Vert2 + Dir1(i,:)*CellEdgeLength(2);
            Vert4 = Vert1 + Dir1(i,:)*CellEdgeLength(4);

            % Insert into FlatGrid
            FlatGrid(Idx(:,1:2)) = [Vert1 ; Vert2 ; Vert3 ; Vert4];

            % Incremental plot?
            %plot(FlatGrid(Idx(:,1)),FlatGrid(Idx(:,2)),'o')
        end
    end
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