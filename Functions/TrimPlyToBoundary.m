function [Dra, Node_trim] = TrimPlyToBoundary(Mold,Dra,OrgNode,Grid,F,d)
% This function trims the draped ply to the mold boundary, i.e. net
% boundary. Nodes in Node outside the boundary are set to NaN. In a loop 
% over all cells, different actions are taken for P:
% - 3 or 4 nodes outside: the cell is set to NaN.
% - 2 nodes outside: the outside nodes are moved to the boundary, i.e. the
%   cell is trimmed.
% - 1 node outside: do nothing

% Mold boundary (flipped left, bottom and right assembled)
MoldBd = [flipud(Mold.MoldEdge{1}(:,1:2)) ; ...
    Mold.MoldEdge{3}(2:end-1,1:2) ; ...
    Mold.MoldEdge{2}(:,1:2)];

% Locate the points inside the polygon defined by the boundary points
[in,on] = inpolygon(Dra.Node(:,:,1),Dra.Node(:,:,2),...
    MoldBd(:,1),MoldBd(:,2));
Outside = ~(in | on);
% Create linear indices for Node (offset to all three pages/slices)
OutsideIdx = find(Outside(:)) + [0 1 2]*(Grid(1)*Grid(2));

% Set the nodes outside the boundary to NaN
Node_trim = Dra.Node; 
Node_trim(OutsideIdx(:)) = NaN;

% Initialize the trimmed area variable
TrimArea = 0.0;

% Loop over all cells analogous to Step 2
Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';
nConNode = Dra.nIniNode([1 2 ; 3 2 ; 3 4 ; 1 4]);
for k = 1:4
    for Row = OrgNode(1) + (0:nConNode(k,1)-1)*(Dir1(k,1)+Dir2(k,1))
        for Col = OrgNode(2) + (0:nConNode(k,2)-1)*(Dir1(k,2)+Dir2(k,2))
            % Get the linear cell indices and cell number
            [Idx, CellNo] = LinearCellIdx(Grid,Row,Col,Dir1,Dir2,k);
            % Find the vertex numbers of nodes (x coordinates) that are 
            % NaN, i.e. outside
            OutsideVertNo = find(isnan(Node_trim(Idx(:,1))));
            if any(length(OutsideVertNo) == [3,4])
                % Three or four nodes are outside: set the cell to NaN
                Dra.P(CellNo,:,:) = NaN;
                % Add the full cell area to TrimArea
                TrimArea = TrimArea + d^2;
            elseif length(OutsideVertNo) == 2
                % Two nodes of the cell are outside
                % Move the two nodes so that they are on the the edge
                
                % Diagonal nodes outside check
                if (all(ConnectedVert(OutsideVertNo(1)) == ...
                    flip(ConnectedVert(OutsideVertNo(2))))) || ...
                    (all(ConnectedVert(OutsideVertNo(1)) == ...
                    ConnectedVert(OutsideVertNo(2))))
                    % Do nothing like with the case of one node outside?
                    continue
                end

                % Initialize the new vertex coordinates and loop over the
                % two outside nodes
                NewVert = NaN(2,3);
                for r = 1:2
                    % Create indices for the intersecting cell edge
                    Idx1 = OutsideVertNo(r);
                    Idx2 = ConnectedVert(OutsideVertNo(r),OutsideVertNo);
                    
                    % Construct a curve with the two vertices
                    CellEdge = squeeze(Dra.P(CellNo,[Idx1,Idx2],1:2));
                    
                    % Set the cyclic options for the IntersecOfTwoCurves
                    % functtion
                    Cycl1 = 0;
                    Cycl2 = 1;
                    
                    NewVert(r,1:2) = IntersecOfTwoCurves(CellEdge,...
                        MoldBd,Cycl1,Cycl2);
                end
                
                % Get z coordinates
                NewVert(1:2,3) = F(NewVert(1:2,1),NewVert(1:2,2));
                
                % Compute and add the trimmed cell area to TrimArea
                V1 = squeeze(Dra.P(CellNo,OutsideVertNo(1),1:3))';
                V2 = NewVert(1,1:3);
                V3 = squeeze(Dra.P(CellNo,OutsideVertNo(2),1:3))';
                V4 = NewVert(2,1:3);
                CellTrim = AreaOfQuad(V1,V2,V3,V4);
                TrimArea = TrimArea + CellTrim;
                
                % Update points in P
                Dra.P(CellNo,OutsideVertNo(1),1:3) = NewVert(1,1:3);
                Dra.P(CellNo,OutsideVertNo(2),1:3) = NewVert(2,1:3);

                if norm(NewVert(2,1:3)-squeeze(Dra.P(CellNo,Idx2,1:3))') > 10*d
                    keyboard
                end

            elseif length(OutsideVertNo) == 1
                % One node is outside
                % Do nothing?
                %Dra.P(CellNo,:,4) = 90;
            end
        end
    end
end

Dra.TrimArea = TrimArea;

% scatter3(Dra.Node(OutsideIdx(:,1)),...
%          Dra.Node(OutsideIdx(:,2)),...
%          Dra.Node(OutsideIdx(:,3)),'mo')
end

function Area = AreaOfQuad(V1,V2,V3,V4)
% Compute the area of a quadrilateral based on the four vertices and cross
% products
Area = 0.5 * norm(cross((V3 - V1),(V4 - V2)),2);
end

function ConnVertNo = ConnectedVert(VertNo,ExcVertNo)
% This function returns the two vertices connected to a vertex in a cell.
% If two input argument are given, i.e. also ExcVertNo, the function
% excludes these vertex numbers from the output

if VertNo == 1
    ConnVertNo = [2 4];
elseif VertNo == 2
    ConnVertNo = [1 3];
elseif VertNo == 3
    ConnVertNo = [2 4];
elseif VertNo == 4
    ConnVertNo = [3 1];
end

if nargin > 1
    %ConnVertNo = setdiff(ConnVertNo,ExcVertNo);
    ConnVertNo = MY_setdiff(ConnVertNo,ExcVertNo);
end

end

function Z = MY_setdiff(X,Y)
check = false(1, max(max(X), max(Y)));
check(X) = true;
check(Y) = false;
Z = X(check(X));
end