function NodeShearMean = AverageShearToGridNodes(Dra,OrgNode,Grid)
% This function takes the output from a course stored in the struct 'Dra'
% along with the 'OrgNode' and the 'Grid' and computes an averaged shear
% value for each node. Recall, that the the cells are computed individually
% and thus each node can have contributions from up to four cells.
% This function was used to check the computed draping pattern against
% Abaqus Composites Modeler, i.e. the nodal shear values.

Node = Dra.Node;
P = Dra.P;

Dir1 = [1 0 -1 0 ; 0 1 0 -1]';
Dir2 = [0 -1 0 1; 1 0 -1 0]';

% Initialize NodeShear as NaN to have the dimensions of the grid and with 
% four slices: one slice for each possible cell vertex contribution to the
% particular node
NodeShear = NaN(Grid(1),Grid(2),4);

% Number of constrained nodes in each quadrant
nConNode = Dra.nIniNode([1 2 ; 3 2 ; 3 4 ; 1 4]);
% Loop over the cells, i.e. in quadrants, rows and columns
for i = 1:4
    for j = OrgNode(1) + (0:nConNode(i,1)-1)*(Dir1(i,1)+Dir2(i,1))
        for k = OrgNode(2) + (0:nConNode(i,2)-1)*(Dir1(i,2)+Dir2(i,2))
            % Def. single idx 
            [Idx,CellNo] = LinearCellIdx(Grid,j,k,Dir1,Dir2,i);
            % Loop over each vertex of the current cell
            for VertNo = 1:4
                % Calculate linear indices for the four slices of NodeShear
                % of the current vertex (node)
                NodeShearIdx = Idx(VertNo,1) + numel(Node(:,:,1))*[0 1 2 3];
                % Find the first of the four entries which is NaN, i.e. not
                % used yet
                NextFreeSlice = find(isnan(NodeShear(NodeShearIdx)),...
                    1,'first');
                % Get the corresponding index
                NextFreeIdx = NodeShearIdx(NextFreeSlice);
                % Store the cell vertex shear angle at that index in
                % Nodeshear
                NodeShear(NextFreeIdx) = P(CellNo,VertNo,4);
            end
        end
    end
end
% Compute the mean across the 3rd dimension
NodeShearMean = mean(NodeShear,3,'omitnan');
end