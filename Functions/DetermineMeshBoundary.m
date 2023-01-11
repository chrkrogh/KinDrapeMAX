function BdNo = DetermineMeshBoundary(NodeConnect,NoOfElemNodes)
% This function determins the boundary of the mesh in the Abaqus inp
% based on the element connectivity data. The output is a vector of indices
% that will give the boundary nodes when indexed in "Node".
% Mixed meshes, i.e. both tri and quad elements are not supported.

% Initialize
nElem = size(NodeConnect,1);
NodeCtr = 0;
ElemBound(1:nElem*NoOfElemNodes,2) = 0.0;

% Loop over all elements
for i = 1:nElem
    
    % Loop over all nodes in element
    for j = 1:NoOfElemNodes
        
        NodeCtr = NodeCtr + 1;
        
        if j < NoOfElemNodes
            
            % Store vector with nodes from each element boundary
            ElemBound(NodeCtr,1:2) = [NodeConnect(i,j) NodeConnect(i,j+1)];
            
        else
            
            ElemBound(NodeCtr,1:2) = [NodeConnect(i,j) NodeConnect(i,1)];
            
        end
        
    end

end

% Sort row-wise for later comparison
ElemBoundSort = sort(ElemBound,2);

% Define vector with all indices of ElemBoundSort
IdxVec(:,1) = 1:length(ElemBoundSort);

% Find the unique indices of ElemeBoundSort
[~, ia, ~] = unique(ElemBoundSort,'rows');

% Find the duplicate indices as the difference 
DuplicateIdx = setdiff(IdxVec,ia);

% Now find all the rows that are only present once in ElemBoundSort
BoundaryNodes = setdiff(ElemBoundSort,ElemBoundSort(DuplicateIdx,:),'rows');

% Concatenate and get unique node numbers
BdNo = unique(BoundaryNodes(:));

end