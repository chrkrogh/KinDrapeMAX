function Mold = CheckMeshAspect(DT,z)

nElem = size(DT.ConnectivityList,1);
ElemAspect(1:nElem) = 0.0;

% Extract the element vertex coordinates and get z
VertIdx = DT.ConnectivityList;
xVertCoor = reshape(DT.Points(VertIdx',1),[3, nElem]);
yVertCoor = reshape(DT.Points(VertIdx',2),[3, nElem]);
zVertCoor = reshape(F_DT(DT,z,xVertCoor,yVertCoor),[3, nElem]);

% Loop over elements and calculate aspect ratio
for i = 1:nElem
    
ElemVertCoor = [xVertCoor(:,i) yVertCoor(:,i) zVertCoor(:,i)];    
    
% Calculate edge vectors
EdgeVec = diff(ElemVertCoor([1 2 3 1],:),1);
% Get the 2D edge lengths
ElemSize2D = vecnorm(EdgeVec(1:3,1:2),2,2);

% Check element aspect
ElemAspect(i) = max(ElemSize2D)/min(ElemSize2D);

% if ElemAspect > 1e12
%     fprintf(2,'Encountered badly scaled element\n\n')
%     keyboard
%     return
% end

end

[~,MaxElemAspect] = max(ElemAspect);

fprintf('Checking 2D element aspect ratios\n')
fprintf('Max ratio: %.2g, min ratio: %.2g \n\n',...
    max(ElemAspect),min(ElemAspect))

Mold.MaxElemAspect = MaxElemAspect;
Mold.ElemAspect = ElemAspect;
end