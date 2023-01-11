function [X,Y,Z,z,F,DT,MoldEdge] = BladeLike___FromMesh_5layer(Inp)
% BladeLike mold defined from surface meshed in Abaqus

FolderExt = '5 Layer';

if ~isfield(Inp,'MoldLayerNo') || Inp(1).MoldLayerNo == 1
    % All layer-variants share the same basis surface
    Filename = 'BladeLikeMoldExt';
    Filename_NetBd = 'BladeLikeMoldExt_NetBd';
elseif Inp(1).MoldLayerNo > 5
    error('Requested layer #%d, but only 5 are available with current mold \nCheck input file',...
        Inp(1).MoldLayerNo)
else
    Filename = ['.\' FolderExt '\BladeLikeMoldExt-Offset' ...
        num2str(Inp(1).MoldLayerNo-1)];
    Filename_NetBd = ['.\' FolderExt '\BladeLikeMoldExt_NetBd-Offset' ...
        num2str(Inp(1).MoldLayerNo-1)];
end

% Get the full surface
MoldEdgeDetec = false;

[X,Y,Z,~,F,DT,z] = ...
    ImportAbaqusMeshAsMold(Filename,MoldEdgeDetec);

% Get the edge from the net boundary
MoldEdgeDetec = true;

[~,~,~,MoldEdge,~,~,~] = ...
    ImportAbaqusMeshAsMold(Filename_NetBd,MoldEdgeDetec);

% Combine the mold edge lines to a net boundary polygon
NetBd(:,1) = [flip(MoldEdge{1}(:,1)) ; ...
    MoldEdge{3}(:,1) ;
    MoldEdge{2}(:,1)];
NetBd(:,2) = [flip(MoldEdge{1}(:,2)) ;
    MoldEdge{3}(:,2) ;
    MoldEdge{2}(:,2)];

% Crop X,Y,Z to net boundary
[In,On] = inpolygon(X,Y,NetBd(:,1),NetBd(:,2));
X(~(In | On)) = NaN;
Y(~(In | On)) = NaN;

end