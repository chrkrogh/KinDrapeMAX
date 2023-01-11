% AddDrapedFiberAnglesToInp
% This script reads-in the FE .inp file defined by MasterInpFilename and 
% the output data from the kinematic draping program stored in the folder
% DrapeDataFolder. For each shell element in the FE model, the closest
% fiber angle calculated by the draping program is found and assigned.
% Lastly a new .inp file with the added fiber angles in the form of
% sections are written to a sub folder of DrapeDataFolder.

clc
clear
close all

addpath('../Functions/')
addpath('../Structural FE models/')

%% Data and settings

% Folder with results of the draping analysis or optimization
DrapeDataFolder = '../Output/2022-December-19_14-52_1layer_baseline/';

% The input file to which sections should be added
MasterInpFilename = 'TestShellv3';

% The common name for the .mat file with the saved workspace
DrapeData = 'Workspace_vars'; 

% The material name used in <MasterInpFilename>.inp
MatName = 'GlassMechanical';

% The (uniform) thickness of each layer (determined by the mold offsets)
MoldLayerOffset = 20e-3; %1e-3; %20e-3;

% Make the layers a different thickness than the mold offset? Can be left
% empty to use MoldLayerOffset as the thickness
OutputLayerThickness = 1e-3;

% Number of nodes in the shell elements of the master model
nShellElemNodes = 4;

% Show the operation, element by element
DebugPlot = false;

WriteNewInp = true;

%% Load variables from draping analysis (optimization)

DrapeDataPath = [DrapeDataFolder DrapeData];

load([DrapeDataPath '.mat'],'Dra_S');

nLayers = length(Dra_S); 

%% Read .inp

tic

% Read .inp and split at line beginning with *Shell Section minus two lines
% and again at *End Part
StartStr = '*Shell Section';
EndStr = '*End Part';
StartOffset = -2;
EndOffset = -1;

[OldInpTop,OldInpBot] = SplitInp(MasterInpFilename,StartStr,EndStr,...
    StartOffset,EndOffset);

% Read the nodes and connectivity from the .inp
[Node, NodeConnect] = ReadInp(MasterInpFilename,nShellElemNodes);

nElem = size(NodeConnect,1);

fprintf('\nFinished reading master .inp...\n')

%% Get fiber orientations from draping analysis

if DebugPlot
    % Load figure, supressing warnings
    warning('off', 'MATLAB:class:LoadInvalidDefaultElement');
    warning('off', 'MATLAB:load:classError')

    openfig([DrapeDataFolder 'FibDevFig.fig']);

    warning('on', 'MATLAB:class:LoadInvalidDefaultElement');
    warning('on', 'MATLAB:load:classError')

end

% Find centroids of all draping elements of all courses in all layers
% Initialize an array with nLayers x nCourses
DrapeElemCentroid = cell(nLayers,max(cellfun(@(x)length(x),Dra_S)));
for LNo = 1:nLayers
    for CNo = 1:length(Dra_S{LNo})
        DrapeElemCentroid{LNo,CNo} = ...
            CentroidOf2DPolygon(Dra_S{LNo}(CNo).P(:,:,1)',...
            Dra_S{LNo}(CNo).P(:,:,2)',Dra_S{LNo}(CNo).P(:,:,3)');
    end
end

% Loop over all finite elements
FEElemFiberAngle(1:nLayers,1:nElem) = 0.0;
for ENo = 1:nElem
    % Get the nodal coordinates of the current element
    ElemNodeCoor = Node(NodeConnect(ENo,:),:);
    % Find the center of the current element
    ElemCenter = CentroidOf2DPolygon(ElemNodeCoor(:,1),...
        ElemNodeCoor(:,2),ElemNodeCoor(:,3))';
    % Find the element normal
    ElementEdge1 = ElemNodeCoor(2,:) - ElemNodeCoor(1,:);
    ElementEdge2 = ElemNodeCoor(3,:) - ElemNodeCoor(1,:);
    ElementNormal = cross(ElementEdge2,ElementEdge1);
    % Normalize the element normal
    ElementNormal = ElementNormal / norm(ElementNormal);
    % Check that z is positive
    if ElementNormal(3) < 0
        ElementNormal = -1 * ElementNormal;
    end
    % Loop over all layers
    for LNo = 1:nLayers
        nCourses_curr = length(Dra_S{LNo});
        % Compute the offset center of the FE in current layer
        ElemOffsetCenter = ElemCenter + ElementNormal * ...
            MoldLayerOffset * (LNo-1);  
        % Find closets center of the drape elements
        Val_min = inf;
        for CNo = 1:nCourses_curr
            [Val,Idx] = min(vecnorm(...
                ElemOffsetCenter'-DrapeElemCentroid{LNo,CNo},2,1));
            if Val < Val_min
                Val_min = Val;
                CourseNo_min = CNo;
                Idx_min = Idx;
            end
        end

        DrapeElemData = squeeze(Dra_S{LNo}(CourseNo_min).P(Idx_min,:,:));
        DrapeElemCentroid_min = DrapeElemCentroid{LNo,CourseNo_min}(:,Idx_min);

        % Get the drape element fiber angle
        FEElemFiberAngle(LNo,ENo) = mean(DrapeElemData(:,5),1);

        if DebugPlot
            % Plotting
            plot3(ElemNodeCoor([1:end, 1],1),ElemNodeCoor([1:end, 1],2),...
                ElemNodeCoor([1:end, 1],3),'m-')
            plot3(ElemCenter(1),ElemCenter(2),ElemCenter(3),'mx')
            plot3(DrapeElemData([1:end, 1],1),DrapeElemData([1:end, 1],2),...
                DrapeElemData([1:end, 1],3),'k-','LineWidth',2)
            plot3(DrapeElemCentroid_min(1),DrapeElemCentroid_min(2),...
                DrapeElemCentroid_min(3),'kx','LineWidth',2)
            % Element normal
            StackHeight = MoldLayerOffset*(nLayers-1);
            ElementNormalEnd = ElemCenter + ElementNormal * StackHeight;
            plot3([ElemCenter(1) ElementNormalEnd(1)],...
                  [ElemCenter(2) ElementNormalEnd(2)],...
                  [ElemCenter(3) ElementNormalEnd(3)],...
                'm--','LineWidth',2)
            ElementOffsetCenters(1,:) = linspace(ElemCenter(1),...
                ElementNormalEnd(1),nLayers);
            ElementOffsetCenters(2,:) = linspace(ElemCenter(2),...
                ElementNormalEnd(2),nLayers);
            ElementOffsetCenters(3,:) = linspace(ElemCenter(3),...
                ElementNormalEnd(3),nLayers);
            plot3(ElementOffsetCenters(1,:),ElementOffsetCenters(2,:),...
                ElementOffsetCenters(3,:),'mo')
          
            pause
        end
    end
    if ENo == round(0.25*nElem)
        fprintf('\nDrilled fiber angles for 25%% of finite elements...\n')
    elseif ENo == round(0.5*nElem)
        fprintf('\nDrilled fiber angles for 50%% of finite elements...\n')
    elseif ENo == round(0.75*nElem)
        fprintf('\nDrilled fiber angles for 75%% of finite elements...\n')
    end
end

fprintf('\nDrilled fiber angles for 100%% of finite elements...\n')

%% Create cell array of strings with element sets

% Create cell array with element set definitions, i.e. one element per set
InpElemSetDef = cell(2*nElem,1);
for i = 1:nElem
    InpElemSetDef{i*2-1,:} = ['*Elset, elset=Elem_' num2str(i) '_set'];
    InpElemSetDef{i*2,:} = ['  ' num2str(i) ','];
end

%% Create cell array of strings with shell sections
if isempty(OutputLayerThickness)
    OutputLayerThickness = MoldLayerOffset;
end

% Create the shell sections corresponding to the element (set)
InpSectionDef = cell(nElem*(nLayers+1),1);
LineCtr = 1;
for ENo = 1:nElem
    ElemSetName = ['Elem_' num2str(ENo) '_set'];
    InpSectionDef{LineCtr} = ...
        ['*Shell Section, elset=' ElemSetName  ', composite, orientation=Ori-1, offset=SNEG'];
    LineCtr = LineCtr + 1;
    for LNo = 1:nLayers
        InpSectionDef{LineCtr} = [num2str(OutputLayerThickness) ', '...
            '3, '...
            MatName ', '...
            num2str(FEElemFiberAngle(LNo,ENo)) ];
        LineCtr = LineCtr + 1;
    end
end

fprintf('\nElement sets and shell sections ready to write...\n')

toc

%% Write the new .inp

if WriteNewInp

    fprintf('\n\nWrite new .inp?\n')

    pause

    OutputFolder = [DrapeDataFolder 'FE_' MasterInpFilename '/'];

    % Check if output folder already exist
    if exist(OutputFolder,'dir')
        fprintf('\nOutput folder already exists. Continue and overwrite?\n')
        pause
    end

    % Create output folder
    mkdir(OutputFolder)

    FileId = fopen([OutputFolder MasterInpFilename '_drapesec.inp'],'w');

    % Write top part
    fprintf(FileId,'%s \r\n',OldInpTop{1,1}{:,1});

    % Write set definitions
    fprintf(FileId,'%s \r\n',InpElemSetDef{:,1});

    % Write composite shell section definitions
    fprintf(FileId,'%s \r\n',InpSectionDef{:,1});

    % Write bottom part
    fprintf(FileId,'%s \r\n',OldInpBot{1,1}{:,1});

    fclose(FileId);

end

%% Aux functions

function C = CentroidOf2DPolygon(x,y,z)
% This function computes the centroid of 2D polygon(s) with vertices given
% by the arrays x, y, and z (vertices x nPolygons)

% First xy-projection
A = x(1:end,:).*y([2:end,1],:)-x([2:end,1],:).*y(1:end,:);
As = sum(A,1)/2;
x_bar = (sum((x([2:end,1],:)+x(1:end,:)).*A,1)*1/6)./As;
y_bar = (sum((y([2:end,1],:)+y(1:end,:)).*A,1)*1/6)./As;

% Then xz-projection
A = x(1:end,:).*z([2:end,1],:)-x([2:end,1],:).*z(1:end,:);
As = sum(A,1)/2;
z_bar = (sum((z([2:end,1],:)+z(1:end,:)).*A,1)*1/6)./As;

% Collect output
C = [x_bar ; y_bar ; z_bar];
end

%% Legacy

% Find the number of elements
% Loop through the top .inp from the last line and backwards and find the
% line where the element defintions end, i.e. consecutive lines with
% nElemenNodes+1 integers. The 1st integer is the element number
% Out = 0;
% for i = length(OldInpTop{1}):-1:1
%     Out_prev = Out;
%     Out = sscanf(OldInpTop{1}{i,1},'%d,');
%     if length(Out) == nElemNodes + 1 && length(Out_prev) == nElemNodes + 1
%         nElem = Out_prev(1);
%         break
%     end
% end


