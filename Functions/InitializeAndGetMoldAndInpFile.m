function [Inp,Set,Mold,F,DT,z,FileList,InpFilePath,nCourses,nLayers] = ...
    InitializeAndGetMoldAndInpFile(Inp,Set)
% This function does miscellaneous start-up operations, initializes 
% variables and loads the input file and mold for the draping program.

% Check the value of Set.Mode
Mode_options = {'opt','opt-seq','baseline','analysis','storedopt','error'};
if ~any(strcmpi(Set.Mode,Mode_options))
    fprintf(2,'Incorrect value for Set.Mode. Please choose from:\n\n')
    fprintf(2,[['''' strjoin(Mode_options, ''', ''') ''''] '\n\n'])
    error('Stopped program')
end

% Check the value of Set.InpFileType
InpFileType_options = ...
    {'course', 'layer', 'stack', 'layer-opt','stack-opt'};
if ~any(strcmpi(Set.InpFileType,InpFileType_options))
    fprintf(2,'Incorrect value for Set.InpFileType. Please choose from:\n\n')
    fprintf(2,[['''' strjoin(InpFileType_options, ''', ''') ''''] '\n\n'])
    error('Stopped program')
end

% Check mismatch between Set.Mode and Set.InpFileType
if ~strcmpi(Set.Mode,'Analysis') && ~contains(Set.InpFileType,'opt')
    fprintf(2,'Invalid combination of Set.Mode and Set.InpFileType\n\n')
    error('Stopped program')
end

% Initialize the command window log file
if exist('./Temporary files/log.txt','file') == 2
    diary('off');
    delete('./Temporary files/log.txt')
end
diary('./Temporary files/log.txt');
diary('on');

% Read input file
if isfield(Inp,'FileName') && ~isempty(Inp.FileName)
    InpMode = 'Parse';
else
    InpMode = 'ReadDispParse';
end
[FileList,Inp,InpFilePath] = ReadInputFiles(Inp,InpMode,Set.InpFileType);

% Process the input data
if isfield(Inp,'nLayersUse') && ~isempty(Inp(1).nLayersUse)
    Inp(1).nLayers = Inp(1).nLayersUse;
end

if isfield(Inp,'nCoursesUse') && all(~isempty([Inp.nCoursesUse]))
    [Inp.nCourses] = Inp.nCoursesUse;
end

% Unpack variables from struct
if isfield(Inp(1), 'nCourses')
    nCourses = [Inp.nCourses];
else
    Inp(1).nCourses = 1;
    nCourses = 1;
end

if isfield(Inp(1),'nLayers')
    nLayers = Inp(1).nLayers;
else
    Inp(1).nLayers = 1;
    nLayers = 1;
end

% Initialize mold variables
Mold = cell(1,nLayers);
F = cell(1,nLayers);
DT = cell(1,nLayers);
z = cell(1,nLayers);

% Get the different molds (+ offsets)
for LNo = 1:nLayers

    Inp(1).MoldLayerNo = LNo;

    [Mold{LNo},F{LNo},DT{LNo},z{LNo}] = ReadMoldData(Inp(1),Set);

    % Calculate (approximate) the mold offset from previous layer (used for
    % calculating the in-plane stagger distance in StackNonlCon)
    if LNo == 1
        Mold{LNo}.OffsetFromPrevLayer = 0;
    else
        Mold{LNo}.OffsetFromPrevLayer = min(sqrt(...
            (Mold{LNo}.X(round(end/2),round(end/2))-Mold{LNo-1}.X(:)').^2 + ...
            (Mold{LNo}.Y(round(end/2),round(end/2))-Mold{LNo-1}.Y(:)').^2 + ...
            (Mold{LNo}.Z(round(end/2),round(end/2))-Mold{LNo-1}.Z(:)').^2));

    end
end
% Go back to root folder
cd('../../')

% Check mold boundary definition
if any(cellfun(@(x)isempty(x.MoldEdge),Mold))
    if Set.TrimPlyToNetBd
        Set.TrimPlyToNetBd = false;
        fprintf(2,'\nNo mold ext boundary defined (Mold.MoldEdge).\n')
        fprintf(2,'Disabeling trimming to boundary...\n\n')
    end
end

Inp = rmfield(Inp,'MoldLayerNo');

if strcmpi(Set.Timing,'time')
    tic
elseif strcmpi(Set.Timing,'profile')
    profile on
end

end