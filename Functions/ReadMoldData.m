function [Mold,F,DT,z] = ReadMoldData(Inp,Set)

% Initialize variables
MoldEdge = [];
F = [];
DT = [];
z = [];

if isfield(Inp,'MoldName') && ~isempty(Inp.MoldName)
    Folder = strsplit(Inp.MoldName,'___');
    Folder = Folder{1};
    if Inp.MoldLayerNo == 1
        % For the first layer get path (for later cd)
        FilePath = ['./Molds/' Folder '/'];
    else
        % For remaining layers stay on same path (for later cd)
        % Trick for subfolder molds (return path after first /)
        [~,RemaingPath] = strtok(Folder,'/');
        FilePath = ['.' RemaingPath '/'];
    end

    % Remove first part of path (if sub folders are used)
    MoldName_use = strsplit(Inp.MoldName,'/');
    MoldName_use = MoldName_use{end};

    % Check if file exist 
    if isfile([FilePath MoldName_use '.m'])
        cd(FilePath)
        try
            [X,Y,Z,z,F,DT,MoldEdge] = feval(MoldName_use,Inp);
        catch ME
           checkcode([MoldName_use '.m'])
           fprintf(2,'\nError in mold file\n\n')
           rethrow(ME)
        end
    else
        error('Specified mold file does not exist')
    end
end

% Check net boundary lines
if ~isempty(MoldEdge)
    MoldEdgeSort(1) = issorted(MoldEdge{1}(:,2));
    MoldEdgeSort(2) = issorted(MoldEdge{2}(:,2));
    MoldEdgeSort(3) = issorted(MoldEdge{3}(:,1));
    if any(~MoldEdgeSort)
        fprintf(2,['Please input mold edge lines such that left and right '...
            'are in ascending y and bottom is in ascending x \n\n'])
        keyboard
    end
end

% Check aspect ratio
if Set.CheckAspect
    Mold = CheckMeshAspect(DT,z);
end

Mold.X = X;
Mold.Y = Y;
Mold.Z = Z;
Mold.MoldEdge = MoldEdge;

end