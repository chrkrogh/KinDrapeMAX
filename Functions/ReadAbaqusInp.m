function [Node, Connectivity] = ReadAbaqusInp(Filename)

% This function reads nodal coordinates and element connecectivity from
% Abaqus .inp file.

%% Nodes

fileID = fopen([Filename '.inp'],'r');

% Find lines where nodal coordinates of ply start and end
lineIDstart = '*Node';
lineIDend = '*Element';

GoOn = 1;
LineNo = 0;
PassedStartLine = 0;

while GoOn == 1
    
    LineNo = LineNo +1;
    
    line = fgetl(fileID);
    
    if strncmp(line,lineIDstart,length(lineIDstart))
        LineStart = LineNo + 1;
        PassedStartLine = 1;
    end
    
    if PassedStartLine == 1 && strncmp(line,lineIDend,length(lineIDend))
        LineEnd = LineNo;
        GoOn = 0;
    end
    
end

% Read nodal coordinates from .inp-file as formatted data
% Go to beginning of file
frewind(fileID)

% Number of lines to read
Nread = LineEnd - LineStart;

% Read and save in array
Node = textscan(fileID,'%*d %f %f %f',Nread,'Delimiter',','...
    ,'Headerlines',LineStart-1);

Node = cell2mat(Node);

%% Read elements for boundary determination

frewind(fileID)

Connectivity = textscan(fileID,'%*d %d %d %d %d','Delimiter',',',...
    'Headerlines',LineEnd);

Connectivity = cell2mat(Connectivity(:,1:3));

fclose(fileID);

end