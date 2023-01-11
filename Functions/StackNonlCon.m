function [C, Ceq, MinDist_cell] = StackNonlCon(x,Mold,Plt,Inp,Set,Opt,F,DT,z)
% This function computes the stack and returns the nonlinear inequality 
% contraints which are the stagger distances 

if isfield(Opt,'SeqLayer') && Opt(1).SeqLayer
    % The input structs only have one row / one cell, i.e. for the 
    % current layer
    [~, OptRes, Dra_L_curr] = LayerObj(x,Mold,Plt,Inp,Set,...
        Opt,F,DT,z);

    nCourses_curr = OptRes.nCourses_trim;

    Dra_L_prev = Opt.Dra_L_prev;
    nCourses_prev = Opt.nCourses_prev;
    MaxLongNodes_prev = Opt.MaxLongNodes_prev;

    nLayers = 1;

    LoopVec = 1;
else
    nLayers = Inp(1).nLayers;
    
    % Initialize
    OptRes = cell(1,nLayers);
    Dra_S = cell(1,nLayers);
    for LNo = 1:nLayers
        if isfield(Set,'OptDebug') && Set.OptDebug
            % In case of errors in the program during optimization, the
            % complete error message is not always printed and the debugger
            % will not stop in the function where the error occured. With this
            % try-catch structure, the vector of design variables can be
            % retrieved and the program can be re-run without optimization.
            try
                [~, OptRes{LNo}, Dra_S{LNo}] = LayerObj(x,Mold{LNo},Plt,Inp(LNo),Set,...
                    Opt(LNo),F{LNo},DT{LNo},z{LNo});
            catch ME
                x
                save('./Temporary files/ErrorWorkSpace');
                fprintf('\n Workspace saved. Rerun the program in error mode\n')
                rethrow(ME)
            end
        else
            % Just run normally
            [~, OptRes{LNo}, Dra_S{LNo}] = LayerObj(x,Mold{LNo},Plt,Inp(LNo),Set,...
                Opt(LNo),F{LNo},DT{LNo},z{LNo});
        end
    end

    % Combine into one struct
    OptRes_S = [OptRes{:}];

    LoopVec = 2:nLayers;
end

%% Check the stagger distance 
% I.e. the distances between ply edges between adjacent layers
% That is, for 2nd layer to last layer check the distances between current
% and previous layer

if strcmpi(Inp(1).LayerPropagation,'left-to-right')
    error(['Layer propagation left-to-right may not be supported. '...
        'Please verify'])
end

% Variable that enables plotting for debugging purposes
DebugDisplay = false;

% Initialize the array with minimum edge distances (nLayers x nEdges)
% NB: each course has two edges, thus nEdges = 2*nCourses
MinDist(1:nLayers,1:2*max([Inp(:).nCourses])) = NaN;

% Loop over layers from 2nd to last
for LNo = LoopVec
    
    if DebugDisplay
        figure
        hold on
        title(['Ply edges in layer #' num2str(LNo)])
    end

    if isfield(Opt,'SeqLayer') && Opt(1).SeqLayer

        % Get max no of initialized longitudinal nodes of prev and curr layer
        MaxLongNodes = max(MaxLongNodes_prev,...
            max(cellfun(@(x) x(2),[Inp.Grid_L])));

    else
        LNo_prev = LNo-1;
        nCourses_prev = OptRes_S(LNo_prev).nCourses_trim;
        LNo_curr = LNo;
        nCourses_curr = OptRes_S(LNo_curr).nCourses_trim;

        % Get max no of initialized longitudinal nodes of prev and curr layer
        MaxLongNodes = max(cellfun(@(x) x(2),[Inp(LNo_prev:LNo_curr).Grid_L]));

        Dra_L_prev = Dra_S{LNo_prev};
        Dra_L_curr = Dra_S{LNo_curr};
    end


    % Initialize arrays for left and right previous edges. The number of
    % rows is set to MaxLongNodes
    LeftEdge_prev = NaN(MaxLongNodes,3,nCourses_prev);
    RightEdge_prev = NaN(MaxLongNodes,3,nCourses_prev);

    % Loop over the number of courses in the previous layer
    for i = 1:nCourses_prev
        % Extract the left and right edges from Node
        LeftEdge_temp = squeeze(Dra_L_prev(i).Node(1,:,1:3));
        RightEdge_temp = squeeze(Dra_L_prev(i).Node(end,:,1:3));

        % Store in array
        LeftEdge_prev(1:size(LeftEdge_temp,1),1:3,i) = LeftEdge_temp;
        RightEdge_prev(1:size(RightEdge_temp,1),1:3,i) = RightEdge_temp;

        if DebugDisplay
            plot3(LeftEdge_temp(:,1),LeftEdge_temp(:,2),...
                LeftEdge_temp(:,3),'k--')
            PltHan1 = plot3(RightEdge_temp(:,1),RightEdge_temp(:,2),...
                RightEdge_temp(:,3),'k--');
        end
    end
    
    % Initialize arrays for left and right current edges. The number of
    % rows is set to MaxLongNodes
    LeftEdge_curr = NaN(MaxLongNodes,3,nCourses_curr);
    RightEdge_curr = NaN(MaxLongNodes,3,nCourses_curr);

    % Loop over the number of courses in the current layer
    for i = 1:nCourses_curr
        % Extract the left and right edges from Node
        LeftEdge_temp = squeeze(Dra_L_curr(i).Node(1,:,1:3));
        RightEdge_temp = squeeze(Dra_L_curr(i).Node(end,:,1:3));

        % Store in array
        LeftEdge_curr(1:size(LeftEdge_temp,1),1:3,i) = LeftEdge_temp;
        RightEdge_curr(1:size(RightEdge_temp,1),1:3,i) = RightEdge_temp;

        if DebugDisplay
            plot3(LeftEdge_temp(:,1),LeftEdge_temp(:,2),...
                LeftEdge_temp(:,3),'r--')
            PltHan2 = plot3(RightEdge_temp(:,1),RightEdge_temp(:,2),...
                RightEdge_temp(:,3),'r--');
        end
    end

    % Loop over all the course edges in the current layer (minus 1st, i.e.
    % right edge of 1st course, and last, i.e. left edge of last course)
    % Make sure the structure is preserved if the layer is cropped before
    % all the initialized courses are draped
    ctr = 0;
    % Loop over all courses
    for j = 1:Inp(LNo).nCourses
        % Loop over right and left course edges
        for k = 1:2
            % Always update the counter to preserve the structure
            ctr = ctr + 1;
            if (j == 1 && k == 1) || (j == nCourses_curr && k == 2)
                % Skip right edge of 1st course and left edge of last
                % course
                continue
            elseif j > nCourses_curr
                % If the current layer has been cropped
                continue
            elseif k == 1
                % Right edge of current course
                RefEdge = RightEdge_curr(:,:,j);
            elseif k == 2
                % Left edge of current course
                RefEdge = LeftEdge_curr(:,:,j);
            end

            % Create 3D arrays with point-to-point distances (euclidean norm)
            % Rows: Ref edge pts. Cols: Prev. edge pts. Page: course # of prev.
            % layer
            % Create array with distances from ref edge to right edges of
            % courses in previous layer
            DistMat_r = sqrt(...
                (RefEdge(:,1)-pagetranspose(RightEdge_prev(:,1,2:end))).^2 + ...
                (RefEdge(:,2)-pagetranspose(RightEdge_prev(:,2,2:end))).^2 + ...
                (RefEdge(:,3)-pagetranspose(RightEdge_prev(:,3,2:end))).^2);

            % Create array with distances from ref edge to left edges of
            % courses in previous layer
            DistMat_l = sqrt(...
                (RefEdge(:,1)-pagetranspose(LeftEdge_prev(:,1,1:end-1))).^2 + ...
                (RefEdge(:,2)-pagetranspose(LeftEdge_prev(:,2,1:end-1))).^2 + ...
                (RefEdge(:,3)-pagetranspose(LeftEdge_prev(:,3,1:end-1))).^2);

            % Concatenate the arrays along the 3rd / page dimension
            DistMat = cat(3,DistMat_r,DistMat_l);

            % If distances cannot be computed, continue (This check was
            % created due to an incident where only the first course in the
            % previous layer could be draped and only partially. Therefore 
            % there were no interior course edges to check and DistMat was
            % empty)
            if isempty(DistMat)
                continue
            end

            % Find the minimum distance for the current ref edge and store
            [MinDist(LNo,ctr),MinIdx] = min(DistMat,[],'all','omitnan','linear');

            if DebugDisplay

                % Find the row,col,page index for the minimum value
                [row,col,pg] = ind2sub(size(DistMat),MinIdx);

                % Calculate the number of pages pertaining to the right edges
                pg_r = size(DistMat_r,3);

                % Plot the right and left edges of the previous layer
                if pg <= pg_r

                    pg_use = pg + 1;

                    plot3(RightEdge_prev(col,1,pg_use),...
                        RightEdge_prev(col,2,pg_use),...
                        RightEdge_prev(col,3,pg_use),'ko')
                    plot3(RefEdge(row,1),RefEdge(row,2),RefEdge(row,3),'ro')
                    plot3([RightEdge_prev(col,1,pg_use) RefEdge(row,1)],...
                        [RightEdge_prev(col,2,pg_use) RefEdge(row,2)],...
                        [RightEdge_prev(col,3,pg_use) RefEdge(row,3)],'b-')

                elseif pg > pg_r

                    pg_use = pg - pg_r;

                    plot3(LeftEdge_prev(col,1,pg_use),...
                        LeftEdge_prev(col,2,pg_use),...
                        LeftEdge_prev(col,3,pg_use),'ko')
                    plot3(RefEdge(row,1),RefEdge(row,2),RefEdge(row,3),'ro')
                    plot3([LeftEdge_prev(col,1,pg_use) RefEdge(row,1)],...
                        [LeftEdge_prev(col,2,pg_use) RefEdge(row,2)],...
                        [LeftEdge_prev(col,3,pg_use) RefEdge(row,3)],'b-')

                end
                legend([PltHan1,PltHan2],...
                    'Edges from previous layer','Edges from current layer')
            end
        end
    end
end

%% Output

% Correct for the distance between layers / mold offsetting
% Here, assume a right-angled triangle where:
% - The 1st cathetus is the sought edge-edge distance
% - The 2nd cathetus is the layer distance / offset
% - The hypotenuse is the calculated direct distance
for LNo = LoopVec
    if isfield(Opt,'SeqLayer') && Opt(1).SeqLayer
        OffsetFromPrevLayer_curr = Mold.OffsetFromPrevLayer;
    else
        OffsetFromPrevLayer_curr = Mold{LNo}.OffsetFromPrevLayer;
    end
   
    MinDist(LNo,:) = real(sqrt(MinDist(LNo,:).^2 - ...
        OffsetFromPrevLayer_curr^2));
end

if isfield(Opt,'SeqLayer') && Opt(1).SeqLayer
    C = MinDist(2:end-1)';
else
    % Put in C, making sure it is not affected by cropped layers
    C = zeros(sum([Inp(2:Inp(1).nLayers).nCourses]*2-2),1);
    EndIdx = 0;
    for LNo = LoopVec
        StartIdx = EndIdx + 1;
        EndIdx = StartIdx + Inp(LNo).nCourses*2-2 -1;
        C(StartIdx:EndIdx) = MinDist(LNo,2:Inp(LNo).nCourses*2-1);
    end
end

% The edge distances of potentially cropped courses are set to inf
C(isnan(C)) = inf;

% Subtract from the minimum stagger distance set in the input file
C = Inp(1).BndAndCon.MinStaggerDist - C;

% No equality constraints
Ceq = [];

MinDist_cell = mat2cell(MinDist,ones(1,nLayers));

if any(~isreal(C))
    keyboard
end

end