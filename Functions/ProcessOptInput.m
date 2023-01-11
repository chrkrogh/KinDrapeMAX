function Opt = ProcessOptInput(Inp,Set)
% This function processes the optimization part of the input file. This
% means creating the indices for the design variables with corresponding
% upper and lower bounds

if strcmpi(Set.Mode,'Analysis')
    % If not doing optimization, set the output to empty and return
    Opt = [];
    return
end

if ~isfield(Inp,'nLayers')
    Inp(1).nLayers = 1;
end

% Unpact Inp struct
nCourses = [Inp.nCourses];
nLayers = Inp(1).nLayers;

% Set the des var ctr to 0
DesVarCtr = 0;

% Initialize struct
Opt(nLayers).lb = [];

% Loop over all layers
for LNo = 1:nLayers
    %% Origin node in the longitudinal direction
    if Inp(1).DesVar.OrgNodeLong
        % If active, there is a total of nCourses design varibles for each
        % layer with bounds between 1 and Grid(2) (can be different for 
        % each course)
        Opt(LNo).lb.OrgNodeLong = 1 * ones(1,nCourses(LNo));
        Opt(LNo).ub.OrgNodeLong = cellfun(@(x)x(2),Inp(LNo).Grid_L);
        Opt(LNo).DesVarIdx.OrgNodeLong = DesVarCtr + (1:nCourses(LNo));
        % Update the des var ctr
        DesVarCtr = max(Opt(LNo).DesVarIdx.OrgNodeLong);
    else
        Opt(LNo).lb.OrgNodeLong = [];
        Opt(LNo).ub.OrgNodeLong = [];
        Opt(LNo).DesVarIdx.OrgNodeLong = [];
    end

    %% First course transverse offsets
    if Inp(1).DesVar.FirstTransOffset
        if isa(Inp(1).DesVar.FirstTransOffset,'logical')
            % If input as 'true', create des vars for all the steering pts for
            % course #1 (for current layer)
            nSteerPtsLong1stDesVar = size(Inp(LNo).SteerPts_L{1,1},2);
        elseif isa(Inp(1).DesVar.FirstTransOffset,'numeric') && ...
                length(Inp(1).DesVar.FirstTransOffset) > 1
            % If input as a vector, create des vars for corresponding steer 
            % pts idx for 1st course (for current layer)
            nSteerPtsLong1stDesVar = length(Inp(1).DesVar.FirstTransOffset);
        end
        % Create upper and lower bounds
        Opt(LNo).lb.FirstTransOffset = Inp(1).BndAndCon.FirstTransOffset(1) * ...
            ones(1,sum(nSteerPtsLong1stDesVar));
        Opt(LNo).ub.FirstTransOffset = Inp(1).BndAndCon.FirstTransOffset(2) * ...
            ones(1,sum(nSteerPtsLong1stDesVar));
        % Create array with des var indices (each col is a trans offset)
        Opt(LNo).DesVarIdx.FirstTransOffset(1:nSteerPtsLong1stDesVar,1) ...
            = DesVarCtr + (1:nSteerPtsLong1stDesVar);
        if ~isempty(Opt(LNo).DesVarIdx.FirstTransOffset)
            % Update the des var ctr
            DesVarCtr = max(Opt(LNo).DesVarIdx.FirstTransOffset(:));
        end
    else
        Opt(LNo).lb.FirstTransOffset = [];
        Opt(LNo).ub.FirstTransOffset = [];
        Opt(LNo).DesVarIdx.FirstTransOffset = [];
    end


    %% Transverse course offsets
    if Inp(1).DesVar.TransOffset
        if isa(Inp(1).DesVar.TransOffset,'logical')
            % If input as 'true', create des vars for all the steering pts for
            % course #2 to end (for current layer)
            nSteerPtsLongDesVar = ...
                cellfun(@(x)size(x,2),Inp(LNo).SteerPts_L(1,2:nCourses(LNo)));
        elseif isa(Inp(1).DesVar.TransOffset,'numeric') && ...
                length(Inp(1).DesVar.TransOffset) > 1
            % If input as a vector, create des vars for corresponding steer pts
            % idx for each course (for current layer)
            nSteerPtsLongDesVar = ...
                repmat(length(Inp(1).DesVar.TransOffset),1,nCourses(LNo)-1);
        end
        % Create upper and lower bounds
        Opt(LNo).lb.TransOffset = Inp(1).BndAndCon.TransOffset(1) * ...
            ones(1,sum(nSteerPtsLongDesVar));
        Opt(LNo).ub.TransOffset = Inp(1).BndAndCon.TransOffset(2) * ...
            ones(1,sum(nSteerPtsLongDesVar));
        % Loop over course 2:end and create array with des var indices.
        % Rows: transverse offsets, cols: courses
        Opt(LNo).DesVarIdx.TransOffset = ...
            NaN(max(nSteerPtsLongDesVar),nCourses(LNo)-1);
        for k = 1:nCourses(LNo)-1
            Opt(LNo).DesVarIdx.TransOffset(1:nSteerPtsLongDesVar(k),k) = ...
                DesVarCtr + (1:nSteerPtsLongDesVar(k));
            DesVarCtr = Opt(LNo).DesVarIdx.TransOffset(nSteerPtsLongDesVar(k),k);
        end
        if ~isempty(Opt(LNo).DesVarIdx.TransOffset)
            % Update the des var ctr
            DesVarCtr = max(Opt(LNo).DesVarIdx.TransOffset(:));
        end
    else
        Opt(LNo).lb.TransOffset = [];
        Opt(LNo).ub.TransOffset = [];
        Opt(LNo).DesVarIdx.TransOffset = [];
    end

    %% Course widths
    if Inp(1).DesVar.CourseWidth
        % If active, there is a total of nCourses for each layer design vars
        % with bounds between 1 and the number of available course widths
        nAvailableCourseWidths = length(Inp(1).BndAndCon.AvailableCourseWidths);
        Opt(LNo).lb.CourseWidth = 1 * ones(1,nCourses(LNo));
        Opt(LNo).ub.CourseWidth = nAvailableCourseWidths * ones(1,nCourses(LNo));
        Opt(LNo).DesVarIdx.CourseWidth = DesVarCtr + (1:nCourses(LNo));
        DesVarCtr = max(Opt(LNo).DesVarIdx.CourseWidth(:));
    else
        Opt(LNo).lb.CourseWidth = [];
        Opt(LNo).ub.CourseWidth = [];
        Opt(LNo).DesVarIdx.CourseWidth = [];
    end

end
end