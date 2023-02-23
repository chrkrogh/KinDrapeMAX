function [FileList,Inp,InpFilePath] = ReadInputFiles(Inp,Mode,Type)
% This function can read, display and parse the input files pertaining to a
% mold (specified through Inp.MoldName). A mold name can have an extension
% indicated by '___', i.e. 'BladeLike___FromMesh'. It is simply a means of
% specifying different variants of the same mold which use the same input
% files. For this reason, the letters before potential '___' are isolated 
% and used for the mold name here.  
% The input argument Mode determines what the function should do and
% return. The options are:
% - 'Read' : Return the available input files in the variable FileList
% - 'Parse': Read and check the input file specified by Inp.Filename
% - 'ReadDispParse': Read, disp., let user choose in cmd window, and parse
% The input argument 'Type' is a filter for the type of input file to look 
% for, indicated by the first capital letters followed by under score:
% - 'course': C_*
% - 'layer' : L_*
% - 'stack' : S_*
% - 'layer-optimization': LO_*
% - 'stack-optimization': SO_*

% Check if '-optimization is specified' as type
OrginalType = Type;
if contains(Type,'opt','IgnoreCase',true)
    % Split at dash/hyphen and keep only the first part
    temp = strsplit(Type,'-');
    Type = temp{1};
    OptimizationFlag = 1;
else
    OptimizationFlag = 0;
end

if strcmpi(Type,'course')
    Size = [];
    RootCourseOffset = [];
end

if strcmpi(Type,'layer')
    nCoursesUse = [];
    nLayers = 1;
end

if strcmpi(Type,'stack')
    nCoursesUse = [];
    nLayersUse = [];
end

% Run the specified input file
MoldName_use = strsplit(Inp.MoldName,'___');
MoldName_use = MoldName_use{1};

InpFilePath = ['./Input files/' MoldName_use '/'];

if ~isfolder(InpFilePath)
    warning('No input file folder exist for current mold.')
    FileList = [];
    Inp = [];
    return
end

if strcmpi(Mode,'ReadDispParse') || strcmpi(Mode,'Read')
    % Read the available input files in folder
    DirInfo = dir(InpFilePath);

    % Find out which are folders
    dirFlags = [DirInfo.isdir];

    % Get the entries that are not folder names
    FileNames = {DirInfo(~dirFlags).name}';

    % Sort based on requested type:
    % - C_* : course
    % - L_* : layer
    % - S_* : stack
    % - LO_* : layer optimization
    % - SO_* : stack optimization
    if strcmpi(Type,'course') && ~OptimizationFlag
        FileNames = FileNames(startsWith(FileNames,'C_'));
    elseif strcmpi(Type,'layer') && ~OptimizationFlag
        FileNames = FileNames(startsWith(FileNames,'L_'));
    elseif strcmpi(Type,'stack') && ~OptimizationFlag
        FileNames = FileNames(startsWith(FileNames,'S_'));
    elseif strcmpi(Type,'layer') && OptimizationFlag
        FileNames = FileNames(startsWith(FileNames,'LO_'));
    elseif strcmpi(Type,'stack') && OptimizationFlag
        FileNames = FileNames(startsWith(FileNames,'SO_'));
    end

    if strcmpi(Mode,'ReadDispParse')
        % Display available file in command window
        FileList = [];

        fprintf('Found the following %s input files: \n',OrginalType)
        % Print to screen
        for i = 1:length(FileNames)
            fprintf([num2str(i),' ',FileNames{i},'\r\n'])
        end

        while true

            FileNo = input('Enter a number to select >> ','s');

            FileNo = str2double(FileNo);

            if ~isnumeric(FileNo) || ~isscalar(FileNo) || isnan(FileNo)...
                    || FileNo < 1 || FileNo > length(FileNames)
                fprintf(2, 'Invalid input number. Try again...\n\n')
            else
                break
            end
        end

        Inp.FileName = FileNames{FileNo};

        % Read selected file
        run([InpFilePath Inp.FileName])

    elseif strcmpi(Mode,'Read')
        % Return the available input files in FileList
        FileList = FileNames;
    end
end

if strcmpi(Mode,'Parse')
    % Read the chosen input file and parse as function outputs
    FileList = [];
    if isfile([InpFilePath Inp.FileName])
        try
            run([InpFilePath Inp.FileName]);
        catch
            checkcode([InpFilePath Inp.FileName])
            error('MATLAB syntax problem with input file')
        end
    else
        fprintf(2,'\nInputfile for current mold/preform was not found\n');
        error('Stopping execution')
    end
end

% Process and check the output
if strcmpi(Mode,'Parse') || strcmpi(Mode,'ReadDispParse')

    VarError = false;

    % Input value checks
    if ~exist('d','var')
        fprintf(2,'\nNo variable d in current input file\n\n');
        VarError = true;
    end
    if ~exist('Grid','var')% && ~isempty(Grid)
        fprintf(2,'\nNo variable Grid in current input file\n\n');
        VarError = true;
    end
    if ~exist('Org','var')
        fprintf(2,'\nNo variable Org in current input file\n\n');
        VarError = true;
    end
    if ~exist('Ang','var') && strcmpi(Type,'course')
        % Ang is only used for courses
        fprintf(2,'\nNo variable Ang in current input file\n\n');
        VarError = true;
    end
    if ~exist('OrgNode','var') && ~exist('OrgNodePct','var')
        fprintf(2,'\nNo variable OrgNode or OrgNodePct in current input file\n\n');
        VarError = true;
    elseif exist('OrgNode','var') && exist('OrgNodePct','var')
        fprintf(2,'\nSpecify either OrgNode or OrgNodePct in current input file\n\n');
        VarError = true;
    elseif exist('OrgNode','var') && ~exist('OrgNodePct','var')
        OrgNodePct = [];
    elseif ~exist('OrgNode','var') && exist('OrgNodePct','var') 
        OrgNode = [];
    end
    if ~exist('DraContr','var')
        fprintf(2,'\nNo variable DraContr in current input file\n\n');
        VarError = true;
    end
    if ~exist('SteerPtsRef','var')
        fprintf(2,'\nNo variable SteerPtsRef in current input file\n\n');
        VarError = true;
    end
    if ~exist('SteerPts','var')
        fprintf(2,'\nNo variable SteerPts in current input file\n\n');
        VarError = true;
    end

    if ~exist('PreShear','var')
        fprintf(2,'\nNo variable PreShear in current input file\n\n');
        VarError = true;
    end

    % Specific checks relating to layer and stack
    if strcmpi(Type,'layer') || strcmpi(Type,'stack')
        if ~exist('nCourses','var')
            fprintf(2,'\nNo variable nCourses in current input file\n\n');
            VarError = true;
        end
        if ~exist('LayerPropagation','var')
            fprintf(2,'\nNo variable LayerPropagation in current input file\n\n');
            VarError = true;
        end
    end

    % Specific check relating to 'stack'
    if strcmpi(Type,'stack')
        if ~exist('nLayers','var')
            fprintf(2,'\nNo variable nLayers in current input file\n\n');
            VarError = true;
        end
    end

    if strcmpi(Type,'course') && ~VarError

        if strcmpi(DraContr{1},'Geo') && strcmpi(DraContr{2},'Steer')
            fprintf(2,...
                ['Longitudinal generator from geodesic in combination with \n' ...
                'transverse generator from steering curve is not supported\n\n'])
            VarError = true;
        end

        % Check values and data types
        ValTypeError = CheckInpFileValues(d,Grid,Org,Ang,OrgNode,OrgNodePct,...
            SteerPts,PreShear,DraContr,SteerPtsRef);

        if ~ValTypeError
            % Check for Grid vs Size definition
            [Grid,d,Size,~] = CourseDimToGridAndd(Size,d,Grid);

            [OrgNode,OrgNodePct] = OrgNodeFromPct(OrgNode,OrgNodePct,Grid);
        end

    elseif strcmpi(Type,'layer') || strcmpi(Type,'stack')

        % Check values and data types
        ValTypeError = false;

        % OrgNode vs OrgNodePct
        OrgNodeFlag = true;
        OrgNodePctFlag = true;
        if isempty(OrgNode) && ~isempty(OrgNodePct)
            OrgNodeFlag = false;
            OrgNode = cell(size(OrgNodePct));
        elseif ~isempty(OrgNode) && isempty(OrgNodePct)
            OrgNodePctFlag = false;
            OrgNodePct = cell(size(OrgNode));
        end

        SizeError = false;
        % Check size of cells
        for LNo = 1:nLayers
            if length(d(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nd must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if length(Grid(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nGrid must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if length(Org(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nOrg must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if OrgNodeFlag && length(OrgNode(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nOrgNode must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if OrgNodePctFlag && length(OrgNodePct(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nOrgNodePct must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if length(PreShear(LNo,:)) ~= nCourses(LNo)
                fprintf(2,'\nPreShear must have %d cells\n\n',nCourses)
                SizeError = true;
            end
            if (strcmpi(Type,'layer') && ...
                    size(SteerPts(LNo,:),2) ~= nCourses)...
                    || ...
                    (strcmpi(Type,'stack') && ...
                    size(SteerPts{LNo},2) ~= nCourses(LNo))     
                fprintf(2,'\nSteerPts must have 2 x %d cells\n\n',nCourses)
                SizeError = true;
            end
        end

        if SizeError
            ValTypeError = true;
        else
            for LNo = 1:nLayers
                for CNo = 1:nCourses
                    Ang_test = 0;
                    DraContr_test = DraContr;
                    SteerPtsRef_test = SteerPtsRef;
                    if strcmpi(Type,'layer')
                        SteerPts_test = SteerPts(:,CNo);
                    elseif strcmpi(Type,'stack')
                        SteerPts_test = SteerPts{LNo}(:,CNo);
                    end
                    % Check each component of the variables in the inp file
                    ValTypeError_loop = CheckInpFileValues(d{LNo,CNo},...
                        Grid{LNo,CNo},Org{LNo,CNo},Ang_test,...
                        OrgNode{LNo,CNo},OrgNodePct{LNo,CNo},...
                        SteerPts_test,PreShear{LNo,CNo},DraContr_test,...
                        SteerPtsRef_test);
                    if ValTypeError_loop && strcmpi(Type,'layer')
                        fprintf(2,'(in the definition for course #%d) \n\n',CNo)
                        ValTypeError = true;
                    elseif ValTypeError_loop && strcmpi(Type,'stack')
                        fprintf(2,['(in the definition for course #%d' ...
                            ', layer #%d) \n\n'],CNo,LNo)
                        ValTypeError = true;
                    end
                end
            end
        end

        % Make sure both OrgNode and OrgNodePct exist as variables
        if ~exist('OrgNodePct','var')
            OrgNodePct = cell(nLayers,max(nCourses));
        elseif ~exist('OrgNode','var')
            OrgNode = cell(nLayers,max(nCourses));
        end
        for q = 1:nLayers
            for p = 1:nCourses(q)
                [OrgNode{q,p},OrgNodePct{q,p}] = ...
                    OrgNodeFromPct(OrgNode{q,p},...
                    OrgNodePct{q,p},Grid{q,p});
            end
        end

        % Check that the origin node is specified correctly for courses 2
        % to end
        OrgNodePct1st = cellfun(@(x)x(1),cellfun(@(x)[x 0],OrgNodePct,'Uniform',0));
        if strcmpi(LayerPropagation,'right-to-left') && ...
                any(OrgNodePct1st(:,2:end) ~= 100,'all')
            fprintf(2,['In the "%s" layer setting, the origin node of course \n'...
                '2:end must be at the right course edge, i.e. with \n'...
                'OrgPctNode(1) = 100. The input will be corrected \n\n'],...
                LayerPropagation)
        end

        % Check that a steering curve is specified as long. generator
        if strcmpi(DraContr{1},'Geo')
            fprintf(2,'For layers and stacks, the longitudinal/lengthwise generator must be a steering curve\n\n')
        end
    end

    % Root course offset from bottom (hidden/undocumented)
    if ~exist('RootCourseOffset','var')
        if strcmpi(Type,'course')
            RootCourseOffset = [];
        elseif strcmpi(Type,'layer')
            RootCourseOffset = cell(1,nCourses);
        elseif strcmpi(Type,'stack')
            RootCourseOffset = cell(nLayers,max(nCourses));
        end
    end

    if OptimizationFlag
        % Check that the optimization settings are included correctly
        % DesVar struct
        if ~exist('DesVar','var')
            fprintf(2,'\nNo variable DesVar in current input file\n\n');
            VarError = true;
        else
            if ~isfield(DesVar,'OrgNodeLong')
                fprintf(2,'\nNo field DesVar.OrgNodeLong in current input file\n\n');
                VarError = true;
            else
                if ~isa(DesVar.OrgNodeLong,'logical')
                    fprintf(2,'\nDesVar.OrgNodeLong must be logical\n\n');
                    VarError = true;
                end
            end
            if ~isfield(DesVar,'TransOffset')
                fprintf(2,'\nNo field DesVar.TransOffset in current input file\n\n');
                VarError = true;
            else
                if ~(isa(DesVar.TransOffset,'logical') || ...
                        (length(DesVar.TransOffset) > 1 && ...
                        all(DesVar.TransOffset == round(DesVar.TransOffset)))...
                        && all(DesVar.TransOffset > 0))
                    fprintf(2,'\nDesVar.TransOffset must be logical or an array of integers > 0\n\n');
                    VarError = true;
                end
            end
            if ~isfield(DesVar,'FirstTransOffset')
                fprintf(2,'\nNo field DesVar.FirstTransOffset in current input file\n\n');
                VarError = true;
            else
                if ~(isa(DesVar.FirstTransOffset,'logical') || ...
                        (length(DesVar.FirstTransOffset) > 1 && ...
                        all(DesVar.FirstTransOffset == round(DesVar.FirstTransOffset)))...
                        && all(DesVar.FirstTransOffset > 0))
                    fprintf(2,'\nDesVar.FirstTransOffset must be logical or an array of integers > 0\n\n');
                    VarError = true;
                end
            end
            if ~isfield(DesVar,'CourseWidth')
                fprintf(2,'\nNo field DesVar.CourseWidth in current input file\n\n');
                VarError = true;
            else
                if ~isa(DesVar.CourseWidth,'logical')
                    fprintf(2,'\nDesVar.CourseWidth must be logical\n\n');
                    VarError = true;
                end
            end
        end
        % Obj struct
        if ~exist('Obj','var')
            fprintf(2,'\nNo variable Obj in current input file\n\n');
            VarError = true;
        else
            if ~isfield(Obj,'Shear')
                fprintf(2,'\nNo field Obj.Shear in current input file\n\n');
                VarError = true;
            else
                if ~isa(Obj.Shear,'logical')
                    fprintf(2,'\nObj.Shear must be logical\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'AngDev')
                fprintf(2,'\nNo field Obj.AngDev in current input file\n\n');
                VarError = true;
            else
                if ~isa(Obj.AngDev,'logical')
                    fprintf(2,'\nObj.AngDev must be logical\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'TrimArea')
                fprintf(2,'\nNo field Obj.TrimArea in current input file\n\n');
                VarError = true;
            else
                if ~isa(Obj.TrimArea,'logical')
                    fprintf(2,'\nObj.TrimArea must be logical\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'p_Shear')
                fprintf(2,'\nNo field Obj.p_Shear in current input file\n\n');
                VarError = true;
            else
                if ~isscalar(Obj.p_Shear) || ~isnumeric(Obj.p_Shear)
                    fprintf(2,'\nObj.p_Shear must be a numeric scalar\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'p_AngDev')
                fprintf(2,'\nNo field Obj.p_AngDev in current input file\n\n');
                VarError = true;
            else
                if ~isscalar(Obj.p_AngDev) || ~isnumeric(Obj.p_AngDev)
                    fprintf(2,'\nObj.p_AngDev must be a numeric scalar\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'MeanOfLayerPNorms')
                fprintf(2,'\nNo field Obj.MeanOfLayerPNorms in current input file\n\n');
                VarError = true;
            else
                if ~isa(Obj.MeanOfLayerPNorms,'logical')
                    fprintf(2,'\nObj.MeanOfLayerPNorms must be logical\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'TargetFiberAng')
                fprintf(2,'\nNo field Obj.TargetFiberAng in current input file\n\n');
                VarError = true;
            else
                if ~isscalar(Obj.TargetFiberAng) || ~isnumeric(Obj.TargetFiberAng)
                    fprintf(2,'\nObj.TargetFiberAng must be a numeric scalar\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'ShearAng_eq')
                fprintf(2,'\nNo field Obj.ShearAng_eq in current input file\n\n');
                VarError = true;
            else
                if ~isscalar(Obj.ShearAng_eq) || ~isnumeric(Obj.ShearAng_eq)
                    fprintf(2,'\nObj.ShearAng_eq must be a numeric scalar\n\n');
                    VarError = true;
                end
            end
            if ~isfield(Obj,'TrimArea_eq')
                fprintf(2,'\nNo field Obj.TrimArea_eq in current input file\n\n');
                VarError = true;
            else
                if ~isscalar(Obj.TrimArea_eq) || ~isnumeric(Obj.TrimArea_eq)
                    fprintf(2,'\nObj.TrimArea_eq must be a numeric scalar\n\n');
                    VarError = true;
                end
            end
        end
        % Bounds and con
        if ~exist('BndAndCon','var')
            fprintf(2,'\nNo variable BndAndCon in current input file\n\n');
            VarError = true;
        else
            if ~isfield(BndAndCon,'TransOffset')
                fprintf(2,'\nNo field BndAndCon.TransOffset in current input file\n\n');
                VarError = true;
            else
                if numel(BndAndCon.TransOffset)~=2 || BndAndCon.TransOffset(1) >= ...
                        BndAndCon.TransOffset(2) || any(~isnumeric(BndAndCon.TransOffset))
                    fprintf(2,['\nBndAndCon.TransOffset must be a two-component numeric \n'...
                        'vector with the 1st component smaller than the 2nd component\n\n']);
                    VarError = true;
                end
            end
            if ~isfield(BndAndCon,'FirstTransOffset')
                fprintf(2,'\nNo field BndAndCon.FirstTransOffset in current input file\n\n');
                VarError = true;
            else
                if numel(BndAndCon.FirstTransOffset)~=2 || BndAndCon.FirstTransOffset(1) >= ...
                        BndAndCon.FirstTransOffset(2) || any(~isnumeric(BndAndCon.FirstTransOffset))
                    fprintf(2,['\nBndAndCon.FirstTransOffset must be a two-component numeric \n'...
                        'vector with the 1st component smaller than the 2nd component\n\n']);
                    VarError = true;
                end
            end
            if ~isfield(BndAndCon,'AvailableCourseWidths')
                fprintf(2,'\nNo field BndAndCon.AvailableCourseWidths in current input file\n\n');
                VarError = true;
            else
                if any(BndAndCon.AvailableCourseWidths <= 0) || ...
                        any(BndAndCon.AvailableCourseWidths ~= ...
                        round(BndAndCon.AvailableCourseWidths)) ...
                        || isempty(BndAndCon.AvailableCourseWidths)
                    fprintf(2,['\nBndAndCon.AvailableCourseWidths must be a vector with integers '...
                        'larger than 0\n\n']);
                    VarError = true;
                end
            end
            if ~isfield(BndAndCon,'MinStaggerDist')
                fprintf(2,'\nNo field BndAndCon.MinStaggerDist in current input file\n\n');
                VarError = true;
            else
                if ~(isempty(BndAndCon.MinStaggerDist) || ...
                        (isscalar(BndAndCon.MinStaggerDist) && ...
                        isnumeric(BndAndCon.MinStaggerDist) && ...
                        BndAndCon.MinStaggerDist > 0))
                        fprintf(2,['\nBndAndCon.MinStaggerDist must be eiter '...
                        'empty or a numeric scalar larger than 0\n\n']);
                    VarError = true;
                end
            end
        end

        % GASet
        if ~exist('GASet','var')
            fprintf(2,'\nNo variable GASet in current input file\n\n');
            VarError = true;
        else
            if ~isa(GASet,'optim.options.GaOptions')
                fprintf(2,'\nGASet must follow the MATLAB optimoptions format\n\n');
                VarError = true;
            end
        end

    end

    % Throw error without a lot of function stack information
    try
        assert(~(VarError || ValTypeError),'Input file error')
    catch ME
        throwAsCaller(ME)
    end

    % Assign variables to output struct
    if strcmpi(Type,'course')
        Inp.Size = Size;
        Inp.Ang = Ang;
        Inp.RootCourseOffset = RootCourseOffset;
    end
    if strcmpi(Type,'course') || strcmpi(Type,'layer')
        Inp.d = d;
        Inp.Grid = Grid;
        Inp.Org = Org;        
        Inp.OrgNode = OrgNode;
        Inp.OrgNodePct = OrgNodePct;
        Inp.SteerPts = SteerPts;
        Inp.PreShear = PreShear;
        Inp.DraContr = DraContr;
        Inp.SteerPtsRef = SteerPtsRef;
    end
    if strcmpi(Type,'layer')
        Inp.nCourses = nCourses;
        Inp.nCoursesUse = nCoursesUse;
        Inp.LayerPropagation = LayerPropagation;
        Inp.Grid_L = Grid;
        Inp.d_L = d;
        Inp.Org_L = Org;
        Inp.OrgNode_L = OrgNode;
        Inp.OrgNodePct_L = OrgNodePct;
        Inp.PreShear_L = PreShear;
        Inp.DraContr_L = DraContr;
        Inp.SteerPts_L = SteerPts;
        Inp.SteerPtsRef_L = SteerPtsRef;
        Inp.RootCourseOffset_L = RootCourseOffset;
    elseif strcmpi(Type,'stack')
        for k = 1:nLayers
            if k == 1
                Inp(k).nLayers = nLayers;
                Inp(k).nLayersUse = nLayersUse;
            end
            Inp(k).nCourses = nCourses(k);
            Inp(k).nCoursesUse = nCoursesUse(k);
            Inp(k).LayerPropagation = LayerPropagation;
            Inp(k).Grid_L = Grid(k,:);
            Inp(k).d_L = d(k,:);
            Inp(k).Org_L = Org(k,:);
            Inp(k).OrgNode_L = OrgNode(k,:);
            Inp(k).OrgNodePct_L = OrgNodePct(k,:);
            Inp(k).PreShear_L = PreShear(k,:);
            Inp(k).DraContr_L = DraContr;
            Inp(k).SteerPts_L = SteerPts{k};
            Inp(k).SteerPtsRef_L = SteerPtsRef;
            Inp(k).RootCourseOffset_L = RootCourseOffset(k,:);
        end
    end

    if OptimizationFlag
        for k = 1:nLayers
            Inp(k).DesVar = DesVar;
            Inp(k).Obj = Obj;
            Inp(k).BndAndCon = BndAndCon;
            Inp(k).GASet = GASet;
        end
    end
end
end