function Error = CheckInpFileValues(d,Grid,Org,Ang,OrgNode,OrgNodePct,SteerPts,...
    PreShear,DraContr,SteerPtsRef)
% This function checks that the input file values are consistent in terms
% of datatype and values
Error = false;

%% d
if ~isscalar(d) || ~isnumeric(d) || d <= 0
   fprintf(2,'\nd must be a scalar larger than 0\n\n'); 
   Error = true;
end

%% Grid
if numel(Grid)~=2 || any(Grid ~= round(Grid)) || any(Grid <= 0)
   fprintf(2,['\nGrid must be a two-component vector with integers '...
       'larger than 0\n\n']); 
   Error = true;
end

%% Org
if numel(Org)~=2 
   fprintf(2,'\nOrg must be a two-component vector\n\n'); 
   Error = true;
end

%% Ang
if ~isscalar(Ang) || ~isnumeric(Ang)
   fprintf(2,'\nAng must be a scalar\n\n'); 
   Error = true;
end

%% OrgNode vs OrgNodePct
OrgNodeFlag = true;
OrgNodePctFlag = true;
if isempty(OrgNode) && ~isempty(OrgNodePct)
    OrgNodeFlag = false;
elseif ~isempty(OrgNode) && isempty(OrgNodePct)
    OrgNodePctFlag = false;
end

%% OrgNode
if OrgNodeFlag && numel(OrgNode)~=2 || any(OrgNode ~= round(OrgNode)) ...
        || any(OrgNode <= 0)
   fprintf(2,['\nOrgNode must be a two-component vector with integers ' ...
       'larger than 0\n\n']); 
   Error = true;
end
if OrgNodeFlag && any(OrgNode > Grid)
   fprintf(2,'\nOrgNode must be less than Grid \n\n'); 
   Error = true;
end

%% OrgNodePct
if OrgNodePctFlag && numel(OrgNodePct)~=2 || any(OrgNodePct < 0) || any(OrgNodePct > 100)
   fprintf(2,['\nOrgNodePct must be a two-component vector with values '...
       'between 0 and 100\n\n']); 
   Error = true;
end

%% PreShear
if ~isscalar(PreShear) || ~isnumeric(PreShear)
   fprintf(2,'\nPreShear must be a scalar\n\n'); 
   Error = true;
end

%% DraContr
Values_DraContr = {'Geo','Steer'};

if ~iscell(DraContr) || numel(DraContr) ~= 2
    fprintf(2,'\nDraContr must be a two-element cell array\n\n'); 
    Error = true;
end
if ~(isstring(DraContr{1}) || ischar(DraContr{1})) || ...
        ~any(strcmpi(DraContr{1},Values_DraContr))
    fprintf(2,['\nThe first cell of DraContr must be a string '...
        'or character array with the following values:\n']);
    fprintf(2,'''%s'' \n',Values_DraContr{:})
    fprintf('\n')
    Error = true;
end
if ~(isstring(DraContr{2}) || ischar(DraContr{2})) || ...
        ~any(strcmpi(DraContr{2},Values_DraContr))
    fprintf(2,['\nThe second cell of DraContr must be a string '...
        'or character array with the following values:\n']);
    fprintf(2,'''%s'' \n',Values_DraContr{:})
    fprintf('\n')
    Error = true;
end

%% SteerPts
if ~iscell(SteerPts) || size(SteerPts,1) ~= 2
    fprintf(2,'\nSteerPts must be a cell array with two rows\n\n'); 
    Error = true;
end
if strcmpi(DraContr{1},'Steer') && (isempty(SteerPts{1}) || ...
    (size(SteerPts{1},1) ~= 2 || ~isnumeric(SteerPts{1})))
        fprintf(2,['\nThe first cell row of SteerPts must be a numeric '...
            'array with two rows (x and y)\n\n']);
        Error = true;
end
if strcmpi(DraContr{2},'Steer') && (isempty(SteerPts{2}) || ...
        (size(SteerPts{2},1) ~= 2 || ~isnumeric(SteerPts{2})))
    fprintf(2,['\nThe second cell row of SteerPts must be a numeric '...
        'array with two rows (x and y)\n\n']);
    Error = true;
end

%% SteerPtsRef
Values_SteerPtsRef1 = {'Right','Left','Abs'};
Values_SteerPtsRef2 = {'Bottom' , 'Abs'};

if ~iscell(SteerPtsRef) || numel(SteerPtsRef) ~= 2
    fprintf(2,'\nSteerPtsRef must be a two-element cell array\n\n'); 
    Error = true;
end
if ~isempty(SteerPtsRef{1})
    if ~(isstring(SteerPtsRef{1}) || ischar(SteerPtsRef{1})) || ...
            ~any(strcmpi(SteerPtsRef{1},Values_SteerPtsRef1))
        fprintf(2,['\nThe first cell of SteerPtsRef '...
            '(longitudinal direction) must be a string '...
            'or character array with the following values:\n']);
        fprintf(2,'''%s'' \n',Values_SteerPtsRef1{:})
        fprintf('\n')
        Error = true;
    end
end
if ~isempty(SteerPtsRef{2})
    if ~(isstring(SteerPtsRef{2}) || ischar(SteerPtsRef{2})) || ...
            ~any(strcmpi(SteerPtsRef{2},Values_SteerPtsRef2))
        fprintf(2,['\nThe second cell of SteerPtsRef '...
            '(transverse direction) must be a string '...
            'or character array with the following values:\n']);
        fprintf(2,'''%s'' \n',Values_SteerPtsRef2{:})
        fprintf('\n')
        Error = true;
    end
end
end